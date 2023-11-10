import click
from Bio import SeqIO
import logging
import gzip
import json
import options as gopt

@click.command()
@click.option('--prefix', required=False, default="Sample", help='Sample name in first column')
@click.option('--deeparg', required=True, help='deepARG *.ARG file')
@click.option('--lfile', required=False, help='File containing the length of the ARGs')
@click.option('--ncounts', default=0, help='to normalize to this counts')
@click.option('--nlength', default=1432, help='Length of the reference gene (default 16S rRNA: 1432)')
@click.option('--output-file', required=True, help='Output file to write the results (*.type, *.subtype)')
@click.option('--rlength', default=100, help='Length of the reads')
@click.option('--gfilter', default='', help='select an specific antibiotic resistance category')

def deeparg_abundance(prefix, deeparg, lfile, ncounts, nlength, rlength, output_file, gfilter):
    '''

    Perform relative abundance for deepARG.

    This script computes the relative abundance of ARGs respect a reference gene. Normally its using the 16S rRNA.

    ABN[ARG] = ( ARG-like-reads*rlength/ARG-gene-length )/( 16S-reads*rlength/16S-gene-length )

    '''

    logging.basicConfig(
        filename=output_file + '.log',
        level=logging.DEBUG,
        filemode="w",
        format="%(levelname)s %(asctime)s - %(message)s"
    )
    log = logging.getLogger()

    # lfile = gopt.path+'/database/v2/features.gene.length'

    # arg_size = {i.split()[0].split("|")[-1].upper(): i.split() for i in open(lfile)}
    arg_size = {i.split()[0].split("|")[2].upper(): i.split() for i in open(lfile)}
    arg_size.update({i.split()[0].split("|")[-1].upper(): i.split() for i in open(lfile)})

    arg_counts = {}
    arg_category = {}
    log.info('counting args in file')
    for i in open(deeparg):
        if "#" in i[0]: continue
        i = i.split()
        try:
            arg_counts[ i[0] ] += 1
        except Exception as e:
            arg_counts[i[0]] = 1

        try:
            arg_category[i[4]].append(i[0])
            arg_category[i[4]] = list(set(arg_category[i[4]]))
        except Exception as e:
            arg_category[ i[4] ] = [ i[0] ]

    log.debug(arg_category)

    fo1 = open(output_file+'.subtype','w')
    fo1.write("#ARG-group\tReadCount\t16s-NormalizedReadCount\n")

    rel_abundances = {}
    abs_abundances = {}
    for arg in arg_counts:
        arg_like_reads = arg_counts[arg]
        arg_length = int(arg_size[arg][1])

        rel_abundance = (float(arg_like_reads)*rlength/arg_length)/(float(ncounts)*rlength/float(nlength))
        # log.debug( (arg, arg_counts[arg], rel_abundance) )

        fo1.write("\t".join([arg, str(arg_counts[arg]), str(rel_abundance)]) + "\n")

        rel_abundances[arg] = rel_abundance
        abs_abundances[arg] = arg_like_reads

    fo = open(output_file+'.type','w')
    fo.write("#ARG-group\tReadCount\t16s-NormalizedReadCount\n")
    fo1 = open(output_file + '.table.tsv', 'w')
    # fo1.write("#prefix\tantibiotic_class\tARG_name\tabsolute_abundance\trelative_abundance\n")

    tree = []

    for antibiotic in arg_category:
        antibiotic_abn = 0
        abs_abundance = 0
        for arg in arg_category[antibiotic]:
            log.debug((antibiotic, arg, rel_abundances[arg]))
            antibiotic_abn += rel_abundances[arg]
            abs_abundance += abs_abundances[arg]
            tree.append({'antibiotic_class': antibiotic, 'ARG': arg, 'absolute_abundance': abs_abundances[arg], 'relative_abundance': rel_abundances[arg]})
            fo1.write("\t".join([prefix, antibiotic, arg, str(abs_abundances[arg]), str(rel_abundances[arg])])+'\n' )

        log.debug(('SUM:', antibiotic, antibiotic_abn))
        fo.write( "\t".join([str(antibiotic), str(abs_abundance), str(antibiotic_abn)])+'\n' )

    json.dump( tree, open( output_file+'.json', 'w' ) )
    fo1.close()

    if gfilter:
        for selected_type in gfilter.split(','):
            fo = open(output_file + '.' + selected_type + '.tsv', 'w')
            fo.write("#ARG-group\tReadCount\t16s-NormalizedReadCount\n")
            try:
                for arg in arg_category[selected_type]:
                    fo.write( "\t".join([arg, '-' , str(rel_abundances[arg]) ])+'\n' )
            except Exception as e:
                print('Error! \nPlease select one of the following list of antibiotics')
                print("\n".join(arg_category.keys()))

if __name__ == '__main__':
    deeparg_abundance()