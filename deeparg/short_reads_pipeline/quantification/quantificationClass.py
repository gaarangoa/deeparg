from deeparg.short_reads_pipeline.tools.deepargClass import dsize
import sys
import os


def merge(inputFile, path_to_data):
    try:
        cmd = "sort -k1,1 -k2,2n "+inputFile + "  | bedtools merge -c 12,5 -o sum,distinct >"+inputFile+".merged"
        print(cmd)

        x = os.popen(cmd).read()

        genes = {}
        for i in open(inputFile+".merged"):
            subtype, start, end, count, Type = i.split()
            start = int(start)
            end = int(end)
            count = int(count)
            try:
                genes[subtype]['count'] += count
                genes[subtype]['length'] += abs(end-start)
                genes[subtype]['type'].append(Type)
                genes[subtype]['type'] = list(set(genes[subtype]['type']))
            except:
                genes[subtype] = {
                    'count': count,
                    'length': abs(end-start),
                    'type': [Type]
                }
        # print genes
        fo = open(inputFile+".merged.quant", 'w')
        gene_len = dsize(path_to_data)
        # print gene_len
        for gene in genes:
            # print gene, gene_len[gene]
            cov = genes[gene]['length']/float(gene_len[gene][1])
            fo.write("\t".join([
                gene,
                "/".join(genes[gene]['type']),
                str(genes[gene]['count']),
                str(genes[gene]['length']),
                gene_len[gene][1],
                str(round(cov, 3))
            ])+"\n")
        fo.close()
        return True
    except Exception as inst:
        print(str(inst))
        return False
