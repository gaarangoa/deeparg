import os


def dsize(ggdata):
    return {i.split()[0]: i.split() for i in open(ggdata+".len")}


class d16sPipe():
    def __init__(self):
        self.info = ""

    def run(self, fi, ggdata="/data/gg13/"):
        self.green = ggdata
        try:
            os.system(" ".join([
                "bowtie2 -f",
                "--fast-local",
                "--no-unal",
                "-x", self.green,
                "-U", fi,
                "-S", fi+'.sam'
            ]))

            os.system(" ".join([
                "samtools view -bS",
                fi+".sam",
                ">",
                fi+".bam"
            ]))

            os.system(" ".join([
                "samtools sort",
                fi+".bam",
                "-o", fi+".sorted.bam"
            ]))

            os.system(" ".join([
                "bedtools merge",
                "-i", fi+".sorted.bam",
                "-c 1 -o count",
                ">", fi+".sorted.bam.merged"
            ]))

            genes = {}
            for i in open(fi+".sorted.bam.merged"):
                subtype, start, end, count = i.split()
                start = int(start)
                end = int(end)
                count = int(count)
                try:
                    genes[subtype]['count'] += count
                    genes[subtype]['length'] += abs(end-start)
                except:
                    genes[subtype] = {
                        'count': count,
                        'length': abs(end-start),
                    }

            fo = open(fi+".sorted.bam.merged.quant", 'w')
            gene_len = dsize(self.green)
            # print gene_len
            for gene in genes:
                # print gene, gene_len[gene]
                cov = genes[gene]['length']/float(gene_len[gene][1])
                fo.write("\t".join([
                    gene,
                    str(genes[gene]['count']),
                    str(genes[gene]['length']),
                    gene_len[gene][1],
                    str(round(cov, 3))
                ])+"\n")
            fo.close()
            return True
        except:
            return False
