def normalize(fi16s, fiArg, covi, parameters={}):
    N16s = sum([int(i.split()[1]) for i in open(fi16s) if float(i.split()[2]) >= 100])
    print( "Total number of 16S Reads in the sample: {}".format(N16s) )
    L16s = 1432
    Atype = {}
    Asubtype = {}

    # print N16s
    fo1 = open(fiArg+'.subtype', 'w')
    fo1.write("#ARG-group\tReadCount\t16s-NormalizedReadCount\n")
    for i in open(fiArg):
        subtype, gtype, count, algLen, geneLen, cov = i.split()
        if float(cov) <= covi:
            continue
        if N16s > 0:
            Asubtype[subtype] = (int(count)/float(geneLen))/(float(N16s)/L16s)
        else:
            Asubtype[subtype] = 0

        fo1.write("\t".join([
            subtype,
            str(count),
            str(Asubtype[subtype])
        ])+"\n")
        try:
            Atype[gtype][0] += int(count)
        except:
            Atype[gtype] = [int(count), geneLen]

    fo2 = open(fiArg+'.type', 'w')
    fo2.write("#ARG-category\tReadCount\t16s-NormalizedReadCount\n")
    Xtype = {}
    for itype in Atype:
        if N16s > 0:
            Xtype[itype] = (Atype[itype][0]/float(Atype[itype][1])) / \
                (float(N16s)/L16s)
        else:
            Xtype[itype] = 0
        fo2.write("\t".join([
            itype,
            str(Atype[itype][0]),
            str(Xtype[itype])
        ])+"\n")

    return True
