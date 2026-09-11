#!/usr/bin/env python3
# parse gene-panel files to get an updated oncoplex-genes.txt file

import re

inFiles = ['../opx-v5-gene-panel.txt','../opx-v6-gene-panel.txt','../opx-v7-gene-panel.txt','../opx-v8a-gene-panel.txt']
outFile1 = '../oncoplex-genes.txt'

with open(outFile1,'w') as out1:
    # write yo output
    outHeader = ['gene','v5','v6','v7','v8']
    out1.write('\t'.join(outHeader) + '\n')

    lookup = dict()
    for inFile in inFiles:
        with open(inFile) as in1:
            parse1 = in1.readline()
            parse1 = parse1.strip()
            version = re.sub('.*OPX','',parse1)
            version = re.sub('a','',version)
            
            in1.readline()
            parse1 = in1.readline()
            parse1 = parse1.strip().split('\t')

            for gene in parse1:
                if gene == 'gene_list:':
                    continue
                
                if gene not in lookup:
                    lookup[gene] = version
                else:
                    lookup[gene] += ',' + version

    # write for each gene
    genes = sorted(list(lookup.keys()))
    for gene in genes:
        versions = lookup[gene]

        v5 = 'TRUE' if 'v5' in versions else 'FALSE' 
        v6 = 'TRUE' if 'v6' in versions else 'FALSE'
        v7 = 'TRUE' if 'v7' in versions else 'FALSE'
        v8 = 'TRUE' if 'v8' in versions else 'FALSE'            

        # write yo data
        lineOut = [gene,v5,v6,v7,v8]
        out1.write('\t'.join(lineOut) + '\n')

