#!/usr/bin/env python3
# generate a list of oncoplex genes

import re

inFile1 = '../input/oncoplexGeneWebCopy.txt'
outFile1 = '../oncoplex-genes.txt'

excludes = ['microsatellite','multigene','next-generation','tumor','somatic','Gynecological','precision','total','paired','BRCA1&2','GYNPTH']
excludes = dict(zip(excludes,excludes))

with open(inFile1) as in1, open(outFile1,'w') as out1:
    # read ...
    parse1 = in1.readlines()[0]

    # filter...
    parse2 = parse1.split(',')
    parse3 = [re.sub('^[ ]','',field) for field in parse2]
    parse3 = [re.sub('[ ].*','',field) for field in parse3]
    parse4 = sorted([field for field in parse3 if field not in excludes])

    # write
    out1.write('gene\n')
    for gene in parse4:
        out1.write(gene + '\n')

    
