#!/usr/bin/env python3
# generate table of mane gene transcripts

import pandas as pd
import re

inFile1 = '../input/MANE.GRCh38.v1.4.summary.txt'
outFile1 = '../mane-gene-transcript.txt'

df1 = pd.read_csv(inFile1,sep="\t")


with open(outFile1,'w') as out1:
    
    # write yo header
    header = ['gene','transcript','entrez','ensemblGene','ensemblTranscript','chrom','start','end','strand']
    out1.write('\t'.join(header) + '\n')

    for idx,row in df1.iterrows():
        rowDict = row.to_dict()

        # for now for simplicity we skip mane clinical transcripts so we avoid dups
        if rowDict['MANE_status'] == 'MANE Plus Clinical':
            continue
        
        rowDict['gene'] = rowDict['symbol']
        rowDict['transcript'] = rowDict['RefSeq_nuc']
        rowDict['entrez'] = rowDict['#NCBI_GeneID'].replace('GeneID:','')
        rowDict['ensemblGene'],rowDict['ensemblTranscript'] = rowDict['Ensembl_Gene'],rowDict['Ensembl_nuc']

        if 'NC_' in rowDict['GRCh38_chr']:
            parse1 = re.findall('NC_[0]+([0-9XY]+)[.][0-9]+$',rowDict['GRCh38_chr'])
            assert len(parse1) == 1
            parse1 = 'chr' + parse1[0]
        else:
            parse1 = rowDict['GRCh38_chr']

        rowDict['chrom'],rowDict['start'] = parse1,rowDict['chr_start']
        rowDict['end'],rowDict['strand'] = rowDict['chr_end'],rowDict['chr_strand']

        # get yo line out and write yo
        lineOut = []
        for field in header:
            lineOut.append(str(rowDict[field]))
        out1.write('\t'.join(lineOut) + '\n')
        
    


