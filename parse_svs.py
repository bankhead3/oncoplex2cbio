#!/usr/bin/env python3
# parse mutations

import os
import pandas as pd
import re
import openpyxl as op

# starting point for parsing reports
def parse(fileDict,outDir,sample,version,pipeline,overwrite=False):
    # make a home for new data tables
    if not os.path.exists(outDir):
        os.makedirs(outDir, exist_ok=True)
    outFile = outDir + '/svs.txt'

    # check if file exists
    if os.path.exists(outFile) and not overwrite:
        print('...svs:skipped',end='')
        return

    with open(outFile,'w') as out1:
        # this is the header we are aiming to parse for
        header = ['id','sample','gene1','gene2','location','caller','eventType','score','VAF','AD','DP','isQuiver','isMitelman','filterRealVariant','chrom1','pos1','chrom2','pos2','filter','reference','group']
        out1.write('\t'.join(header) + '\n')

        if pipeline == 'TGC2' and version in ['V7','V8']:
            tmp = parseExcelv8(fileDict, 'interGeneSvsFeather', sample, out1, header)
        else:
            print('i no understant cnvs for ' + pipeline + ' ' + version)
            raise
    print('...svs:parsed',end='')

# *** parse clinically flagged feather ***
def parseExcelv8(fileDict, fileType, sample, out1, header):
    
#    df1 = pd.read_excel(fileDict[fileType],sheet_name = 'Intergene SVs', keep_default_na = False, na_values = ['NA'], skiprows = 2)

    df1 = pd.read_feather(fileDict[fileType])
    
    # read each row
    for idx,row in df1.iterrows():
        rowDict = row.to_dict()        
            
        # assemble line from row read
        record = {'sample':sample}

        record['gene1'],record['gene2'] = rowDict['Gene1'],rowDict['Gene2']
        record['eventType'] = rowDict['Event_Type']
        record['score'] = rowDict['QUAL']
        record['group'] = rowDict['Group']
        record['VAF'],record['AD'] = rowDict[sample + '.VAF'],int(rowDict[sample + '.ALT_Fragments'])
        record['DP'] = record['AD'] + int(rowDict[sample + '.REF_Fragments'])
        record['isQuiver'] = True if rowDict['QUIVER'] != '' else False
        record['isMitelman'] = True if rowDict['MITELMAN'] != '' else False            
        record['filterRealVariant'] = rowDict['FILTER:_REAL_VARIANT']
        record['caller'],record['filter'] = rowDict['caller'],rowDict['FILTER']

        # source from hg19 to make cbio happy
        record['reference'] = 'GRCh37'
        if rowDict['hg19_Coordinates1']:
            record['chrom1'],record['pos1'] = rowDict['hg19_Coordinates1'].split(':')
        else:
            record['chrom1'],record['pos1'] = 'NA','NA'
        if rowDict['hg19_Coordinates2']:
            record['chrom2'],record['pos2'] = rowDict['hg19_Coordinates2'].split(':')
        else:
            record['chrom2'],record['pos2'] = 'NA','NA'
                
        # fix the order if needed
        # if tmprss2 erg then make sure order is correct
        if record['gene1'] in ['TMPRSS2','ERG'] and record['gene2'] in ['TMPRSS2','ERG']:
            if record['gene1'] != 'TMPRSS2':
                record['gene1'] = 'TMPRSS2'
                record['gene2'] = 'ERG'

                tmp = record['pos1']
                record['pos1'] = record['pos2']
                record['pos2'] = tmp
                    
        # otherwise alphabetical
        elif record['gene2'] < record['gene1'] and record['gene2'] != 'Intergenic':
            tmp = record['gene1']
            record['gene1'] = record['gene2']
            record['gene2'] = tmp

            tmp = record['chrom1']
            record['chrom1'] = record['chrom2']
            record['chrom2'] = tmp
                
            tmp = record['pos1']
            record['pos1'] = record['pos2']
            record['pos2'] = tmp

        record['id'] = '__'.join([record['sample'],record['gene1'],record['gene2']])
        record['location'] = '_'.join([record['chrom1'],record['pos1'],record['chrom2'],record['pos2']])
            
        # write it yo
        lineOut = []
        for field in header:
            lineOut.append(str(record[field]))
        out1.write('\t'.join(lineOut) + '\n')

