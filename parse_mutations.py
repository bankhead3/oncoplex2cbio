#!/usr/bin/env python3
# parse mutations

import os
import pandas as pd
import re
import openpyxl as op
# import pyarrow.feather as feather

# starting point for parsing reports
def parse(fileDict,outDir,sample,version,pipeline,overwrite=False):
    # make a home for new data tables
    if not os.path.exists(outDir):
        os.makedirs(outDir, exist_ok=True)

    outFile = outDir + '/mutations.txt'

    # check if file exists
    if os.path.exists(outFile) and not overwrite:
        print('...mutations:skipped',end='')
        return

    with open(outFile,'w') as out1:
        # this is the header we are aiming to parse for
        header = ['id','sample','sheet','mutationID','filter','gene','chrom','pos','ref','alt','hgvsP','hgvsC','consequence','transcript','VAF','AD','DP','gnomadFreq','uwFreq','COSMIC','clinvar','filterRealVariant','reference']
        out1.write('\t'.join(header) + '\n')
    
        if pipeline == 'TGC2' and version in ['V7','V8']:
            tmp = parseClinicallyFlaggedFeather(fileDict, 'clinicallyFlaggedFeather', sample, out1, header)
            tmp = parseClinicallyFlaggedFeather(fileDict, 'smallVariantsFeather', sample, out1,header)            
        else:
            print('i no understant mutations for ' + pipeline + ' ' + version)
            raise
    print('...mutations:parsed',end='')

# *** parse clinically flagged feather ***
def parseClinicallyFlaggedFeather(fileDict, fileType, sample, out1, header):
    
#    print(sample)
#    print(header)
    sheet = fileType.replace('Feather','')
    reference = 'GRCh37'
    
    df1 = pd.read_feather(fileDict[fileType])
    df1 = df1.fillna('NA')    
    for idx,row in df1.iterrows():

        rowDict = row.to_dict()
#        print(rowDict)

        record = {'sample':sample,'sheet':sheet,'reference':reference}
        # get hg19 pos
        tmp = rowDict['Position_(hg19)']
        record['chrom'],record['pos'] = re.sub('[:].*','',tmp['text']),re.sub('.*[:]','',tmp['text']).replace(',','')
        record['ref'],record['alt'] = rowDict['ref'],rowDict['alt']
        record['gene'],record['transcript'] = rowDict['Gene'],re.sub('[:].*','',rowDict['HGVSc'])
                
        record['mutationID'] = '_'.join([record['chrom'],record['pos'],record['ref'],record['alt']])
        record['hgvsP'],record['hgvsC'] = rowDict['HGVSp'],re.sub('.*[:]','',rowDict['HGVSc'])
        
        record['consequence'],record['VAF'] = rowDict['Consequence'],rowDict[sample + '.VAF']
        if 'gnomAD_AF' in rowDict:
            record['gnomadFreq'],record['COSMIC'] = rowDict['gnomAD_AF'],rowDict['COSMIC']['text']
        elif 'gnomADg_AF' in rowDict:
            record['gnomadFreq'],record['COSMIC'] = rowDict['gnomADg_AF'],rowDict['COSMIC']['text']            

        if sheet == 'smallVariants':
            record['filter'],record['filterRealVariant'] = rowDict['FILTER'],rowDict['FILTER:_REAL_VARIANT']
        elif sheet == 'clinicallyFlagged':
            record['filter'],record['filterRealVariant'] = 'NA',rowDict['FILTER:_REAL_VARIANT']
        else:
            print('i no understand')
            raise
        
        # parse cosmic
        if ';' in record['COSMIC']:
            record['COSMIC'] = record['COSMIC'].split(';')[0].replace('ID=','')

        record['uwFreq'],record['clinvar'] = rowDict['UW_FREQ'],rowDict['ClinVar']['text']
        rowDict['Alt Reads'] = 0 if rowDict[sample + '.Alt_Reads'] in ['','NA'] else rowDict[sample + '.Alt_Reads']
        rowDict['Ref Reads'] = 0 if rowDict[sample + '.Ref_Reads'] in ['','NA'] else rowDict[sample + '.Ref_Reads']

        record['AD'],record['DP'] = int(rowDict['Alt Reads']),int(rowDict['Alt Reads']) + int(rowDict['Ref Reads'])
        record['VAF'] = 0 if record['VAF'] == '' or record['VAF'] == 'NA' else record['VAF']

        if record['hgvsC'] == '' or record['hgvsC'] == 'NA':
            record['id'] = '__'.join([record['sample'],record['mutationID']])
        else:
            record['id'] = '__'.join([record['sample'],record['gene'],record['hgvsC']])

        # write yo line out
        lineOut = []
        for field in header:
            lineOut.append(str(record[field]))
        lineOut = ['NA' if field == '' else field for field in lineOut]
        out1.write('\t'.join(lineOut) + '\n')
        # **

#        print()
#        print()
#        print(rowDict.keys())
#        print()
#        print()
#        print(record)

