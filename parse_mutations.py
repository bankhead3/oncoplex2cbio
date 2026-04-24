#!/usr/bin/env python3
# parse mutations

import os
import pandas as pd
import re
import openpyxl as op

import csv
from io import StringIO

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
    
#    print(fileDict)
#    print(fileDict['snpAnalysisTxt'])
#    print(fileDict['genotypeAnalysisTxt'])    
#    print(outFile)
    with open(outFile,'w') as out1:
        # this is the header we are aiming to parse for
        header = ['id','sample','sheet','mutationID','filter','gene','chrom','pos','ref','alt','hgvsP','hgvsC','consequence','transcript','VAF','AD','DP','humanFreq','uwFreq','COSMIC','clinvar','filterRealVariant','reference']
        out1.write('\t'.join(header) + '\n')
    
        if pipeline == 'TGC2' and version in ['V7','V8']:
            tmp = parseClinicallyFlaggedFeather(fileDict, 'clinicallyFlaggedFeather', sample, out1, header)
            tmp = parseClinicallyFlaggedFeather(fileDict, 'smallVariantsFeather', sample, out1,header)
        elif pipeline == 'TGC' and version in ['V5','V6','V7']:
            # we have to keep track of mutations to avoid dups
            mutations = {}
            tmp = parseAnalysisTxt(fileDict, 'snpAnalysisTxt', sample, out1, header, mutations)
            tmp = parseAnalysisTxt(fileDict, 'genotypeAnalysisTxt', sample, out1, header, mutations)  
#            print(fileDict)
#            print(fileDict.keys())
#            print(outFile)
        else:
            print('i no understant mutations for ' + pipeline + ' ' + version)
            raise
    print('...mutations:parsed',end='')

# *** parse analysis txt file ***
def parseAnalysisTxt(fileDict, fileType, sample, out1, header, mutations):
    
#    print(sample)
#    print('FIELDS WE WANT TO PARSE:')
#    print(header)
    sheet = fileType
    reference = 'GRCh37'

#    print(fileDict[fileType])
    
    with open(fileDict[fileType]) as in1:
        inHeader = in1.readline()
        inHeader = inHeader.strip().split('\t')

        for line in in1:
            # use don't use split because there are quoted strings with fun tabs
            parse1 = next(csv.reader(StringIO(line),delimiter='\t',quotechar='"'))
#            parse1 = line.split('\t')

            
            assert len(inHeader) == len(parse1)
            lineDict = dict(zip(inHeader,parse1))

            # build out record to include fields to parse
            record = {'sample':sample,'sheet':sheet,'reference':reference}
            
            record['chrom'],record['pos'] = lineDict['Position'].split(':')
            record['pos'] = record['pos'].split('-')[0]
            
            record['ref'],record['alt'] = lineDict['Ref_Base'],lineDict['Var_Base']
            record['mutationID'] = '_'.join([record['chrom'],record['pos'],record['ref'],record['alt']])
            record['filter'] = 'NA'
            if fileType == 'snpAnalysisTxt':
                record['gene'],record['transcript'] = lineDict['Gene'].replace('"',''),re.sub('[:].*','',lineDict['c.'])
                record['hgvsP'],record['hgvsC'] = lineDict['p.'],re.sub('.*[:]','',lineDict['c.'])                                
                record['consequence'] = lineDict['Variant_Type']
                record['VAF'],record['AD'] = lineDict['Allele_Frac'],lineDict['Var_Reads']
                record['DP'] = str(int(lineDict['Ref_Reads']) + int(lineDict['Var_Reads']))

                # take whatevery frequency is higher to be conservative
                if 'HiSeq_Freq' in lineDict:
                    record['uwFreq'] = max(float(lineDict['HiSeq_Freq']),float(lineDict['NextSeq_Freq']))
                elif 'UW_Freq' in lineDict:
                    record['uwFreq'] = lineDict['UW_Freq']
                else:
                    raise 'no understand freq'
                if record['uwFreq'] == -1:
                    record['uwFreq'] = 'NA'

                # v5 and v6 had different headers
                if 'Cosmic' in lineDict:
                    record['COSMIC'],record['clinvar'] = lineDict['Cosmic'],lineDict['ClinVar']
                elif 'Cosmic_ID' in lineDict:
                    record['COSMIC'],record['clinvar'] = lineDict['Cosmic_ID'],lineDict['CLNSIG']
                else:
                    raise 'i no understand how to get cosmic'
                    
                record['filterRealVariant'] = 'NA'
                record['humanFreq'] = lineDict['EXAC']
            else:
#                print(lineDict['Clinically_Flagged'])
                flag = lineDict['Clinically_Flagged'].split(' ')
                gene = flag[0].replace('"','')
                # look for p. or amino acid start
                hgvsP = [field for field in flag if 'p.' in field or re.findall('^[ARNDCQEGHILKMFPSTWYV]',field) and field != field[0]]
                hgvsP = hgvsP[0] if len(hgvsP) == 1 else 'NA'
                

                hgvsC = [field for field in flag if 'c.' in field]
                hgvsC = hgvsC[0] if len(hgvsC) == 1 else 'NA'
                assert gene != '-'
                record['gene'],record['hgvsP'],record['hgvsC'] = gene,hgvsP,hgvsC
                record['transcript'] = 'NA'

                record['consequence'] = 'NA'
                record['DP'],record['AD'] = lineDict['Valid_Reads'],lineDict['Variant_Reads']
                if int(record['DP']) != 0:
                    record['VAF'] = str(round(float(record['AD'])/float(record['DP']),3))
                else:
                    record['VAF'] = str(0)

                record['uwFreq'] = 'NA'
                record['COSMIC'],record['clinvar'] = 'NA','NA'
                record['filterRealVariant'] = 'NA'
                record['humanFreq'] = 'NA'
                

            if record['hgvsC'] == '' or record['hgvsC'] == 'NA':
                record['id'] = '__'.join([record['sample'],record['mutationID']])
            else:
                record['id'] = '__'.join([record['sample'],record['gene'],record['hgvsC']])


            # ignore duplicate mutations (they present in genotyping
            if record['id'] in mutations:
                continue
            else:
                mutations[record['id']] = record['id']
                
            # write yo line out
            lineOut = []
            for field in header:
                lineOut.append(str(record[field]))
            lineOut = [field.strip() for field in lineOut]                
            lineOut = ['NA' if field == '' else field for field in lineOut]
            out1.write('\t'.join(lineOut) + '\n')
            # **

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
        record['pos'] = record['pos'].split('-')[0]
        
        record['ref'],record['alt'] = rowDict['ref'],rowDict['alt']
        record['gene'],record['transcript'] = rowDict['Gene'].replace('"',''),re.sub('[:].*','',rowDict['HGVSc'])
                
        record['mutationID'] = '_'.join([record['chrom'],record['pos'],record['ref'],record['alt']])
        record['hgvsP'],record['hgvsC'] = rowDict['HGVSp'],re.sub('.*[:]','',rowDict['HGVSc'])
        
        record['consequence'],record['VAF'] = rowDict['Consequence'],rowDict[sample + '.VAF']
        if 'gnomAD_AF' in rowDict:
            record['humanFreq'],record['COSMIC'] = rowDict['gnomAD_AF'],rowDict['COSMIC']['text']
        elif 'gnomADg_AF' in rowDict:
            record['humanFreq'],record['COSMIC'] = rowDict['gnomADg_AF'],rowDict['COSMIC']['text']            

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

