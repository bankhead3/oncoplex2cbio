#!/usr/bin/env python3
# parse macro workbooks provided by dlmp

import os
import pandas as pd
import re
import openpyxl as op

# starting point for parsing reports
def parse(inFile,outDir,sample):
    # make a home for new data tables
    if not os.path.exists(outDir):
        os.makedirs(outDir, exist_ok=True)

    # read first tab to infer workbook flavor
    df1 = pd.read_excel(inFile,keep_default_na = False, na_values = ['NA'])
    general = dict(zip(list(df1.iloc[:,1]),list(df1.iloc[:,2])))

    general['pipeline'] = '__'.join([general['pipeline_name'],general['pipeline_version'],general['assay']])
    print(general['pipeline'])

    if 'tgc2' in general['pipeline']:
        tmp = parseV8_0(inFile,outDir,sample)

# *** v8.0 oncoplex (arbitrarily 0) *** 
def parseV8_0(inFile, outDir, sample):

    # ** parse svs **
    print('Intergene SVs')
    outFile = outDir + '/svs-' + sample + '.txt'
    with open(outFile,'w') as out1:
        
        # write yo header
        header = ['id','sample','gene1','gene2','location','caller','eventType','score','VAF','AD','DP','isQuiver','isMitelman','filterRealVariant','chrom1','pos1','chrom2','pos2','filter']
        out1.write('\t'.join(header) + '\n')

        df1 = pd.read_excel(inFile,sheet_name = 'Intergene SVs', keep_default_na = False, na_values = ['NA'], skiprows = 2)

        for idx,row in df1.iterrows():
            rowDict = row.to_dict()        
            
            # assemble line from row read
            record = {'sample':sample}

            record['gene1'],record['gene2'] = rowDict['Gene1'],rowDict['Gene2']
            record['eventType'] = rowDict['Event Type']
            record['score'] = rowDict['QUAL']
            record['VAF'],record['AD'] = rowDict['VAF'],int(rowDict['ALT Fragments'])
            record['DP'] = record['AD'] + int(rowDict['REF Fragments'])
            record['isQuiver'] = True if rowDict['QUIVER'] != '' else False
            record['isMitelman'] = True if rowDict['MITELMAN'] != '' else False            
            record['filterRealVariant'] = rowDict['FILTER: REAL VARIANT']
            record['caller'],record['filter'] = rowDict['caller'],rowDict['FILTER']
            record['chrom1'],record['pos1'] = rowDict['Position1'].replace(',','').split(':')
            record['chrom2'],record['pos2'] = rowDict['Position2'].replace(',','').split(':')

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
                
    # ** parse cnvs **
    print('Copy Number by Gene')
    outFile1 = outDir + '/cnvs-' + sample + '.txt'
    with open(outFile1,'w') as out1:
        
        # write yo header
        header = ['sample','gene','avgLogRatio','maxAbsLogRatio']
        out1.write('\t'.join(header) + '\n')
        
        df1 = pd.read_excel(inFile,sheet_name = 'Copy Number by Gene', keep_default_na = False, na_values = ['NA'])
        for idx,row in df1.iterrows():
            rowDict = row.to_dict()

            # ** only include genes that are listed on the panel **
            if not rowDict['FILTER: GENE ORDERED']:
                continue
            # **
            
            record = {'sample':sample}

            record['gene'] = rowDict['gene']
            record['avgLogRatio'] = float(rowDict['gene average'])

            if rowDict['exon abs max'] != '':
                record['maxAbsLogRatio'] = float(rowDict['exon abs max'])
            else:
                record['maxAbsLogRatio'] = 'NA'
                
            # write yo line out yo
            lineOut = []
            for field in header:
                lineOut.append(str(record[field]))
            lineOut = ['NA' if field == '' else field for field in lineOut]
            out1.write('\t'.join(lineOut) + '\n')
            # **
    # **

    # ** parse mutations **
    outFile1 = outDir + '/mutations-' + sample + '.txt'
    with open(outFile1,'w') as out1:
        # write yo header
        header = ['id','sample','sheet','mutationID','filter','gene','chrom','pos','ref','alt','hgvsP','hgvsC','consequence','transcript','VAF','AD','DP','gnomadFreq','uwFreq','COSMIC','clinvar','filterRealVariant']
        out1.write('\t'.join(header) + '\n')

        sheets = ['Small Variants','Clinically Flagged']
        for sheet in sheets:
            print(sheet)
            # ** parse mutations ** probably needs its own function **
            df1 = pd.read_excel(inFile,sheet_name = sheet, keep_default_na = False, na_values = ['NA'],skiprows=2)
            for idx,row in df1.iterrows():
                rowDict = row.to_dict()
                record = {'sample':sample,'sheet':sheet}
                record['chrom'],record['pos'] = re.sub('[:].*','',rowDict['Position (hg38)']),re.sub('.*[:]','',rowDict['Position (hg38)']).replace(',','')
                record['ref'],record['alt'] = rowDict['ref'],rowDict['alt']
                record['gene'],record['transcript'] = rowDict['Gene'],re.sub('[:].*','',rowDict['HGVSc'])
                record['mutationID'] = '_'.join([record['chrom'],record['pos'],record['ref'],record['alt']])
                record['hgvsP'],record['hgvsC'] = rowDict['HGVSp'],re.sub('.*[:]','',rowDict['HGVSc'])
                record['consequence'],record['VAF'] = rowDict['Consequence'],rowDict['VAF']
                record['gnomadFreq'],record['COSMIC'] = rowDict['gnomADg AF'],rowDict['COSMIC']
                if sheet == 'Small Variants':
                    record['filter'],record['filterRealVariant'] = rowDict['FILTER'],rowDict['FILTER: REAL VARIANT']
                elif sheet == 'Clinically Flagged':
                    record['filter'],record['filterRealVariant'] = 'NA',rowDict['FILTER: REAL VARIANT']
                else:
                    print('i no understand')
                    raise
                
                record['uwFreq'],record['clinvar'] = rowDict['UW FREQ'],rowDict['ClinVar']
                rowDict['Alt Reads'] = 0 if rowDict['Alt Reads'] == '' else rowDict['Alt Reads']
                rowDict['Ref Reads'] = 0 if rowDict['Ref Reads'] == '' else rowDict['Ref Reads']            
                record['AD'],record['DP'] = rowDict['Alt Reads'],int(rowDict['Alt Reads']) + int(rowDict['Ref Reads'])
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
    # ** 
