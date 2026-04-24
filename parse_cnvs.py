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
    outFile = outDir + '/cnvs.txt'

    # check if file exists
    if os.path.exists(outFile) and not overwrite:
        print('...cnvs:skipped',end='')
        return

    with open(outFile,'w') as out1:
        # this is the header we are aiming to parse for
        header = ['sample','gene','avgLogRatio','maxAbsLogRatio','del3ExonValue','gain3ExonValue']
        out1.write('\t'.join(header) + '\n')

        if pipeline == 'TGC2' and version in ['V7','V8']:
            tmp = parseExcelv8(fileDict, 'excel', sample, out1, header)
        elif pipeline == 'TGC' and version in ['V5','V6','V7']:
            tmp = parseTxt(fileDict, 'cnvGeneTxt', sample, out1, header)
        else:
            print('i no understant cnvs for ' + pipeline + ' ' + version)
            raise
    print('...cnvs:parsed',end='')

# *** parse cnvGeneTxt file ***
def parseTxt(fileDict, fileType, sample, out1, header):
    
    with open(fileDict[fileType]) as in1:
        inHeader = in1.readline()
        inHeader = inHeader.strip().split('\t')

        for line in in1:
            parse1 = line.split('\t')
            parse1 = [field.strip() for field in parse1]

            assert len(inHeader) == len(parse1)
            lineDict = dict(zip(inHeader,parse1))

            record = {'sample':sample}
            record['gene'] = lineDict['Gene']
            record['avgLogRatio'] = lineDict['Ave_Adjusted_Log_Ratio']

            # write yo line out yo
            lineOut = []
            for field in header:
                if field in record:
                    lineOut.append(str(record[field]))
                else:
                    lineOut.append('NA')
            lineOut = ['NA' if field == '' else field for field in lineOut]
            out1.write('\t'.join(lineOut) + '\n')
            # **

# *** parse clinically flagged feather ***
def parseExcelv8(fileDict, fileType, sample, out1, header):
    
    df1 = pd.read_excel(fileDict[fileType],sheet_name = 'Copy Number by Gene', keep_default_na = False, na_values = ['NA'])
    for idx,row in df1.iterrows():

        rowDict = row.to_dict()

        # ** only include genes that are listed on the panel **
        if not rowDict['FILTER: GENE ORDERED']:
            continue
        # **

        # get exons into a list numeric form
        # warning: top exon number is hard coded
        exons = [field for field in rowDict.keys() if field.isdigit()]
        values = [rowDict[exon] for exon in exons if rowDict[exon] != '']

        # ** sliding window to identify focal events spanning 3 consecutive exons **
        def slide(values,num):

            if num > len(values) or len(values) == 0:
                return ['NA','NA']
                
            gainWindowGlobal = min(values)
            delWindowGlobal = max(values)
                
            for i in range(len(values) - num + 1):
                window = values[i:(i + num)]
                maxWindow = max(window)
                minWindow = min(window)

                if minWindow > gainWindowGlobal:
                    gainWindowGlobal = minWindow
                if maxWindow < delWindowGlobal:
                    delWindowGlobal = maxWindow
                        
            return [delWindowGlobal,gainWindowGlobal]
            # ** 
        del3ExonValue,gain3ExonValue = slide(values,3)
            
        record = {'sample':sample,'del3ExonValue':del3ExonValue,'gain3ExonValue':gain3ExonValue}
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
