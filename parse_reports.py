#!/usr/bin/env python3

import os
import pandas as pd
import re
import openpyxl as op

# starting point for parsing reports
def parse(inFile,outDir):
    # make a home for new data tables
    if not os.path.exists(outDir):
        os.makedirs(outDir, exist_ok=True)
    
    # try to infer from the header which version of reports were are reading
    # call appopriate parse function
    df1 = pd.read_excel(inFile,keep_default_na = False, na_values = ['NA'])
    header = list(df1.columns)

    # initialize output files
    initOutput(outDir)

    # use arbitrary version marker for now...
    if header[3] == 'Position (hg38)':
        tmp = parseV8_0(df1,outDir)

    return(tmp)

# v8.0 oncoplex (arbitrarily 0)
def parseV8_0(df1, outDir):
    # update column names
    cols = list(df1.columns)
    cols[1] = 'sample'; cols[3] = 'position'; cols[4] = 'gene'
    df1.columns = cols
    header = list(df1.columns)

    # iterate through rows
    sample = None
    for idx,row in df1.iterrows():
        rowDict = row.to_dict()

        # identify the last non-empty column and turn all fields to strings
        lastCol = None
        for col in header:
            rowDict[col] = str(rowDict[col])
            if rowDict[col] != '':
                lastCol = col
                
        # get sample
        if rowDict['sample'] != '':
            sample = str(int(rowDict['sample']))
        else:
            rowDict['sample'] = sample
        assert sample is not None

        # determine variant type
        # if a mutation ...
        muts = [field for field in header if isinstance(rowDict[field], str) and not rowDict[field] == 'p.K802Rfs*11' and (('c.' in rowDict[field] or 'p.' in rowDict[field]) and field != 'clinical report') or rowDict[field] in ['chrX:67,546,506']]
        fusions = [field for field in header if isinstance(rowDict[field], str) and ('rearrangement' in rowDict[field] or 'fusion' in rowDict[field] or '::' in rowDict[field]) and field != 'clinical report']        
        isMutation = True if len(muts) > 0 else False
        isFusion = True if ':' in rowDict['gene'] or len(fusions) > 0 else False
        cnvs = [field for field in header if isinstance(rowDict[field], str) and (rowDict[field] in ['amplification','focal copy loss','single copy copy loss', 'single copy loss','focal copy gain','loss','biallelic loss'] or 'del' in rowDict[field] or 'LOH' in rowDict[field] or 'bi-allelic' in rowDict[field] or 'gain' in rowDict[field] or 'copy loss' in rowDict[field]) and field != 'clinical report']
        isCNV = True if len(cnvs) > 0 else False
        isIgnore = True if rowDict['gene'] in ['Numerous CNVs','Additional copy number alterations','Numerous additional mutations in the setting of high TMB (not listed here)','Numerous copy number alterations','Low tumor content','Additional copy number changes','Numerous mutations in the context of high TMB','Low tumor content limits the study'] or 'additional' in rowDict['gene'] else False
        tmbs = [field for field in header if 'TMB' in rowDict[field]]
        isTMB = True if len(tmbs) > 0 else False
        isBlank = True if rowDict['gene'] == '' and rowDict['HGVSp'] == '' and rowDict['position'] == '' else False
        
#        print(rowDict)
        if isMutation:
            status = 'mutation'
#            print(status)
            parseMutation(rowDict,outDir)
        elif isCNV:
            status = 'cnv'
#            print(status)
            parseCNV(rowDict,outDir)
        elif isFusion:
            status = 'fusion'
#            print(status)             
            parseFusion(rowDict,outDir)
        elif isIgnore:
            status = 'ignore'
#            print(status)                               
        elif isTMB:
            status = 'tmb'
#            print(status)
            parseTMB(rowDict,outDir)
        elif isBlank:
            status = 'blank'
#            print(status)
        else:
            print(rowDict)
            print('i no understand')
            raise
        
        outFile1 = outDir + '/report-parse-log.txt'    
        with open(outFile1,'a') as out1:
            lineOut = [sample,str(idx+1),status]
            out1.write('\t'.join(lineOut) + '\n')
        
    return(rowDict)

# parse and write mutation from report
def parseMutation(rowDict,outDir):

    header = ['id','sample','gene','mutationID','chrom','pos','ref','alt','hgvsP','hgvsC','type','transcript','vaf']

    # update position
    rowDict['chrom'] = re.sub('[:].*','',rowDict['position'])
    rowDict['pos'] = re.sub('.*[:]','',rowDict['position']).replace(',','')
    rowDict['hgvsP'] = rowDict['HGVSp']
    
    # PARSE C. annotation to get ref and alt and transcript
    # duplicate copies occur...
    tmp1 = list(set([rowDict[field] for field in rowDict.keys() if 'c.' in rowDict[field] and field != 'clinical report' and 'p.' in rowDict[field]]))
    tmp2 = list(set([rowDict[field] for field in rowDict.keys() if 'c.' in rowDict[field] and field != 'clinical report' and 'VAF' in rowDict[field]]))
    tmp3 = list(set([rowDict[field] for field in rowDict.keys() if 'c.' in rowDict[field] and 'NM_' in rowDict[field] and field != 'clinical report']))
    tmp4 = [rowDict[field] for field in rowDict.keys() if rowDict[field] == 'chrX:67,546,506']
    tmp5 = list(set([rowDict[field] for field in rowDict.keys() if field != 'clinical report' and rowDict[field].count('-') == 3]))

    if len(tmp1) == 1:
        tmp = tmp1
    elif len(tmp2) == 1:
        tmp = tmp2
    elif len(tmp3) == 1:
        tmp = tmp3
    elif len(tmp4) == 1:
        tmp = tmp4
    elif len(tmp5) == 1:
        tmp = tmp5
    else:
        print('i no understand mutation')
        raise
    assert len(tmp) == 1
    tmp = tmp[0]
    """
    if rowDict['sample'] == '78230' and rowDict['gene'] == '':
        print(tmp1)
        print(tmp2)
        print(tmp3)
        print(tmp4)               
        print(tmp5)       
        print(tmp)
        print(rowDict)
        raise
    """
    
    # get flavor
    # flavor 1
    if '(' in tmp and ')' in tmp:
        parse1= re.findall('.*[(].*(NM_[0-9.]+)[:](c.*), VAF ([0-9]+)[%)]', tmp)
        assert len(parse1[0]) == 3

        if '>' in parse1[0][1]:
            parse2 = re.findall('.*[0-9]+(.*)[>](.*)',parse1[0][1])
            assert len(parse2[0]) == 2
            rowDict['ref'] = parse2[0][0]
            rowDict['alt'] = parse2[0][1]
        else:
            rowDict['ref'] = 'NA'
            rowDict['alt'] = 'NA'
            
        rowDict['transcript'] = parse1[0][0]
        rowDict['hgvsC'] = parse1[0][1]
        rowDict['vaf'] = str(round(float(parse1[0][2])/100,3))

        type1 = None
        for key in rowDict.keys():
            if rowDict[key] in ['missense','frameshift','start/stop','splicing-canonical']:
                type1 = rowDict[key]
        rowDict['type'] = type1 if type1 is not None else 'NA'
    # flavor 2
    elif 'c.' in tmp:
        parse1 = re.findall('(NM_[0-9.]+)[:](c.*)',tmp)
        assert len(parse1[0]) == 2
        
        if '>' in parse1[0][1]:
            parse2 = re.findall('.*[0-9]+(.*)[>](.*)',parse1[0][1])
            assert len(parse2[0]) == 2
            rowDict['ref'] = parse2[0][0]
            rowDict['alt'] = parse2[0][1]
        else:
            rowDict['ref'] = 'NA'
            rowDict['alt'] = 'NA'
            
        rowDict['transcript'] = parse1[0][0]
        rowDict['hgvsC'] = parse1[0][1]

        type1 = None
        for key in rowDict.keys():
            if rowDict[key] in ['missense','frameshift','start/stop','splicing-canonical']:
                type1 = rowDict[key]
        rowDict['type'] = type1 if type1 is not None else 'NA'
        rowDict['vaf'] = 'NA'

        # ** check to make sure gene is not a single letter - issue with #78230 ** 
        # ** re-assign gene and hgvsP **
        if len(rowDict['gene']) == 1:
            key_hgvsP = [key for key in rowDict.keys() if re.findall('^p.',rowDict[key]) and key != 'clinical report']
            assert len(key_hgvsP) == 1
            key_hgvsP = key_hgvsP[0]
            
            # get previous key (assuming order is consistent)
            tmp = list(rowDict.keys())
            key_gene = tmp[tmp.index(key_hgvsP)-1]

            rowDict['hgvsP'] = rowDict[key_hgvsP]
            rowDict['gene'] = rowDict[key_gene]
        # **
        
    # flavor 3
    elif tmp in ['chrX:67,546,506']:
        rowDict['ref'] = rowDict['alt'] = rowDict['transcript'] = 'NA'
        rowDict['hgvsC'] = rowDict['type'] = rowDict['vaf'] = 'NA'
    # flavor 4
    elif tmp.count('-') == 3:
        parse1 = tmp.split('-')
        rowDict['chrom'] = 'chr' + parse1[0]
        rowDict['pos'] = parse1[1]
        rowDict['ref'] = parse1[2]
        rowDict['alt'] = parse1[3]
        rowDict['transcript'] = rowDict['hgvsC'] = rowDict['type'] = 'NA'
        rowDict['vaf'] = 'NA'

        # ** check to make sure gene is not a single letter - issue with #78230 ** 
        # ** re-assign gene and hgvsP **
        if len(rowDict['gene']) in [0,1]:
            key_hgvsP = [key for key in rowDict.keys() if re.findall('^p.',rowDict[key]) and key != 'clinical report']
            assert len(key_hgvsP) == 1
            key_hgvsP = key_hgvsP[0]
            
            # get previous key (assuming order is consistent)
            tmp = list(rowDict.keys())
            key_gene = tmp[tmp.index(key_hgvsP)-1]

            rowDict['hgvsP'] = rowDict[key_hgvsP]
            rowDict['gene'] = rowDict[key_gene]
        # **
    else:
        print('no understand mutation parse flavor')
        raise

    rowDict['mutationID'] = '_'.join([rowDict['chrom'],rowDict['pos'],rowDict['ref'],rowDict['alt']])
    
    if rowDict['hgvsC'] == '' or rowDict['hgvsC'] == 'NA':
        rowDict['id'] = '__'.join([rowDict['sample'],rowDict['mutationID']])
    else:
        rowDict['id'] = '__'.join([rowDict['sample'],rowDict['gene'],rowDict['hgvsC']])

    # write the parsed line
    outFile1 = outDir + '/report-mutations.txt'    
    with open(outFile1,'a') as out1:
        lineOut = []
        for field in header:
            lineOut.append(rowDict[field])
        lineOut = ['NA' if field == '' else field for field in lineOut]
        out1.write('\t'.join(lineOut) + '\n')        

# parse cnv row
def parseCNV(rowDict,outDir):
    header = ['sample','gene','type','call']    

    # find the call by considering calls that have been previously used
    tmpCall = [rowDict[field] for field in rowDict.keys() if (rowDict[field] in ['amplification','focal copy loss','single copy copy loss','single copy loss','focal copy gain','loss','biallelic loss'] or 'LOH' in rowDict[field] or 'del' in rowDict[field] or 'bi-allelic' in rowDict[field] or 'gain' in rowDict[field] or 'copy loss' in rowDict[field]) and field != 'clinical report']
    tmpCall = tmpCall if len(tmpCall) == 1 else [call for call in tmpCall if 'biallelic' in call]
    assert len(tmpCall) == 1
    call = tmpCall[0]

    # categorize call
    if call in ['amplification']:
        rowDict['call'] = 'AMP'
    elif call in ['focal copy gain'] or 'gain' in call:
        rowDict['call'] = 'GAIN'        
    elif call in ['focal copy loss','single copy copy loss','single copy loss','loss']:
        rowDict['call'] = 'DEL'
    elif ('del' in call and 'LOH' in call) or 'bi-allelic' in call or 'biallelic' in call:
        rowDict['call'] = 'HOMODEL'
        
    elif call == 'chr13 copy loss incl BRCA2':
        rowDict['call'] = 'DEL'
        rowDict['gene'] = 'BRCA2'
    else:
        print('i no understand cnv call')
        raise

    # deal with misplaced gene fun
    if rowDict['gene'] == '':
        key_call = [key for key in list(rowDict) if rowDict[key] == call]
        idx_key_call = list(rowDict).index(key_call[0])
        key_gene = list(rowDict)[idx_key_call-1]
        rowDict['gene'] = rowDict[key_gene]

    # assume gene level call until we learn otherwise
    rowDict['type'] = 'gene'
    if len(rowDict['gene']) == 1 or ',' in rowDict['gene']:
        rowDict['type'] = 'chrom'
        rowDict['gene'] = rowDict['gene'].replace(' ','')
        rowDict['gene'] = re.sub('[(].*[)]','',rowDict['gene'])

    # write the parsed line
    outFile1 = outDir + '/report-cnvs.txt'    
    with open(outFile1,'a') as out1:
        lineOut = []
        for field in header:
            lineOut.append(rowDict[field])
        out1.write('\t'.join(lineOut) + '\n')        

# parse fusion row
def parseFusion(rowDict,outDir):
    header = ['sample','gene1','gene2']

    # rearrangement?
    gene = None
    rearrangement = False
    for key in rowDict.keys():
        if 'rearrangement' in rowDict[key]:
            rearrangement = True
            break
        elif ':' in rowDict[key] and 'chr' not in rowDict[key]:
            rowDict['gene'] = rowDict[key]
        gene = rowDict[key]

    assert ':' in rowDict['gene'] or rearrangement

    if '::' in rowDict['gene']:
        rowDict['gene1'] = re.sub('[:][:].*','',rowDict['gene'])
        rowDict['gene2'] = re.sub('.*[:][:]','',rowDict['gene'])
    elif ':' in rowDict['gene']:
        rowDict['gene1'] = re.sub('[:].*','',rowDict['gene'])
        rowDict['gene2'] = re.sub('.*[:]','',rowDict['gene'])
    elif gene:
        rowDict['gene1'] = gene
        rowDict['gene2'] = 'NA'
    else:
        # find gene after coords - only a few cases like this
        key_coords = [key for key in rowDict.keys() if key != 'clinical report' and re.findall('^chr',rowDict[key])]
        assert len(key_coords) == 1
        key_coords = key_coords[0]
        idx_coords = list(rowDict.keys()).index(key_coords)
        key_gene = list(rowDict.keys())[idx_coords+1]
        rowDict['gene1'] = rowDict[key_gene]
        rowDict['gene2'] = 'NA'
    # write the parsed line
    outFile1 = outDir + '/report-fusions.txt'    
    with open(outFile1,'a') as out1:
        lineOut = []
        for field in header:
            lineOut.append(rowDict[field])
        lineOut = ['NA' if field == '' else field for field in lineOut]
        lineOut = [field.replace(' rearrangement','') for field in lineOut]        
        out1.write('\t'.join(lineOut) + '\n')        

# parse tmb field
def parseTMB(rowDict,outDir):
    header = ['sample','property','value']    

    # find tmb field
    for key in rowDict.keys():
        if 'TMB' in rowDict[key]:
            parse1 = rowDict[key]

    outFile1 = outDir + '/report-other.txt'
    with open(outFile1,'a') as out1:

        # if available get MSI status
        if 'MSI' in parse1 or 'MSS' in parse1:
            value = True if 'MSI' in parse1 else False
            lineOut = [rowDict['sample'], 'MSI', str(value)]
            out1.write('\t'.join(lineOut) + '\n')

        # get TMB status
        parse2 = re.findall('.*TMB[ -=]+([0-9]+)',parse1)
        parse2 = re.findall('([0-9]+)[)]*$',parse1)        
        
        assert len(parse2) == 1
        value = parse2[0]
        lineOut = [rowDict['sample'],'TMB',value]
        out1.write('\t'.join(lineOut) + '\n')
    
# create output files with headers
def initOutput(outDir):

    # mutations file
    outFile1 = outDir + '/report-mutations.txt'    
    with open(outFile1,'w') as out1:
        header = ['id','sample','gene','mutationID','chrom','pos','ref','alt','hgvsP','hgvsC','type','transcript','vaf']        
        out1.write('\t'.join(header) + '\n')

    # cnv file
    outFile1 = outDir + '/report-cnvs.txt'    
    with open(outFile1,'w') as out1:
        header = ['sample','gene','type','call']
        out1.write('\t'.join(header) + '\n')

    # fusion file
    outFile1 = outDir + '/report-fusions.txt'    
    with open(outFile1,'w') as out1:
        header = ['sample','gene1','gene2']
        out1.write('\t'.join(header) + '\n')

    # other file
    outFile1 = outDir + '/report-other.txt'    
    with open(outFile1,'w') as out1:
        header = ['sample','property','value']
        out1.write('\t'.join(header) + '\n')

    # log file - record what was done with each line
    outFile1 = outDir + '/report-parse-log.txt'    
    with open(outFile1,'w') as out1:
        header = ['sample','row','status']
        out1.write('\t'.join(header) + '\n')

        
