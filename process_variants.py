#!/usr/bin/env python3
# on a sample by sample basis perform post-processing
# this is where we convert genomic coords if needed, filter variants

import os,pysam,re
import pandas as pd

inFile1 = 'oncoplex2cbio/annotation/oncoplex-genes.txt'
df = pd.read_csv(inFile1,sep="\t")
genes = list(df['gene'])
genes = dict(zip(genes,genes))
genesX = list(df[(df['chrom'] == 'chrX')]['gene'])
genesX = dict(zip(genesX,genesX))
chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrMT','chrM','chrY']

# *** process ***
def process(inFile,outFile,pipeline,overwrite = True):

    # check if file exists
    if os.path.exists(outFile) and not overwrite:
        print(' ... skipped',end='')
        return
    
#    pipeline = 'TGC2' if 'v7' in inFile or 'v8' in inFile else 'TGC'
    
    # ** mutation processing here **
    if 'mutation' in inFile:
        # create output directory if one is not already 
        outDir = os.path.dirname(outFile)
        if not os.path.exists(outDir):
            os.makedirs(outDir, exist_ok=True)

        consequences = ["coding", "frameshift", "inframe", "missense", "splicing-canonical", "start/stop",'frameshift deletion','frameshift insertion','nonframeshift deletion','nonsynonymous SNV','splicing','nonframeshift insertion']

            
        # read variant file
        df1 = pd.read_csv(inFile,sep="\t",keep_default_na = False)
        # ** filtering criteria **

        if pipeline == 'TGC2':
            idx0 = (df1['sheet'] == 'smallVariants')

            idx1 = (df1['humanFreq'] == 'NA') | (pd.to_numeric(df1['humanFreq'],errors='coerce') < 0.001)

            idx2 = (pd.to_numeric(df1['uwFreq'],errors='coerce') < 0.05)  # intwernal frequency database        

            idx3 = (df1['consequence'].isin(consequences))
            idx4 = (~df1['clinvar'].str.contains('benign',case = False))
            idx5 = (pd.to_numeric(df1['VAF'],errors='coerce') > 0.03)
            idx6 = (df1['filterRealVariant'] == True)
            idx10 = (pd.to_numeric(df1['AD'],errors='coerce') > 12)  
            # **
        
            idx7 = (df1['sheet'] == 'clinicallyFlagged') & (df1['clinvar'].str.contains('pathogenic',case=False)) & (~df1['clinvar'].str.contains('conflicting',case=False))
            idx8 = (pd.to_numeric(df1['AD'],errors='coerce') > 12)  # was 5 based on macro filters
            idx9 = (pd.to_numeric(df1['uwFreq'],errors='coerce') < 0.1)  # based on most frequent mutations in cancer ~5%
            df1a = df1[(idx0 & idx1 & idx2 & idx3 & idx4 & idx5 & idx6 & idx10) | (idx5 & idx7 & idx8 & idx9)].copy()
    
            # get rid of duplicates from small variants and clinically flagged
            df1b = df1a.drop_duplicates(subset=['id'],keep = 'first')
            df1b = df1b.replace('','NA')
            df1b.to_csv(outFile,index = False, sep="\t")
        # different filtering for older pipelines
        elif pipeline == 'TGC':
            idx0 = (df1['sheet'] == 'snpAnalysisTxt')
            idx1 = (df1['humanFreq'] == 'NA') | (pd.to_numeric(df1['humanFreq'],errors='coerce') < 0.001)

            idx2 = (df1['uwFreq'] == 'NA') | (pd.to_numeric(df1['uwFreq'],errors='coerce') < 0.05)  # internal frequency database

            # look for matches to consequence list
            regex = '|'.join(map(re.escape,consequences))
            idx3 = (df1['consequence'].str.contains(regex,na=False))

            idx4 = (~df1['clinvar'].str.contains('benign',case = False))
            idx5 = (df1['VAF'] == 'NA') | (pd.to_numeric(df1['VAF'],errors='coerce') > 0.03)
            idx10 = (pd.to_numeric(df1['AD'],errors='coerce') > 12)  

            df1a = df1[(idx0 & idx1 & idx2 & idx3 & idx4 & idx5 & idx10)].copy()
            
            # get rid of duplicates from small variants and clinically flagged
            df1b = df1a.drop_duplicates(subset=['id'],keep = 'first')
            df1b = df1b.replace('','NA')            
            df1b.to_csv(outFile,index = False, sep="\t")
        else:
            raise 'i no understand mutations'
            
        return df1b  # for troubleshooting

    # ** cnvs **
    elif 'cnv' in inFile:
        # create output directory if one is not already 
        outDir = os.path.dirname(outFile)
        if not os.path.exists(outDir):
            os.makedirs(outDir, exist_ok=True)

        # read variant file
        df1 = pd.read_csv(inFile,sep="\t",keep_default_na = False)

        df1['call'] = 'NA'
        df1['id'] = 'NA'
        df1['umbrellaID'] = 'NA'                
        df1['gisticValue'] = -10

        for idx,row in df1.iterrows():

            value = df1.iloc[idx]['avgLogRatio']
            value = value if value == 'NA' else float(value)
            del3 = df1.iloc[idx]['del3ExonValue']
            del3 = del3 if del3 == 'NA' else float(del3)
            gain3 = df1.iloc[idx]['gain3ExonValue']
            gain3 = gain3 if gain3 == 'NA' else float(gain3)
            absMax = df1.iloc[idx]['maxAbsLogRatio']
            absMax = absMax if absMax == 'NA' else float(absMax)
            sample = str(df1.iloc[idx]['sample'])
            gene = df1.iloc[idx]['gene']

            # thresholds based on interactions with Saji Hassan
            if value >= 1:
                call = 'AMP'
                gisticValue = 2
            elif value >= 0.6:   
                call = 'GAIN'
                gisticValue = 1
            elif value <= -1.3 or (pipeline == 'TGC2' and value <= -1 and absMax > 3):  # focal homodels
                call = 'HOMODEL'
                gisticValue = -2
            elif value <= -0.4:
                call = 'DEL'
                gisticValue = -1
            elif gain3 != 'NA' and value > 0 and gain3 >= 0.6:  # catch focal gain or amp events 
                call = 'GAIN'
                gisticValue = 1
            elif del3 != 'NA' and value < 0 and del3 <= -0.4:  # catch focal deletion events
                call = 'DEL'
                gisticValue = -1
            else:
                call = 'intact'
                gisticValue = 0
            df1.at[idx,'call'] = call
            df1.at[idx,'gisticValue'] = gisticValue
            df1.at[idx,'id'] = '__'.join([sample,gene,call])

            # get umbrellaID
            if call in ['AMP','GAIN']:
                umbrellaCall = 'gain'
            elif call in ['DEL','HOMODEL']:
                umbrellaCall = 'loss'
            else:
                umbrellaCall = 'intact'
            df1.at[idx,'umbrellaID'] = '__'.join([sample,gene,umbrellaCall])                 
        df1a = df1[['id','umbrellaID','sample','gene','avgLogRatio','maxAbsLogRatio','call','gisticValue','del3ExonValue','gain3ExonValue']]
        df1a.to_csv(outFile,index = False, sep='\t',na_rep = 'NA')
    # **

    # ** svs **
    elif 'svs' in inFile:
        # create output directory if one is not already 
        outDir = os.path.dirname(outFile)
        if not os.path.exists(outDir):
            os.makedirs(outDir, exist_ok=True)

        # read variant file
        df1 = pd.read_csv(inFile,sep="\t",keep_default_na = False,na_values = ['NA','nan'])
        df1['svLength'] = -1

        if pipeline == 'TGC2':
            # ** filtering criteria **
            idx1 = (df1['VAF'] > 0.008)
            idx2 = (df1['score'] > 1500) | ((df1['score'] > 400) & (df1['isQuiver'] | df1['isMitelman']))
            idx3 = (df1['filterRealVariant'] == True)
            # group 0 are consider more likely to be FP by dlmp
            idx4 = ((df1['group'] == 0) & (df1['isQuiver'] | df1['isMitelman'] | ~df1['id'].str.contains('Intergenic'))) | (df1['group'] != 0)
            df1a = df1[idx1 & idx2 & idx3 & idx4].copy()

            df1a = df1a.reset_index(drop = True)

        # **
        elif pipeline == 'TGC':
            idx1 = (df1['chrom1'].isin(chroms)) & (df1['chrom2'].isin(chroms))
            idx2 = (df1['gene1'].isin(genes)) | (df1['gene2'].isin(genes))
            idx3 = (df1['gene1'] != df1['gene2'])

            df1a = df1[idx1 & idx2 & idx3].copy()
            df1a = df1a.reset_index(drop = True)

        # ** add svLength **
        if len(df1a) > 0:        
            for idx,row in df1a.iterrows():
                chrom1 = df1a.iloc[idx]['chrom1']
                chrom2 = df1a.iloc[idx]['chrom2']            
                pos1 = df1a.iloc[idx]['pos1']
                pos2 = df1a.iloc[idx]['pos2']
                
                if pd.isna(pos1) or pos1 == 'NA' or pd.isna(pos2) or pos2 == 'NA':
                    svLength = pd.NA
                elif int(pos2) > int(pos1):
                    svLength = float(pos2) - float(pos1) + 1
                else:
                    svLength = float(pos1) - float(pos2) + 1
                df1a.at[idx,'svLength'] = -1 if chrom1 != chrom2 else svLength
        # ** 
        
        # get rid of duplicates
        df1b = df1a.drop_duplicates(subset=['id'],keep = 'first')
        df1b = df1b.replace('','NA')                    
        df1b.to_csv(outFile,index = False, sep="\t",na_rep = 'NA')
            
        return df1b # for troubleshooting
    # **
        
    # **
    else:
        print('i no understand')
        raise
    # **
