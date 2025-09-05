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

# *** process ***
def process(inFile,outFile):

    # ** mutation processing here ** 
    if 'mutation' in inFile:
        # create output directory if one is not already 
        outDir = os.path.dirname(outFile)
        if not os.path.exists(outDir):
            os.makedirs(outDir, exist_ok=True)

        # read variant file
        df1 = pd.read_csv(inFile,sep="\t",keep_default_na = False)

        # ** filtering criteria **
        idx0 = (df1['sheet'] == 'Small Variants')
        idx1 = (df1['gnomadFreq'] == 'NA') | (pd.to_numeric(df1['gnomadFreq'],errors='coerce') < 0.001)
        idx2 = (df1['uwFreq'] < 0.05)  # internal frequency database
        consequences = ["coding", "frameshift", "inframe", "missense", "splicing-canonical", "start/stop"]    
        idx3 = (df1['consequence'].isin(consequences))
        idx4 = (~df1['clinvar'].str.contains('benign',case = False))
        idx5 = (df1['VAF'] > 0.03)
        idx6 = (df1['filterRealVariant'] == True)
        idx10 = (df1['AD'] > 12)  
        # **

        idx7 = (df1['sheet'] == 'Clinically Flagged') & (df1['clinvar'].str.contains('pathogenic',case=False)) & (~df1['clinvar'].str.contains('conflicting',case=False))
        idx8 = (df1['AD'] > 12)  # was 5 based on macro filters
        idx9 = (df1['uwFreq'] < 0.1)  # based on most frequent mutations in cancer ~5%
        df1a = df1[(idx0 & idx1 & idx2 & idx3 & idx4 & idx5 & idx6 & idx10) | (idx5 & idx7 & idx8 & idx9)].copy()
    
        # get rid of duplicates from small variants and clinically flagged
        df1b = df1a.drop_duplicates(subset=['id'],keep = 'first')
        df1b.to_csv(outFile,index = False, sep="\t")

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
            elif value <= -1.3 or (value <= -1 and absMax > 3):  # focal homodels
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
        df1a.to_csv(outFile,index = False, sep='\t')
    # **

    # ** svs **
    elif 'svs' in inFile:
        # create output directory if one is not already 
        outDir = os.path.dirname(outFile)
        if not os.path.exists(outDir):
            os.makedirs(outDir, exist_ok=True)

        # read variant file
        df1 = pd.read_csv(inFile,sep="\t",keep_default_na = False)
        df1['svLength'] = -1
            
        # ** filtering criteria **
        idx1 = (df1['VAF'] > 0.008)
        idx2 = (df1['score'] > 1500) | ((df1['score'] > 400) & (df1['isQuiver'] | df1['isMitelman']))
        idx3 = (df1['filterRealVariant'] == True)
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
                    svLength = 'NA'
                elif pos2 > pos1:
                    svLength = float(pos2) - float(pos1) + 1
                else:
                    svLength = float(pos1) - float(pos2) + 1
                df1a.at[idx,'svLength'] = -1 if chrom1 != chrom2 else svLength
        # **
        
        # get rid of duplicates
        df1b = df1a.drop_duplicates(subset=['id'],keep = 'first')
        df1b.to_csv(outFile,index = False, sep="\t")
            
        return df1b # for troubleshooting
    # **
        
    # **
    else:
        print('i no understand')
        raise
    # **
    
