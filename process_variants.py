#!/usr/bin/env python3
# on a sample by sample basis perform post-processing
# this is where we convert genomic coords if needed, filter variants

import os,pysam,re
import pandas as pd

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
        idx1 = (df1['gnomadFreq'] == 'NA') | (pd.to_numeric(df1['gnomadFreq'],errors='coerce') < 0.001)
        idx2 = (df1['uwFreq'] < 0.05)  # internal frequency database
        consequences = ["coding", "frameshift", "inframe", "missense", "splicing-canonical", "start/stop"]    
        idx3 = (df1['consequence'].isin(consequences))
        idx4 = (~df1['clinvar'].str.contains('benign',case = False))
        idx5 = (df1['VAF'] > 0.03)
        idx6 = (df1['filterRealVariant'] == True)
        df1a = df1[idx1 & idx2 & idx3 & idx4 & idx5 & idx6].copy()
        # **

        idx1 = (df1a['sheet'] == 'Clinically Flagged') & (df1a['clinvar'].str.contains('pathogenic',case=False)) & (~df1a['clinvar'].str.contains('conflicting',case=False))
        idx2 = (df1a['sheet'] == 'Small Variants')
        df1b = df1a[idx1 | idx2].copy()
    
        # get rid of duplicates from small variants and clinically flagged
        df1c = df1b.drop_duplicates(subset=['id'],keep = 'first')
        df1c.to_csv(outFile,index = False, sep="\t")

        return df1c  # for troubleshooting

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
        df1['gisticValue'] = -10

        for idx,row in df1.iterrows():

            value = df1.iloc[idx]['avgLogRatio']
            sample = str(df1.iloc[idx]['sample'])
            gene = df1.iloc[idx]['gene']

            if value >= 1:
                call = 'AMP'
                gisticValue = 2
            elif value >= 0.4:
                call = 'GAIN'
                gisticValue = 1
            elif value <= -1.3:
                call = 'HOMODEL'
                gisticValue = -2
            elif value <= -0.4:
                call = 'DEL'
                gisticValue = -1
            else:
                call = 'intact'
                gisticValue = 0
            df1.at[idx,'call'] = call
            df1.at[idx,'gisticValue'] = gisticValue
            df1.at[idx,'id'] = '_'.join([sample,gene,call])

        df1a = df1[['id','sample','gene','avgLogRatio','call','gisticValue']]
        df1a.to_csv(outFile,index = False, sep='\t')
    # **

    # **
    else:
        print('i no understand')
        raise
    # **
    
