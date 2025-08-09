#!/usr/bin/env python3
# generates cbioportal formatted files

import os,re
import pandas as pd

def build(fileDict,outDir):

    # create output directory if one is not already 
    if not os.path.exists(outDir):
        os.makedirs(outDir, exist_ok=True)

    if 'params' in fileDict:
        writeMetas(fileDict['params'],outDir)
    else:
        print('ERROR: params file is required')
        raise

    if 'mutations' in fileDict:
        writeMutations(fileDict['mutations'],outDir)
    if 'metadata' in fileDict:
        writeMetaDatas(fileDict['metadata'],outDir)
    if 'cnvs' in fileDict:
        writeCnvs(fileDict['cnvs'],outDir)
    if 'svs' in fileDict:
        writeSvs(fileDict['svs'],outDir)
    if 'panel' in fileDict:
        writePanel(fileDict,outDir)            

    # write case lists
    if 'caseLists' in fileDict:
        writeCaseLists(fileDict,outDir)

# *** write gene panel matrix data file ***
def writePanel(fileDict,outDir):
    # get all samples
    df1 = pd.read_csv(fileDict['caseListIds'],sep="\t")
    samples = sorted(list(df1[df1['caseList'] == 'all']['sample']))

    lines = open(fileDict['panel']).read().splitlines()
    panel = re.sub('.*[ ]','',lines[0])

    outFile = outDir + '/data_gene_panel_matrix.txt'
    with open(outFile,'w') as out1:
        header = ['SAMPLE_ID','mutations','structural_variants']
        out1.write('\t'.join(header) + '\n')

        # for each sample
        for sample in samples:
            lineOut = [str(sample),panel,panel]
            out1.write('\t'.join(lineOut) + '\n')
    
# *** write case list files ***
def writeCaseLists(fileDict,outDir):
    
    # create output directory if one is not already 
    outDir1 = outDir + 'case_lists/'
    if not os.path.exists(outDir1):
        os.makedirs(outDir1, exist_ok=True)

    # get study id
    df1 = pd.read_csv(fileDict['caseLists'],sep="\t")
    df2 = pd.read_csv(fileDict['caseListIds'],sep="\t")

    for name,groupDf in df1.groupby('caseList'):

        # ** get case list ids to write below  **
        samples = list(df2[df2['caseList'] == name]['sample'])
        case_list_ids = [str(sample) for sample in samples]
        case_list_ids[0] = 'case_list_ids: ' + case_list_ids[0]
        case_list_ids = '\t'.join(case_list_ids)
        # **
        
        outFile = outDir1 + 'cases_' + name + '.txt'
        with open(outFile,'w') as out1:        
            out1.write('cancer_study_identifier: ' + groupDf['cancer_study_identifier'].iloc[0] + '\n')
            out1.write('stable_id: ' + groupDf['stable_id'].iloc[0] + '\n')
            out1.write('case_list_name: ' + groupDf['case_list_name'].iloc[0] + '\n')
            out1.write('case_list_description: ' + groupDf['case_list_description'].iloc[0] + '\n')
            out1.write('case_list_category: ' + groupDf['case_list_category'].iloc[0] + '\n')
            out1.write(case_list_ids + '\n')
        
# *** write metadata data header ***
def writeMetadataHeader(patientSample,df1a,outFile):
    with open(outFile,'w') as out1:
        if patientSample == 'patient':
            out1.write('\t'.join(['#Patient Identifier'] + list(df1a['name'])) + '\n')
            out1.write('\t'.join(['#Patient identifier'] + list(df1a['description'])) + '\n')
            # datatype row
            row3 = ['#STRING'] + list(df1a['dataType'])            
            out1.write('\t'.join(row3) + '\n')
            
            row4 = ['1' if tmp > 0 else '#1' for tmp in range(df1a.shape[0] + 1)] 
            out1.write('\t'.join(row4) + '\n')
            out1.write('\t'.join(['PATIENT_ID'] + list(df1a['property'])) + '\n')
        elif patientSample == 'sample':
            out1.write('\t'.join(['#Patient Identifier'] + list(df1a['name'])) + '\n')
            out1.write('\t'.join(['#Patient identifier'] + list(df1a['description'])) + '\n')
            # datatype row
            row3 = ['#STRING'] + list(df1a['dataType'])
            out1.write('\t'.join(row3) + '\n')

            row4 = ['1' if tmp > 0 else '#1' for tmp in range(df1a.shape[0] + 1)] 
            out1.write('\t'.join(row4) + '\n')
            out1.write('\t'.join(['PATIENT_ID'] + list(df1a['property'])) + '\n')
            
# *** write metadata data files ***
def writeMetaDatas(inFile,outDir):
    df1 = pd.read_csv(inFile,sep="\t")
    df1a = df1[df1['patientSample'] == 'patient'].copy()

    # *** write data_clinical_patient ***
    df1b = df1a.drop_duplicates(subset=['property'])
    outFile = outDir + 'data_clinical_patient.txt'
    writeMetadataHeader('patient',df1b,outFile)

    with open(outFile,'a') as out1:
        for name,groupDf in df1a.groupby('id'):
            properties = list(groupDf['property'])
            values = list(groupDf['value'])
            lineDict = dict(zip(properties,values))
            lineDict['id'] = name

            # assemble and writ yo line out
            lineOut = []
            for field in ['id'] + properties:
                lineOut.append(lineDict[field])
            out1.write('\t'.join(lineOut) + '\n')
    # ***

    # *** write data_clinical_sample ***
    df1b = df1[df1['patientSample'] == 'sample'].copy()
    outFile = outDir + 'data_clinical_sample.txt'    
    df1c = df1b.drop_duplicates(subset=['property'])
    outFile = outDir + 'data_clinical_sample.txt'
    writeMetadataHeader('sample',df1c,outFile)

    with open(outFile,'a') as out1:
        for name,groupDf in df1b.groupby('id'):
            properties = list(groupDf['property'])
            values = list(groupDf['value'])
            lineDict = dict(zip(properties,values))
            lineDict['patientID'] = list(groupDf['patientID'])[0]

            # assemble and writ yo line out
            lineOut = []
            for field in ['patientID'] + properties:
                lineOut.append(lineDict[field])
            out1.write('\t'.join(lineOut) + '\n')
    # ***
    
# *** write mutations file ***
def writeMutations(inFile,outDir):
    
    outFile = outDir + 'data_mutations.txt'
    with open(inFile) as in1, open(outFile,'w') as out1:

        header = in1.readline()
        header = header.strip().split('\t')

        # ** write out header **
        outHeader = ['Hugo_Symbol','NCBI_Build','Chromosome','Start_Position','End_Position','Variant_Classification','Reference_Allele','Tumor_Seq_Allele2','Tumor_Sample_Barcode','HGVSp_Short','t_alt_count','t_ref_count']
        out1.write('\t'.join(outHeader) + '\n')
        
        for line in in1:
            parse1 = line.strip().split('\t')
            lineDict = dict(zip(header,parse1))

            record = {'Chromosome':lineDict['chrom'].replace('chr','')}
            record['Start_Position'] = lineDict['pos']
            record['End_Position'] = lineDict['pos']

            if len(lineDict['ref']) > 1 and lineDict['pos'] != 'NA':
                record['End_Position'] = str(int(lineDict['pos']) + len(lineDict['ref']) - 1)
                
            record['Reference_Allele'] = lineDict['ref']
            record['Tumor_Seq_Allele2'] = lineDict['alt']
            record['Tumor_Sample_Barcode'] = lineDict['sample']
            record['HGVSp_Short'] = lineDict['hgvsP']
            record['HGVSc'] = lineDict['hgvsC']
            record['t_ref_count'] = str(int(lineDict['DP']) - int(lineDict['AD']))
            record['t_alt_count'] = lineDict['AD']
            record['Hugo_Symbol'] = lineDict['gene']
            record['NCBI_Build'] = lineDict['reference']  # this needs to be added to the mutation file for confirmatory purposes.
            record['Variant_Classification'] = lineDict['consequence']

            # ** write yo output **
            lineOut = []
            for field in outHeader:
                lineOut.append(record[field])
            out1.write('\t'.join(lineOut) + '\n')
            # **

# ***    

# *** write svs file ***
def writeSvs(inFile,outDir):
    
    outFile = outDir + 'data_sv.txt'
    df1 = pd.read_csv(inFile,sep="\t")

    with open(outFile,'w') as out1:

        # write yo header
        header = ['Sample_Id','SV_Status','Site1_Hugo_Symbol','Site1_Chromosome','Site1_Position','Site2_Hugo_Symbol','Site2_Chromosome','Site2_Position','NCBI_Build','Class','SV_Length','Tumor_Variant_Count']
        out1.write('\t'.join(header) + '\n')

        # for each gene right a line
        for idx,row in df1.iterrows():
            lineDict = row.to_dict()

            record = {'SV_Status':'SOMATIC','NCBI_Build':'GRCh38'}
            record['Sample_Id'] = lineDict['sample']
            record['Site1_Hugo_Symbol'],record['Site2_Hugo_Symbol'] = lineDict['gene1'],lineDict['gene2']
            record['Site1_Chromosome'],record['Site1_Position'] = lineDict['chrom1'],lineDict['pos1']
            record['Site2_Chromosome'],record['Site2_Position'] = lineDict['chrom2'],lineDict['pos2']
            record['Class'] = lineDict['eventType']
            record['SV_Length'],record['Tumor_Variant_Count'] = lineDict['svLength'],lineDict['AD']

            lineOut = []
            for field in header:
                lineOut.append(record[field])
            lineOut = [str(field) for field in lineOut]
            out1.write('\t'.join(lineOut) + '\n')
# ***    

# *** write cnvs file ***
def writeCnvs(inFile,outDir):
    
    outFile = outDir + 'data_cna.txt'
    df1 = pd.read_csv(inFile,sep="\t")

    with open(outFile,'w') as out1:

        # write yo header
        samples = sorted(set(list(df1['sample'])))
        header = ['Hugo_Symbol'] + samples
        header = [str(field) for field in header]        
        out1.write('\t'.join(header) + '\n')

        # for each gene right a line
        for gene,df1a in df1.groupby('gene'):

            # order by samples
            df1b = df1a.set_index('sample')
            df1b = df1b.loc[samples]
            values = list(df1b['gisticValue'])
            
            lineOut = [gene] + values
            lineOut = [str(field) for field in lineOut]
            out1.write('\t'.join(lineOut) + '\n')
# ***    

# *** write meta files ***
def writeMetas(inFile,outDir):

    df1 = pd.read_csv(inFile,sep="\t")
    
    groups = sorted(list(set(df1['group'].tolist())))

    for group in groups:
        # write meta_study.txt
        if group == 'study':
            outFile = outDir + 'meta_study.txt'            
            writeMeta(group,df1,outFile)
                    
        # write meta_clinical_sample.txt
        elif group == 'clinicalSample':
            outFile = outDir + 'meta_clinical_sample.txt'            
            writeMeta(group,df1,outFile)
            
        # write meta_clinical_sample.txt
        elif group == 'clinicalPatient':
            outFile = outDir + 'meta_clinical_patient.txt'            
            writeMeta(group,df1,outFile)
            
        # write meta_clinical_sample.txt
        elif group == 'mutations':
            outFile = outDir + 'meta_mutations.txt'
            writeMeta(group,df1,outFile)        

        # write meta_clinical_sample.txt
        elif group == 'cnvs':
            outFile = outDir + 'meta_cna.txt'
            writeMeta(group,df1,outFile)
        elif group == 'svs':
            outFile = outDir + 'meta_sv.txt'
            writeMeta(group,df1,outFile)
        elif group == 'panel':
            outFile = outDir + 'meta_gene_panel_matrix.txt'
            writeMeta(group,df1,outFile)        
        else:
            raise 'i no understand'
# ***
            
# *** write metafile ***
def writeMeta(group,df1,file):
    with open(file,'w') as out1:
        for idx,row in df1.iterrows():
            lineDict = row.to_dict()

            if lineDict['group'] == group:
                lineOut = lineDict['parameter'] + ': ' + lineDict['value']
                out1.write(lineOut + '\n')
# ***
