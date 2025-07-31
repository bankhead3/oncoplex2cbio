#!/usr/bin/env python3
# parser for metadata recieved 

import os

def parse(files,outFile):
    
    # create a table with the following
    # patient, sample, field 
    with open(outFile,'w') as out1:
        outHeader = ['id','patientSample','property','value','description','name','patientID']
        out1.write('\t'.join(outHeader) + '\n')
    
    records = {}
    for file in files:
        basename = os.path.basename(file)

        if basename == 'data_clinical_patient.txt':
            parseClinicalPatient(file,outFile,outHeader,'patient')
        elif basename == 'data_clinical_sample.txt':
            parseClinicalPatient(file,outFile,outHeader,'sample')            
        else:
            print('i no unerstand file')
            raise
            
        

# parse data_clinical_patient.txt data
def parseClinicalPatient(inFile,outFile,outHeader,patientSample):

    with open(inFile) as in1,open(outFile,'a') as out1:
        names = in1.readline()
        names = names.strip().split('\t')
        descs = in1.readline()
        descs = descs.strip().split('\t')
        in1.readline(); in1.readline();
        header = in1.readline()
        header = header.strip().split('\t')

        assert len(descs) == len(header) and len(names) == len(header)
        nameDict = dict(zip(header,names))        
        descDict = dict(zip(header,descs))
        
        for line in in1:
            parse1 = line.strip().split('\t')
            assert len(parse1) == len(header)
            lineDict = dict(zip(header,parse1))

            if patientSample == 'patient':
                id = lineDict['PATIENT_ID']
            elif patientSample == 'sample':
                id = lineDict['SAMPLE_ID']
            else:
                print('i no understand')
                raise 
                
            id = id.replace('_ptID','')
            patientID = 'patient' + id
            id = patientSample + id
            properties = list(lineDict.keys())
            properties = [field for field in properties if field != 'PATIENT_ID']

            for property in properties:
                lineDict['property'] = property
                lineDict['value'] = lineDict[property]
                lineDict['id'] = id
                lineDict['name'] = nameDict[property]
                lineDict['description'] = descDict[property]                
                lineDict['patientSample'] = patientSample
                lineDict['patientID'] = patientID
                lineOut = []
                for field in outHeader:
                    lineOut.append(lineDict[field])
                out1.write('\t'.join(lineOut) + '\n')
