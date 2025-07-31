#!/usr/bin/python
# combine together tables

# read and combine files
def merge(files,outFile):

    # write yo file
    with open(outFile,'w') as out1:
        # write header
        file = files[0]
        with open(file) as in1:
            header = in1.readline()
            out1.write(header)

        for file in files:
            with open(file) as in1:
                in1.readline()
                for line in in1:
                    out1.write(line)
                

            
