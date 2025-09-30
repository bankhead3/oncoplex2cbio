# oncoplex2cbio
Converts OncoPlex variants to cBioPortal format.

## Important Scripts:
1. parse_variants.py - converts oncoplex workbooks to table delimited files (supports v8 oncoplex only)
    1. parse(inFile, outDir, sample) - takes as input oncoplex excel workbook and generates tab-delimited variant files for mutations, cnvs, svs
2. process_variants.py - filters variants and re-formats tab-delimited variant files
    1. process(inFile, outFile) - takes as input a mutation, cnv, or sv tab-delimited variant file and outputs an updated file
3. build_cbio.py - writes mutation, metadata, cnv, sv, panel, and case list cbioportal files
    1. build(fileDict,outDir) - takes as input a python dictionary specifying file locations for mutation, metadata, cnv, sv, panel, and caseLists 

## Requirements:
1. oncoplex v8 excel workbooks
2. python3
