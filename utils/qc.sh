#!/bin/bash

# make sure we have 3 total arguments
if [ $# != 2 ]; then echo "USAGE: $0 <file> <uniqColumns>"; exit; fi 

inFile=$1
uniqColumns=$2

echo -n 'QC '

# make sure file exists
if [ ! -e $inFile ]; then echo "ERROR: FILE DOES NOT EXIST! $inFile"; exit; fi

grep -e '[0-9]' <(echo $uniqColumns) > /dev/null
if [ $? == 1 ]; then echo "ERROR: INVALID COLUMN #'S SPECIFIED! $uniqColumns"; exit; fi

# test 1: general dups
numDups=$(sort $inFile | sort | uniq -d | wc -l)
if [ $numDups != 0 ]; then echo "FAIL: GENERAL DUPLICATES PRESENT IN $inFile !!"; exit; fi

# test 2: uniq fields are uniq
numDups=$(cut -f$uniqColumns $inFile | sort | uniq -d | wc -l)
if [ $numDups != 0 ]; then echo "FAIL: DUPLICATES PRESENT IN $inFile col $uniqColumns!!"; exit; fi

# test 3: consistent # of fields
numErrors=$(awk 'BEGIN{OFS=FS="\t"} NR == 1 { headerNF=NF; header=$0 } headerNF != NF { print "ERROR!"}' $inFile | wc -l)
if [ $numErrors != 0 ]; then
    awk 'BEGIN{OFS=FS="\t"} NR == 1 { headerNF=NF; header=$0 } headerNF != NF { print "FAIL: NUMBER OF FIELDS DOES NOT MATCH HEADER!!"; print header; print $0; exit }' $inFile
  exit; 
fi

# test 4: check for blank fields
numErrors=$(awk 'BEGIN{OFS=FS="\t"} { for(i = 1; i <= NF; i++) if(length($i) == 0) print "ERROR!"}' $inFile | wc -l)
if [ $numErrors != 0 ]; then
  awk 'BEGIN{OFS=FS="\t"} { for(i = 1; i <= NF; i++) if(length($i) == 0) { print "FAIL: FIELD",i," IS BLANK!!"; print $0; exit } }' $inFile
  exit;
fi

lines=$(wc -l $inFile | awk '{ print $1 }')
echo "PASS: $inFile $lines lines..."
