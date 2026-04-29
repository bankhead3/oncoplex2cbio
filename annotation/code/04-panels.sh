#!/bin/bash
# create panel files

inFile1=../input/opx-v8-genes-web.txt
outFile1=../opx-v8-gene-panel-web.txt

echo "stable_id: OPXv8" > $outFile1
echo "description: OncoPlex V8 gene list." >> $outFile1
echo -n -e "gene_list:\t" >> $outFile1
tail -n+2 $inFile1 | awk '{printf "%s\t", $0} END {print ""}' | sed 's/\t$//' >> $outFile1

inFile1=../input/opx-v8-genes.txt
outFile1=../opx-v8-gene-panel.txt
echo "stable_id: OPXv8" > $outFile1
echo "description: OncoPlex V8 gene list." >> $outFile1
echo -n -e "gene_list:\t" >> $outFile1
cat $inFile1 | sed 's/["]//g' | sed 's/[ ]//g' | sed 's/[,]/\n/g' | sort | uniq | sed '/^$/d' | awk '{printf "%s\t", $0} END {print ""}' | sed 's/\t$//' >> $outFile1
echo $outFile1

inFile1=../input/opx-v7-genes.txt
outFile1=../opx-v7-gene-panel.txt
echo "stable_id: OPXv7" > $outFile1
echo "description: OncoPlex V7 gene list." >> $outFile1
echo -n -e "gene_list:\t" >> $outFile1
cat $inFile1 | sed 's/["]//g' | sed 's/[ ]//g' | sed 's/[,]/\n/g' | sort | uniq | sed '/^$/d' | awk '{printf "%s\t", $0} END {print ""}' | sed 's/\t$//' >> $outFile1
echo $outFile1

inFile1=../input/opx-v6-genes.txt
outFile1=../opx-v6-gene-panel.txt
echo "stable_id: OPXv6" > $outFile1
echo "description: OncoPlex V6 gene list." >> $outFile1
echo -n -e "gene_list:\t" >> $outFile1
cat $inFile1 | sed 's/["]//g' | sed 's/[ *]//g' | sed 's/[,]/\t/g' | xargs -n1 -d'\t' | sort | uniq | sed '/^$/d' | awk '{printf "%s\t", $0} END {print ""}' | sed 's/\t$//' >> $outFile1
echo $outFile1

inFile1=../input/opx-v5-genes.txt
outFile1=../opx-v5-gene-panel.txt
echo "stable_id: OPXv5" > $outFile1
echo "description: OncoPlex V5 gene list." >> $outFile1
echo -n -e "gene_list:\t" >> $outFile1
cat $inFile1 | sed 's/["]//g' | sed 's/[ *]//g' | sed 's/[,]/\t/g' | xargs -n1 -d'\t' | sort | uniq | sed '/^$/d' | awk '{printf "%s\t", $0} END {print ""}' | sed 's/\t$//' >> $outFile1
echo $outFile1
