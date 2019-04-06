#!/bin/bash
# phylo tree building script
# place maketree.sh in directory of .txt 16s sequencing files and run

echo "======
1)file formatting
2)mafft alignment
3)fasttree build
======"

# concatenate the 16s files into 1 text file fasta format with newlines

sed -e '$s/$/\n/' -s *.txt > all16s.txt

# trimming file to remove from ".abi" to number of nt reads

sed -e 's/\..*\d*//' all16s.txt > trim16s.txt
#rm all16s.txt

# run mafft, --auto -nuc

mafft --nuc --quiet --auto trim16s.txt > rrna_align
#rm trim16s.txt

# run fasttree

mv rrna_align ~/Desktop/fasttree/rrna_align
cd ~/Desktop/fasttree
./FastTree -quiet -nt rrna_align > rrna_tree
#rm rrna_align

echo "EOF"

#DEBUG: IF SCRIPT IS RUN WITH FILES MATCHING THE NAMES ALREADY THERE, IT WILL HANG, ONLY REMOVE THE RM COMMANDS FOR DEBUGGING OR SPECIAL CASES
