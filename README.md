# bioinf
Misc. helpful scripts for bioinformatics. May require installation of BioPython, HMMER3, and/or MUSCLE.

- drawOrf.py - basic ORF diagram drawer, takes a nucleotide GI or genbank file as input. NOT VERY USEFUL YET. PROBABLY DON'T WANT TO USE.
- phmmerall.py - tool to take a list of protein GIs and perform all-vs-all phmmer, without output directly importable into Cytoscape as a sequence similarity network
- auto_hmm.py - takes a list of protein GIs, downloads them, aligns with MUSCLE, creats an HMM, and optionally presses it for hmmscan
- hmmerformat.py	- format tblout output from hmmscan in rich HTML
- hmmerformatdom.py	- format a domtbl output from hmmscan in rich HTML
- phmmerall.py	- takes a list of GIs and does an all-vs-all phmmer, then generates an SSN
- phmmerallf.py - takes a FASTA file and does an all-vs-all phmmer, then generates an SSN
