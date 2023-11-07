# ProteinSequencing_CLI

In this project, I analyze DNA sequences for the gene p53, which is used to suppress cancer in organisms. Specifically, I compare the p53 genes in humans and elephants, to identify what they have in common and how they are different.

To do this, I interpret DNA data files from the NIH, convert them to RNA, then convert that to a sequence of generated proteins. Next, I compare the protein sequences of the two different genes, automatically generate text reports of your findings, and graph the results, which are displayed at the bottom of this notebook.

The data that is used for this project is retrieved from the NIH Gene Database.

Human p53 Gene:
https://www.ncbi.nlm.nih.gov/nuccore/NC_000017.11?report=fasta&from=7668402&to=7687550&strand=true

Elephant p53 Gene:
https://www.ncbi.nlm.nih.gov/nuccore/NW_003573467.1?report=fasta&from=11687413&to=11699835&strand=true

DNA is composed of four bases: cytosine (C), guanine (G), adenine (A), and thymine (T). These bases, when read together, produce instructions that the organism can follow in order to create useful proteins. The connected sequence of bases is represented in a text file with a single letter per base.

This program is designed to work on the command line and can be run with one command. The two DNA sequences are stored in the same directory as the prot_seq.py file and the command to run is simply:

python prot_seq {file_path_1} {file_path_2}