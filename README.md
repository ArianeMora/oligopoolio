# Code and protocols for making oligopools for ordering proteins

### Setting up the code

```
conda create --name oligo python=3.10
```

```
conda activate oligo
pip install -r requirements.txt 
```

### Dependencies 
None outside python for the primer generation.  

For the demultiplexing and annotating reads you'll need the following tools added to your path.  
1. clustal-omega  
2. samtools  
3. minimap2  

## Generic primers order:

Generic primers:  
a.	5’: `gaaataattttgtttaactttaagaaggagatatacat`
b.	3’: `CTCGAGCACCACCACCACCACCACTGAGATCCGGCTGCTAACAAAGC`  (i.e. rev comp of 5’ to 3’ region of the end) Histag + 006 primer commonly used in the Arnold lab   


## Step 1: 

Codon optimize your gene sequences for whichever organism you are planning on expressing with. 

I do this using IDT. You can bulk upload your sequences and then download them as a fasta file. **This is your input to this tool!**

#### Getting primers
To generate a primer for a single sequence you would do:

```
# Placeholder gene ends (replace with your actual gene sequences)
gene = "ATGAGCGATCTGCATAACGAGTCCATTTTTATTACCGGCGGCGGATCGGGATTAGGGCTGGCGCTGGTCGAGCGATTTAT\
CGAAGAAGGCGCGCAGGTTGCCACGCTGGAACTGTCGGCGGCAAAAGTCGCCAGTCTGCGTCAGCGATTTGGCGAACATA\
AATAATCGCGCATTAAGCGGTGTGATGATCAACGCTGATGCGGGTTTAGCGATTCGCGGCATTCGCCACGTAGCGGCTGG\
GCTGGATCTTTAA"

# Standard pET-22b(+) primer sequences
forward_plasmid_primer = "GGAGATATACATATG"
reverse_plasmid_primer = "GCTTTGTTAGCAGCCGGATCTCA"

# Desired Tm range for optimization
desired_tm = 62.0  # Target melting temperature in °C
tm_tolerance = 5.0  # Allowable deviation from the desired Tm

# Generate and optimize forward primer
min_length = 13
max_length = 20
forward_gene_primer, forward_tm = optimize_primer(forward_plasmid_primer, gene, desired_tm, 'forward',
                                                  min_length, max_length, tm_tolerance)
reverse_gene_primer, reverse_tm = optimize_primer(reverse_plasmid_primer, gene, desired_tm, 'reverse',
                                                  min_length, max_length, tm_tolerance)
```

To generate a primer for a fasta file of sequences you would do the following (this will output it for IDT):

```
make_primers_IDT(fasta_file, remove_stop_codon=True, his_tag='',
                     max_length=60, min_length=15, tm_tolerance=30, desired_tm=62.0,
                     forward_primer='gaaataattttgtttaactttaagaaggagatatacat',
                     reverse_primer='ctttgttagcagccggatc')
```


#### Get flanking primers e.g. you want to PCR a gene out of E.coli

You need the GFF file, you can download this from NCBI, and the gene sequence (here you don't need to 
do the optimization from IDT since it is already in the nucleotide sequence).
```
gff_file = f'{test_data_dir}genome_NEB_B/genomic.gff'
reference_fasta = f'{test_data_dir}genome_NEB_B/GCF_001559615.2_ASM155961v2_genomic.fna'
gene_name = 'ysgA'
seqid, start, end, strand, upstream_flank, downstream_flank, gene_seq = get_flanking_primers(gene_name,
                                                                                             gff_file,
                                                                                             reference_fasta)
```

### Simple pools:
Here you have a simple pool with only DNA sequences < 350nt.

When you have short sequences you just need to ensure that you create the gene sequences with an overhang to the backbone. 
Then you can use universal primers to amplify the pool (see Generic primers above).


command:
```

```

### Chimera pools:
Here you have pools with sequences between 350nt and 700nt (e.g. the protoglobins.) In this case the pool will have two sequences for each one, and we need an "overhang" to join 
the two parts of the sequence together. So essentially, you want 1) an overhang with the backbone, 2) an overlap with the other sequence, 3) you need to order primers for both the 5-->3 and 3--> for the overhang 
to ensure it amplifies correctly.

Again this is just your codon optimized DNA sequences for your genes you want to order, note they should be less than 640 base pairs! 

command:
```

```
##### Chimera pool explanation:

Part 1 is the start of the gene:  
-	15bp upstream + gene to 50% of gene  
-	Forward primer of the start: gaaataattttgtttaactttaagaaggagatatacat  


Part 2 is the end of the gene:  
-	3’ PCR reverse complement: gcagccaactcagcttcctttcgggctttgttagcagccggatc  
-	5’3’ 15bp overlap + the Part1 of the gene + last 50% of sequence + overlap with 3’ primer: e.g. end of this is CCCAATCCACGTCTTgatccggctgctaac  


#### Example of a chimera
Gaaggagatatacat = overlap with backbone  
Gcagcgtgttcgtcgttt = overlap between the two oligos  
Gatccggctgctaac = Overlap with the 3’ primer backbone  

**Oligo for part 1:**  
gaaggagatatacatATGGACGACCTGGAACGTGCAGGCAAAGATGCGTGGACATTTGAAAAGGCATTAGCGCGCCTGGAAGAAGTAGTAGAACGTCTGGAGAGTGCAGACCTGCCATTGGATAAGGCATTAAGTCTTTACGAGGAGGGCACCCGCCTTGTTCGTTATCTGAACGGTGAATTGAATCGTTTTGAgcagcgtgttcgtcgttt

**Oligo for part 2:**  
gcagcgtgttcgtcgtttGCGCGAAGAGGAGGTATCCCCGGAACCTAAAGTCAGTGAGGGGTTTGCTCCCGCGTCAGAAAATGAGTTGTTTCCCTTCGAGGGAGAGGAAGATTTCGCGGAGTGGGAGGATGAAATCGAATTTGAGGAGGAGTTCCCCGGCGAAGAGGAAGAGGGTGATGATCCCAATCCACGTCTTgatccggctgctaac

**Middle of the gene primers:**  
-	Overhang 5’--> 3’: GCAGCGTGTTCGTCGTTT    
-	Overhang reverse complement 3’-->5’: AAACGACGAACACGCTGC    


### Primer design
All the primer design stuff here is just appended so you should have checked your primers before. I do check the primers using primer3 and provide feedback about the melting temp.


### Debugging bench stuff
This is all use as is, if you find that some of the sequneces didn't amplify, the best thing you can do is to try debugging whether it was the designs or ya messed up a wetlab step.

Some simple things to try:  
-	PCR amplify specifically the shortest one and the middle one and the longest one  
-	Clone them in individually and check manually to ensure that there is not a design issue   
-	Be aware of ones which have a secondary structure hairpin  


