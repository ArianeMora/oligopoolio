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
a.	5’: gaaataattttgtttaactttaagaaggagatatacat   
b.	3’: CTCGAGCACCACCACCACCACCACTGAGATCCGGCTGCTAACAAAGC  (i.e. rev comp of 5’ to 3’ region of the end) Histag + 006 primer commonly used in the Arnold lab   


## Step 1: 
### Simple pools:
Here you have a simple pool with only sequences < 350nt.
When you have short sequneces you just need to ensure that you create the gene sequnence with an overhanmg to the backbone. Then you can use universal primers to amplify the pool (see Generic primers above).

### Chimera pools:
Here you have pools with sequences between 350nt and 700nt (e.g. the protoglobins.) In this case the pool will have two sequences for each one, and we need an "overhang" to join 
the two parts of the sequence together. So essentially, you want 1) an overhang with the backbone, 2) an overlap with the other sequence, 3) you need to order primers for both the 5-->3 and 3--> for the overhang 
to ensure it amplifies correctly.

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
All the primer design stuff here is ai


### Debugging bench stuff
This is all use as is, if you find that some of the sequneces didn't amplify, the best thing you can do is to try debugging whether it was the designs or ya messed up a wetlab step.

Some simple things to try:  
-	PCR amplify specifically the shortest one and the middle one and the longest one  
-	Clone them in individually and check manually to ensure that there is not a design issue   
-	Be aware of ones which have a secondary structure hairpin  


