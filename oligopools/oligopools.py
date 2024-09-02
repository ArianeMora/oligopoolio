from Bio import SeqIO
import pandas as pd
from sciutil import SciUtil
from Bio.Seq import Seq
import os
import pandas as pd
from Bio import SeqIO
import pysam
import numpy as np
from collections import defaultdict

u = SciUtil()

def make_oligo_single(codon_optimized_fasta, forward_primer='gaaataattttgtttaactttaagaaggagatatacat', 
                      forward_primer_len=15, reverse_primer='gatccggctgctaacaaag', reverse_primer_len=15, max_len=320):
    records = list(SeqIO.parse(codon_optimized_fasta, "fasta"))
    rows = []
    for r in records:
        seq = str(r.seq)
        if seq[:3] == 'ATG' and len(seq) < max_len:
            rows.append([str(r.id), seq])
        else:
            u.dp([r.id, f'Issue: either did not have a methionine start or was longer than {max_len}. This has been omitted'])
                    
    # Now we want to drop any that are > 320 lengths 
    single_df = pd.DataFrame(rows)
    single_df.columns = ['id', 'seq']
    single_df['forward_primer'] = forward_primer[-1*forward_primer_len:]
    single_df['reverse_primer'] = reverse_primer[:reverse_primer_len]
    single_df['oligo'] = [f'{f}{seq}{r}' for f, seq, r in single_df[['forward_primer', 'seq', 'reverse_primer']].values]
    return single_df

def make_oligo_double(codon_optimized_fasta, forward_primer='gaaataattttgtttaactttaagaaggagatatacat', 
                      forward_primer_len=15, reverse_primer='gatccggctgctaacaaag', reverse_primer_len=15, max_len=640, 
                      overlap_len=9):

    records = list(SeqIO.parse(codon_optimized_fasta, "fasta"))
    rows = []
    fwd = forward_primer[-1*forward_primer_len:]
    bwd = reverse_primer[:reverse_primer_len]
    for seq_id, seq in records.values:
        if seq[:3] == 'ATG' and len(seq) < max_len:
            cut = int(len(seq)/2)
            first_half = seq[:cut-overlap_len]
            second_half = seq[-cut+overlap_len:]
            overlap = seq[-cut-overlap_len:-cut+overlap_len].lower()
            # Also reverse complement the overlap
            rev_comp_overlap = str(Seq(overlap).reverse_complement())
            rows.append([seq_id + '_1', seq, 'first', overlap, rev_comp_overlap, f'{fwd}{first_half}{overlap}'])
            rows.append([seq_id + '_2', seq, 'second', overlap, rev_comp_overlap, f'{overlap}{second_half}{bwd}'])
        else:
            u.dp([seq_id, f'Issue: either did not have a methionine start or was longer than {max_len}. This has been omitted'])
                    
    # Now we want to drop any that are > 320 lengths 
    double_df = pd.DataFrame(rows)
    double_df.columns = ['id', 'seq', 'label', 'overlap', 'rev_comp_overlap', 'oligo']
    return double_df


# Take the demultiplexed files from a LevSeq run and then demultiplex oligoPools from these

def alignment_from_cigar(cigar: str, alignment: str, ref: str, query_qualities: list):
    """
    Generate the alignment from the cigar string.
    Operation	Description	Consumes query	Consumes reference
    0 M	alignment match (can be a sequence match or mismatch)	yes	yes
    1 I	insertion to the reference	yes	no
    2 D	deletion from the reference	no	yes
    3 N	skipped region from the reference	no	yes
    4 S	soft clipping (clipped sequences present in SEQ)	yes	no
    5 H	hard clipping (clipped sequences NOT present in SEQ)	no	no
    6 P	padding (silent deletion from padded reference)	no	no
    7 =	sequence match	yes	yes
    8 X	sequence mismatch	yes	yes
    """
    new_seq = ''
    ref_seq = ''
    qual = []
    inserts = []
    pos = 0
    ref_pos = 0
    for op, op_len in cigar:
        if op == 0:  # alignment match (can be a sequence match or mismatch)
            new_seq += alignment[pos:pos + op_len]
            qual += query_qualities[pos:pos + op_len]

            ref_seq += ref[ref_pos:ref_pos + op_len]
            pos += op_len
            ref_pos += op_len
        elif op == 1:  # insertion to the reference
            inserts.append(alignment[pos - 1:pos + op_len])
            pos += op_len
        elif op == 2:  # deletion from the reference
            new_seq += '-' * op_len
            qual += [-1] * op_len
            ref_seq += ref[ref_pos:ref_pos + op_len]
            ref_pos += op_len
        elif op == 3:  # skipped region from the reference
            new_seq += '*' * op_len
            qual += [-2] * op_len
            ref_pos += op_len
        elif op == 4:  # soft clipping (clipped sequences present in SEQ)
            inserts.append(alignment[pos:pos + op_len])
            pos += op_len
        elif op == 5:  # hard clipping (clipped sequences NOT present in SEQ)
            continue
        elif op == 6:  # padding (silent deletion from padded reference)
            continue
        elif op == 7:  # sequence mismatch
            new_seq += alignment[pos:pos + op_len]
            ref_seq += ref[ref_pos:ref_pos + op_len]
            qual += query_qualities[pos:pos + op_len]
            pos += op_len
            ref_pos += op_len
    return new_seq, ref_seq, qual, inserts


def annotate_to_wells(plate_path, ref_fasta):
    """ 
    This is for annotating the demuiltiplexed files from LevSeq to the wells that it came from 
    Assuming that you've demultiplexed the long read reads from a 96 well plate using the levSeq formatting.
    """
    files = os.listdir(plate_path)
    # Then just go through each well and annotate these using minimap2
    # make a mapping between the read id and the sequence so that this can be used as the reference string when aligning
    ref_map = {}
    seqs = [str(record.seq) for record in SeqIO.parse(ref_fasta, "fasta")]
    ref_read_ids = [str(record.id) for record in SeqIO.parse(ref_fasta, "fasta")]
    for i, seq in enumerate(seqs):
        ref_map[ref_read_ids[i]] = seq
    rows_all = []
    for well in files:
        try:
            msa_path = f'{plate_path}/{well}/{well}.msa'
            # get the file from there
            os.system(f'minimap2 -ax map-ont -A 4 -B 2 -O 10,24 {ref_fasta} {plate_path}/{well}/demultiplexed_RB07_{well}_000.fastq > {plate_path}/{well}/{well}.sam')
            os.system(f'samtools view -bS {plate_path}/{well}/{well}.sam > {plate_path}/{well}/{well}.bam')
            os.system(f'samtools sort {plate_path}/{well}/{well}.bam -o {plate_path}/{well}/{well}.bam')
            os.system(f'samtools index {plate_path}/{well}/{well}.bam')
            bam_file_path = f'{plate_path}/{well}/{well}.bam'
            bam = pysam.AlignmentFile(bam_file_path, "rb")
            # Ensure the BAM file is indexed
            if not os.path.exists(bam_file_path + ".bai"):
                pysam.index(bam_file_path)

            cramHeader = bam.header.to_dict()
            
            seqs = []
            read_ids = []
            read_quals = []
            err = 0
            read_counts = defaultdict(int)

            for read in bam.fetch(until_eof=True):
                try:
                    ref_str = ref_map.get(read.reference_name)
                    if read.query_sequence is not None and len(read.query_sequence) > 0.9*len(ref_str) and read.cigartuples is not None:
                        seq, ref, qual, ins = alignment_from_cigar(read.cigartuples, read.query_sequence, ref_str,
                                                                read.query_qualities)
                        # Make it totally align
                        seq = "-" * read.reference_start + seq + "-" * (len(ref_str) - (read.reference_start + len(seq)))
                        seqs.append(seq)
                        #seqs.append(read.query_sequence)
                        read_ids.append(f'{read.reference_name}')
                        read_quals.append(read.qual)
                        # keep track of how many reads there were
                        read_counts[read.reference_name] += 1
                except:
                    err += 1
            print(err)
            # Check if we want to write a MSA
            if msa_path is not None:
                print("Writing MSA")
                with open(msa_path, 'w+') as fout:
                    # Write the reference first
                    fout.write(f'>{read.query_name}\n{ref_str}\n')
                    for i, seq in enumerate(seqs):
                        fout.write(f'>{read_ids[i]}\n{"".join(seq)}\n')
                # Align using clustal for debugging if you need the adapter! Here you would change above to use a different version
                print(f'clustal-omega --force -i {msa_path} -o {msa_path.replace(".fa", "_msa.fa")}')
                os.system(f'clustal-omega --force -i "{msa_path}" -o "{msa_path.replace(".fa", "_msa.fa")}"')
                seqs = [str(record.seq) for record in SeqIO.parse(msa_path.replace(".fa", "_msa.fa"), "fasta")]
                read_ids = [str(record.id) for record in SeqIO.parse(msa_path.replace(".fa", "_msa.fa"), "fasta")]
        except:
            print(well)
        #Again check that we actually had enough reads for this to be considered a good well
        # Count how many unique read IDs mapped to this (we hope that it is all the same )
        bam.close()
        # Get the number for each read in there
        row = []
        most_read = 0
        well_id = None
        for read in ref_read_ids:
            row.append(read_counts.get(read))
            if read_counts.get(read) and read_counts.get(read) > most_read:
                most_read = read_counts.get(read)
                well_id = read
        rows_all.append([well, len(read_ids), well_id, most_read] + row) #, np.mean(read_quals)])

    if len(rows_all) > 1:  # Check if we have anything to return
        seq_df = pd.DataFrame(rows_all, columns=['Well', 'Number of reads', 'Assigned ID', 'Reads for Assigned ID'] + ref_read_ids)
        seq_df.to_csv('cutadapt_demultiplex_seqs.csv', index=False)



def optimize_primer(plasmid_sequence, gene_sequence, desired_tm, direction, min_length=10, max_length=60,
                    tm_tolerance=5, max_temp_deviation = 20, his_tag=False):
    """Optimize a primer sequence to achieve a desired Tm by adjusting its length."""

    best_primer = 'None'
    temp = 0
    if direction == 'forward':
        gene_sequence = gene_sequence
    elif direction == 'reverse':
        # Reverse comp it
        gene_sequence = str(Seq(gene_sequence).reverse_complement())
    else:
        u.dp([f'Direction: {direction} is not an option, must be: forward or reverse.'])
    if his_tag and direction == 'reverse':
        min_length = 33  # i.e. the 18 of his tag + the seq
    for primer_length in range(0, 15, 3):
        for length in range(min_length, max_length, 3):
            plas_seq = plasmid_sequence[-1*(primer_length + 9):]
            temp_primer = plas_seq + gene_sequence[:length]
            tm = primer3.bindings.calcTm(temp_primer)
            temp_deviation = abs(desired_tm - tm)
            if temp_deviation < max_temp_deviation and temp_deviation < tm_tolerance and (temp_primer[-1] == 'C' or temp_primer[-1] == 'G'):
                max_temp_deviation = temp_deviation
                best_primer = plas_seq.lower() + gene_sequence[:length]
                temp = tm
    return best_primer, temp

def make_primer(gene, max_length=60, min_length=15, tm_tolerance=30, desired_tm=62.0, 
                forward_primer='gaaataattttgtttaactttaagaaggagatatacat', reverse_primer='ctttgttagcagccggatc'):
    # Standard pET-22b(+) primer sequences
    forward_plasmid_primer = forward_primer.upper()  # Cut off a bit and then do the same for the tail
    # Play with the length to get the melting but you can optimise for around 38 length --> Full GCTTTGTTAGCAG  --> CCGGATCTCA
    reverse_plasmid_primer = reverse_primer.upper()
    # Desired Tm range for optimization
    # Target melting temperature in °C
    # Generate and optimize forward primer
    forward_gene_primer, forward_tm = optimize_primer(forward_plasmid_primer, gene, desired_tm, 'forward',
                                                      min_length, max_length, tm_tolerance)
    reverse_gene_primer, reverse_tm = optimize_primer(reverse_plasmid_primer, gene, desired_tm, 'reverse',
                                                      min_length, max_length, tm_tolerance, his_tag=True)

    print(f"Forward Gene Primer: 5'-{forward_gene_primer}-3' (Tm: {forward_tm:.2f} °C)")
    print(f"Reverse Gene Primer: 3'-{reverse_gene_primer}-5' (Tm: {reverse_tm:.2f} °C)")
    return forward_gene_primer, forward_tm, reverse_gene_primer, reverse_tm 

def make_primers_IDT(fasta_file):
    seqs = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        try:
            seq_id = str(record.id)
            seqs[seq_id] = str(record.seq)
        except:
            print(record.description)
    rows = []
    for gene in seqs:
        gene_seq = seqs.get(gene)
        if gene_seq:
            forward_gene_primer, forward_tm, reverse_gene_primer, reverse_tm = make_primer(gene_seq[:-3] + 'CACCACCACCACCACCAC') # Cut off the stop codon and add in the his Tag
            rows.append([f'5F_{gene}', forward_gene_primer, '25nm', 'STD'])
            rows.append([f'3R_{gene}', reverse_gene_primer, '25nm', 'STD'])
            print(len(forward_gene_primer), len(reverse_gene_primer))
        else:
            print(gene)
    primer_df = pd.DataFrame(rows)
    primer_df.columns = ['Name', 'Sequence', 'Scale', 'Purification']
    return primer_df