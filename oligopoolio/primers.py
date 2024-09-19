###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

"""
Author: Ariane Mora
Date: 24th March 2024
"""

from Bio.Seq import Seq
import primer3
from sciutil import SciUtil
import gffutils
import pandas as pd
from Bio import SeqIO
import pandas as pd

u = SciUtil()


def reverse(seq):
    """ Just reverse a sequence. """
    return reversed(seq)


def check_aa(seq, allow_gaps=False):
    """ Check that the sequence starts with a methionine and ends with a stop codon."""
    bases = 'ACDEFGHIKLMNPQRSTVWY'
    if seq[0] != 'M':
        u.warn_p(['Your sequence did not start with Methionine...'])
        return False
    if 'M' in seq[1:]:
        u.warn_p(['Your sequence had multiple Methionines...'])
        return False
    for b in seq:
        if allow_gaps:
            bases = 'ACDEFGHIKLMNPQRSTVWY-'
        if b not in bases:
            u.warn_p(['Your sequence contained non-canonical bases...', b])
            return False
    return True


def check_nt(seq, allow_gaps = False):
    """ Check that the sequence starts with a methionine and ends with a stop codon."""
    bases = 'ATGC'
    if seq[:3] != 'ATG':
        u.warn_p(['Your sequence did not start with Methionine ATG...'])
        return False
    if seq[-3:] != 'TAA':
        u.warn_p(['Your sequence did not end with stop codon TAA...'])
        return False
    if 'TAA' in seq[:-3]:
        u.warn_p(['Your sequence had multiple stop codons...'])
        return False
    if 'ATG' in seq[3:]:
        u.warn_p(['Your sequence had multiple start codons...'])
        return False
    for b in seq:
        if allow_gaps:
            bases = 'ATGC-'
        if b not in bases:
            u.warn_p(['Your sequence contained non-canonical bases...', b])
            return False
    return True


def clean_aa_seq(seq, allow_gaps=True, remove=False, replace='N', amino_acids=None):
    """
        Remove any non-canonical bases and swap them with N. amino_acids can be
        specified if the user wants non-canonical amino acid bases.
    """
    if allow_gaps:
        bases = amino_acids or list('ACDEFGHIKLMNPQRSTVWY-')  # Allow gaps
    else:
        bases = amino_acids or list('ACDEFGHIKLMNPQRSTVWY')
    seq = seq.strip().replace(' ', '')  # Remove any spaces
    if remove:
        replace = ''  # Swap the other base of N for -
    return ''.join([s if s in bases else replace for s in seq])


def clean_nt_seq(seq, allow_gaps=True, remove=False, replace='N'):
    """ Remove any non-canonical bases and swap them with N."""
    if allow_gaps:
        bases = ['A', 'T', 'G', 'C', '-']  # Allow gaps
    else:
        bases = ['A', 'T', 'G', 'C']
    seq = seq.strip().replace(' ', '')  # Remove any spaces
    if remove:
        replace = ''  # Swap the other base of N for -
    return ''.join([s if s in bases else replace for s in seq])


def reverse_complement(seq):
    """Generate the reverse complement of a DNA sequence."""
    seq = seq.strip().replace(' ', '') # Do a little cleaning
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join([complement[base] for base in reversed(seq)])


def gff_to_dataframe(db):
    # Function to convert GFF database to pandas DataFrame
    records = []

    # Iterate over each feature in the database
    for feature in db.all_features():
        record = {
            'seqid': feature.seqid,
            'source': feature.source,
            'feature_type': feature.featuretype,
            'start': feature.start,
            'end': feature.end,
            'score': feature.score,
            'strand': feature.strand,
            'phase': feature.frame,
        }

        # Add attributes as separate columns
        for key, value in feature.attributes.items():
            record[key] = value[0] if len(value) == 1 else value
        # Let's also do this for each of the dbx references making it a new column
        db_entries = record.get('Dbxref')
        if db_entries:
            if isinstance(db_entries, str):
                record[db_entries.split(':')[0]] = db_entries.split(':')[1]
            elif isinstance(db_entries, list):
                for entry in db_entries:
                    record[entry.split(':')[0]] = entry.split(':')[1]

        records.append(record)

    return pd.DataFrame(records)


def create_db_from_gff(gff_file):
    # Function to create a database from the GFF file
    db_file = gff_file + '.db'
    return gffutils.create_db(gff_file, dbfn=db_file, force=True, keep_order=False, merge_strategy='merge',
                              sort_attribute_values=True)


def get_flanking_primers(gene_id, gff_file, fasta_reference, upstream_flank_number=50, downstream_flank_number=50):
    """
    Get the flanking primers for a given gene from a reference sequence.
https://benchling.com/arnold_lab/f/lib_yyMnf2lS-trpb_landscape/prt_kpSpRW0e-1-red-mediated-e-coli-knockout/edit
    Note we expect the file format to be in that from NCBI.
    """
    # First get the gene location from the gff file.
    db = create_db_from_gff(gff_file)

    # Convert the database to a pandas DataFrame
    gff_df = gff_to_dataframe(db)
    gff_df.to_csv('test.csv')
    # Then get the position in the fasta file and ensure that the upstream starts at TAG and ends with TAA
    # i.e. start and stop codon.
    gene_df = gff_df[gff_df['Name'] == gene_id]
    start = gene_df['start'].values[0]
    end = gene_df['end'].values[0]
    strand = gene_df['strand'].values[0]
    seqid = gene_df['seqid'].values[0]
    # Get this from the fasta file now
    # There is probably a nicer way to do this with the sequence package but given I'm on a plane gonna do this dumb
    # dumb way and I'm really tired
    upstream_flank = ''
    gene_seq = ''
    downstream_flank = ''
    with open(fasta_reference, 'r+') as fin:
        seq_found = False
        sequence = ''
        for line in fin:
            if line[0] == '>':
                if seqid in line[1:].strip():  # This is our sequence so we just go to that position
                    # Need to get the next line
                    seq_found = True
                elif seq_found:
                    # This is now our sequnece and it is complete
                    if strand == '+':  # This is normal direction
                        gene_seq = sequence[start - 1:end]
                        upstream_flank = sequence[start - 1 - upstream_flank_number: start - 1]
                        downstream_flank = sequence[end: end + downstream_flank_number + 1]
                        return seqid, start, end, strand, upstream_flank, downstream_flank, gene_seq

                    else:  # reverse stranded so the upstream is actually from the end
                        gene_seq = sequence[start - (upstream_flank_number + 1):end + downstream_flank_number]
                        # Reverse complement this
                        gene_seq = reverse_complement(gene_seq)
                        upstream_flank = gene_seq[:upstream_flank_number]
                        downstream_flank = gene_seq[-downstream_flank_number:]
                        gene_seq = gene_seq[upstream_flank_number: -downstream_flank_number]
                        return seqid, start, end, strand, upstream_flank, downstream_flank, gene_seq
            elif seq_found:
                sequence += line.strip() # Build the sequence.
    if seq_found:
        # This is now our sequnece and it is complete
        if strand == '+':  # This is normal direction
            gene_seq = sequence[start - 1:end]
            upstream_flank = sequence[start - 1 - upstream_flank_number: start - 1]
            downstream_flank = sequence[end: end + downstream_flank_number]
        else:  # reverse stranded so the upstream is actually from the end
            # Possibly also reverset complement this???
            gene_seq = sequence[start - (upstream_flank_number + 1):end + downstream_flank_number]
            # Reverse complement this
            gene_seq = reverse_complement(gene_seq)
            upstream_flank = gene_seq[:upstream_flank_number]
            downstream_flank = gene_seq[-downstream_flank_number:]
            gene_seq = gene_seq[upstream_flank_number: -downstream_flank_number]

    return seqid, start, end, strand, upstream_flank, downstream_flank, gene_seq


def optimize_primer(plasmid_sequence, gene_sequence, desired_tm, direction, min_length=10, max_length=60,
                    tm_tolerance=5, max_temp_deviation=20, his_tag=False):
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
        min_length = min_length + 18  # 33  # i.e. the 18 of his tag + the seq
    for primer_length in range(0, min_length, 3):
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


def make_primers_IDT(fasta_file, remove_stop_codon=True, his_tag='',
                     max_length=60, min_length=15, tm_tolerance=30, desired_tm=62.0,
                     forward_primer='gaaataattttgtttaactttaagaaggagatatacat',
                     reverse_primer='ctttgttagcagccggatc'):
    """
    optionally would set the his_tag to be: 'CACCACCACCACCACCAC'
    # Note expects the stop codon to still be there since this gets removed and then the his tag added
    """
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
        if remove_stop_codon:
            gene_seq = gene_seq[: -3]
        if his_tag:
            gene_seq = gene_seq + his_tag
        if gene_seq:
            forward_gene_primer, forward_tm, reverse_gene_primer, reverse_tm = make_primer(gene_seq,
                                                                                           max_length=max_length,
                                                                                           min_length=min_length,
                                                                                           tm_toleranc=tm_tolerance,
                                                                                           desired_tm=desired_tm,
                                                                                           forward_primer=forward_primer,
                                                                                           reverse_primer=reverse_primer)
            # Cut off the stop codon and add in the his Tag
            rows.append([f'5F_{gene}', forward_gene_primer, '25nm', 'STD'])
            rows.append([f'3R_{gene}', reverse_gene_primer, '25nm', 'STD'])
            print(len(forward_gene_primer), len(reverse_gene_primer))
        else:
            print(gene)
    primer_df = pd.DataFrame(rows)
    primer_df.columns = ['Name', 'Sequence', 'Scale', 'Purification']
    return primer_df

