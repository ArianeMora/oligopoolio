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


u = SciUtil()


def optimize_primer3(gene_sequence, upstream_seq='GAATTC', downstream_seq='CTCGAG', label='gene_cloning'):
    """
    https://www.qiagen.com/us/knowledge-and-support/knowledge-hub/bench-guide/pcr/introduction/pcr-primer-design
    Args:
        gene_sequence:
        upstream_seq:
        downstream_seq:
        label:

    Returns:

    """
    # Add restriction site sequences to your primers
    forward_primer_seq = upstream_seq + gene_sequence[:18]  # EcoRI site by default
    reverse_primer_seq = gene_sequence[-18:] + downstream_seq  # XhoI site by default, considering reverse complement

    primer_config = {
        'SEQUENCE_ID': label,
        'SEQUENCE_TEMPLATE': forward_primer_seq + gene_sequence + reverse_primer_seq,
        'SEQUENCE_INCLUDED_REGION': [0, len(forward_primer_seq + gene_sequence + reverse_primer_seq)],
        'PRIMER_TASK': 'generic',
        'PRIMER_PICK_LEFT_PRIMER': 1,
        'PRIMER_PICK_INTERNAL_OLIGO': 0,
        'PRIMER_PICK_RIGHT_PRIMER': 1,
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 58.0,
        'PRIMER_MAX_TM': 62.0,
        'PRIMER_MIN_GC': 40.0,
        'PRIMER_MAX_GC': 60.0,
        'PRIMER_MAX_NS_ACCEPTED': 0,
        'PRIMER_PRODUCT_SIZE_RANGE': [[len(gene_sequence) + 40, len(gene_sequence) + 100]]
    }

    # Design primers
    results = primer3.bindings.designPrimers(primer_config)
    return results


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