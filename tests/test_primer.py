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

from oligopoolio.oligopools import *
from oligopoolio.primers import *
from oligopoolio.oligos import *

import unittest
from sciutil import *

u = SciUtil()

test_data_dir = os.path.join(os.path.dirname(__file__), 'data/')


class TestPrimers(unittest.TestCase):

    def test_primer_making(self):
        # Placeholder gene ends (replace with your actual gene sequences)
        gene = "ATGAGCGATCTGCATAACGAGTCCATTTTTATTACCGGCGGCGGATCGGGATTAGGGCTGGCGCTGGTCGAGCGATTTAT\
        CGAAGAAGGCGCGCAGGTTGCCACGCTGGAACTGTCGGCGGCAAAAGTCGCCAGTCTGCGTCAGCGATTTGGCGAACATA\
        TTCTGGCGGTGGAAGGTAACGTGACCTGTTATGCCGATTATCAACGCGCGGTCGATCAGATCCTGACTCGTTCCGGCAAG\
        CTGGATTGTTTTATCGGCAATGCAGGCATCTGGGATCACAATGCCTCACTGGTTAATACTCCCGCAGAGACGCTCGAAAC\
        CGGCTTCCACGAGCTGTTTAACGTCAACGTACTCGGTTACCTGCTGGGCGCAAAAGCCTGCGCTCCGGCGTTAATCGCCA\
        GTGAAGGCAGCATGATTTTCACACTGTCAAATGCCGCCTGGTATCCTGGCGGCGGTGGCCCGCTGTACACCGCCAGTAAA\
        CATGCCGCAACCGGACTTATTCGCCAACTGGCTTATGAACTGGCACCGAAAGTGCGGGTGAATGGCGTCGGCCCGTGTGG\
        TATGGCCAGCGACCTGCGCGGCCCACAGGCGCTCGGGCAAAGTGAAACCTCGATAATGCAGTCTCTGACGCCGGAGAAAA\
        TTGCCGCCATTTTACCGCTGCAATTTTTCCCGCAACCGGCGGATTTTACGGGGCCGTATGTGATGTTGGCATCGCGGCGC\
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

        print(f"Forward Gene Primer: 5'-{forward_gene_primer}-3' (Tm: {forward_tm:.2f} °C)")
        print(f"Reverse Gene Primer: 5'-{reverse_gene_primer}-3' (Tm: {reverse_tm:.2f} °C)")

    def test_get_flanking(self):
        gff_file = f'{test_data_dir}genome_NEB_B/genomic.gff'
        reference_fasta = f'{test_data_dir}genome_NEB_B/GCF_001559615.2_ASM155961v2_genomic.fna'
        gene_name = 'udp'
        seqid, start, end, strand, upstream_flank, downstream_flank, gene_seq = get_flanking_primers(gene_name,
                                                                                                     gff_file,
                                                                                                     reference_fasta)

        assert gene_seq[:3] == 'ATG'  # i.e. methionine
        assert gene_seq[-3:] == 'TAA'
        assert len(upstream_flank) == 50
        assert len(downstream_flank) == 50
        assert 'AAAAC' == downstream_flank[:5]
        assert 'AAGT' == upstream_flank[:4]

    def test_get_flanking_reverse(self):
        gff_file = f'{test_data_dir}genome_NEB_B/genomic.gff'
        reference_fasta = f'{test_data_dir}genome_NEB_B/GCF_001559615.2_ASM155961v2_genomic.fna'
        gene_name = 'ysgA'
        seqid, start, end, strand, upstream_flank, downstream_flank, gene_seq = get_flanking_primers(gene_name,
                                                                                                     gff_file,
                                                                                                     reference_fasta)

        actual_gene_seq = 'ATGGCAACAACACAACAATCTGGATTTGCACCTGC\
            TGCATCGCCTCTCGCTTCGACCATCGTTCAGACTCCGGACGACGC\
            GATTGTGGCGGGCTTCACCTCTATCCCTTCACAAGGGGATAACATGCCTGCTTACCATGCCAGACCAAAGCAAAGCGATG\
            GCCCACTGCCAGTGGTCATTGTAGTGCAGGAAATTTTTGGCGTGCATGAACATATCCGCGATATTTGTCGCCGTCTGGCG\
            CTGGAGGGGTATCTGGCTATCGCACCTGAACTTTACTTCCGCGAAGGCGATCCGAATGATTTTGCCGATATCCCTACGCT\
            GCTTAGCGGTCTGGTAGCAAAAGTGCCTGACTCGCAGGTGCTGGCCGATCTCGATCATGTCGCCAGTTGGGCGTCACGCA\
            ACGGCGGCGATGTTCATCGTTTAATGATCACCGGATTCTGCTGGGGTGGACGTATCACCTGGCTGTATGCCGCGCATAAT\
            CCACAGCTAAAAGCCGCAGTGGCGTGGTACGGCAAGCTGACAGGCGATAAGTCGCTGAATTCACCGAAACAACCTGTTGA\
            TATCGCAACCGATCTTAACGCGCCGGTTCTCGGCTTATATGGCGGCCAGGATAACAGCATTCCGCAAGAGAGCGTGGAAA\
            CGATGCGCCAGGCGCTGCGGGCAGCAAACGCGAAAGCAGAGATTATCGTCTACCCGGATGCCGGGCATGCATTCAACGCC\
            GATTATCGCCCGAGCTATCATGCCGAGTCTGCGAAAGACGGCTGGCAGCGAATGTTGGAATGGTTTACACAGTATGGTGTTAAAAAGTAA'

        print(gene_seq[:20])
        print(actual_gene_seq[:20])
        print(gene_seq[-20:])
        print(actual_gene_seq[-20:])
        assert gene_seq[:20] == actual_gene_seq[:20]
        assert gene_seq[-20:] == actual_gene_seq[-20:]

        assert gene_seq[:3] == 'ATG'  # i.e. methionine
        assert gene_seq[-3:] == 'TAA'
        assert len(upstream_flank) == 50
        print(downstream_flank, upstream_flank)
        assert len(downstream_flank) == 50
        assert 'TACC' == upstream_flank[-4:]
        assert 'CGCC' == downstream_flank[:4]
        # At the moment I'm unsure if the flanks need to be reversed?

    def test_from_fasta(self):
        fasta_file = f'{test_data_dir}example_fasta.fasta'
        df = make_primers_IDT(fasta_file, remove_stop_codon=True, his_tag='',
                              max_length=60, min_length=15, tm_tolerance=30, desired_tm=62.0,
                              forward_primer='gaaataattttgtttaactttaagaaggagatatacat',
                              reverse_primer='ctttgttagcagccggatc')
        assert len(df) == 8
        assert df['Sequence'].values[0] == 'ctttaagaaggagatatacatATGACCATAGACAAAAATTGG'

    def test_single_oligo(self):
        fasta_file = f'{test_data_dir}example_fasta.fasta'
        df = make_oligo_single(fasta_file, forward_primer='gaaataattttgtttaactttaagaaggagatatacat',
                          forward_primer_len=15, reverse_primer='gatccggctgctaacaaag', reverse_primer_len=15,
                          max_len=320)
        assert df['forward_primer'].values[0] == 'gaaggagatatacat'


    def test_double_oligo(self):
        fasta_file = f'{test_data_dir}example_fasta.fasta'

        df = make_oligo_double(fasta_file, forward_primer='gaaataattttgtttaactttaagaaggagatatacat',
                          forward_primer_len=15, reverse_primer='gatccggctgctaacaaag', reverse_primer_len=15,
                          max_len=640,
                          overlap_len=9)

    def test_double_oligo_run(self):
        fasta_file = f'/Users/arianemora/Documents/code/degradeo/manuscript/data/oligopool_random_forest.fasta' #oligopool_random_forest_NT_SwissProt_ActiveSite_prediction_DEHP_25092024.txt'

        df = make_oligo_double(fasta_file, forward_primer='gaaataattttgtttaactttaagaaggagatatacat',
                               forward_primer_len=15, reverse_primer='gatccggctgctaacaaag', reverse_primer_len=15,
                               max_len=670,
                               overlap_len=24)
        df.to_csv('/Users/arianemora/Documents/code/degradeo/manuscript/data/order_oligopool_random_forest_NT_SwissProt_prediction_DEHP_25092024.csv', index=False)

    def test_double_oligo_run_v2(self):
        # Read in a csv and then convert to fasta
        df = pd.read_csv('~/Documents/oligo_seq_order_metagenomics.csv')
        # save the seqs to a fasta file here
        with open('oligo_seq_order_metagenomics.fasta', 'w+') as fout:
            for seq_name, seq in df[['sequence_name', 'codon_optiimized']].values:
                fout.write(f'>{seq_name}\n{seq}\n')
        fasta_file = f'oligo_seq_order_metagenomics.fasta' #oligopool_random_forest_NT_SwissProt_ActiveSite_prediction_DEHP_25092024.txt'

        df = make_oligo_double(fasta_file, forward_primer='cccctctagaaataattttgtttaactttaagaaggagatatacat',
                               forward_primer_len=46, reverse_primer='CTCGAGCACCACCACCACCACCACTGAgatccggctgctaacaaag',
                               reverse_primer_len=46,
                               max_len=670,
                               overlap_len=24)
        df.to_csv('oligo_seq_order_metagenomics_order.csv', index=False)

    def test_double_oligo_run_v3(self):
        # # Read in a csv and then convert to fasta
        # df = pd.read_csv('~/Documents/oligo_seq_order_metagenomics.csv')
        # # save the seqs to a fasta file here
        # with open('oligo_seq_order_metagenomics.fasta', 'w+') as fout:
        #     for seq_name, seq in df[['sequence_name', 'codon_optiimized']].values:
        #         fout.write(f'>{seq_name}\n{seq}\n')
        fasta_file = f'oligo_seq_order_metagenomics.fasta' #oligopool_random_forest_NT_SwissProt_ActiveSite_prediction_DEHP_25092024.txt'

        df = make_oligo_double(fasta_file, forward_primer='cccctctagaaataattttgtttaactttaagaaggagatatacat',
                               forward_primer_len=46, reverse_primer='CTCGAGCACCACCACCACCACCACTGAgatccggctgctaacaaag',
                               reverse_primer_len=46,
                               max_len=670,
                               overlap_len=24)
        df.to_csv('oligo_seq_order_metagenomics_order_rev_comp_2.csv', index=False)

    def test_double_oligo_TWIST_thermofiles(self):
        # # Read in a csv and then convert to fasta
        df = pd.read_csv('~/Documents/optimized_thermofiles.csv')
        # save the seqs to a fasta file here
        with open('optimized_thermofiles.fasta', 'w+') as fout:
            for seq_name, seq in df[['Name', 'Optimized Sequence']].values:
                print(len(seq))
                fout.write(f'>{seq_name}\n{seq}\n')
        fasta_file = f'optimized_thermofiles.fasta' #oligopool_random_forest_NT_SwissProt_ActiveSite_prediction_DEHP_25092024.txt'

        df = make_oligo_double(fasta_file, forward_primer='cccctctagaaataattttgtttaactttaagaaggagatatacat',
                               forward_primer_len=46, reverse_primer='CTCGAGCACCACCACCACCACCACTGAgatccggctgctaacaaag',
                               reverse_primer_len=46,
                               max_len=900,
                               overlap_len=24)
        df.to_csv('oligo_seq_order_optimized_thermofiles_order_rev_comp_2.csv', index=False)

    def test_double_oligo_TWIST(self):
        fasta_file = f'oligo_seq_order_metagenomics.fasta' #oligopool_random_forest_NT_SwissProt_ActiveSite_prediction_DEHP_25092024.txt'

        df = make_oligo_double(fasta_file, forward_primer='cccctctagaaataattttgtttaactttaagaaggagatatacat',
                               forward_primer_len=38, reverse_primer='CTCGAGCACCACCACCACCACCACTGAgatccggctgctaacaaag',
                               reverse_primer_len=46,
                               max_len=670,
                               overlap_len=18)
        df.to_csv('oligo_seq_order_metagenomics_double.csv', index=False)
    def test_double_oligo_many(self):
        fasta_file = f'oligo_seq_order_metagenomics.fasta' #oligopool_random_forest_NT_SwissProt_ActiveSite_prediction_DEHP_25092024.txt'

        df = make_splitty_oligo(fasta_file, forward_primer='cccctctagaaataattttgtttaactttaagaaggagatatacat',
                               forward_primer_len=38, reverse_primer='CTCGAGCACCACCACCACCACCACTGAgatccggctgctaacaaag',
                               reverse_primer_len=46,
                               max_len=670,
                               overlap_len=18)
        df.to_csv('oligo_seq_order_metagenomics_splitty.csv', index=False)

    def test_multi_oligo_rxnfp_chai(self):
        df = pd.read_csv('../predictions_rxnfp_chai/all-sequences.csv')
        # save the seqs to a fasta file here
        with open('../predictions_rxnfp_chai/optimized_rxnfp.fasta', 'w+') as fout:
            for seq_name, seq in df[['Name', 'Full Optimized Sequence']].values:
                fout.write(f'>{seq_name}\n{seq}\n')
        fasta_file = '../predictions_rxnfp_chai/optimized_rxnfp.fasta' #oligopool_random_forest_NT_SwissProt_ActiveSite_prediction_DEHP_25092024.txt'

        df = make_splitty_oligo(fasta_file, forward_primer='cccctctagaaataattttgtttaactttaagaaggagatatacat',
                               forward_primer_len=38, reverse_primer='CTCGAGCACCACCACCACCACCACTGAgatccggctgctaacaaag',
                               reverse_primer_len=46,
                               max_len=670,
                               overlap_len=18)
        df.to_csv('../predictions_rxnfp_chai/oligo_seq_order_rxnfp_splitty.csv', index=False)


    def test_multi_oligo_rxnfp_chai_thermo(self):
        df = pd.read_csv('../predictions_rxnfp_chai/codon_optimized_creep_predictions_thermofiles/optimized.csv')
        # save the seqs to a fasta file here
        with open('../predictions_rxnfp_chai/codon_optimized_creep_predictions_thermofiles/optimized_chemberta.fasta', 'w+') as fout:
            for seq_name, seq in df[['Name', 'Optimized Sequence']].values:
                fout.write(f'>{seq_name}\n{seq}\n')
        fasta_file = '../predictions_rxnfp_chai/codon_optimized_creep_predictions_thermofiles/optimized_chemberta.fasta' #oligopool_random_forest_NT_SwissProt_ActiveSite_prediction_DEHP_25092024.txt'

        df = make_splitty_oligo(fasta_file, forward_primer='cccctctagaaataattttgtttaactttaagaaggagatatacat', # CTCGAGCACCACCACCACCACCACTGA
                               forward_primer_len=38, reverse_primer='gatccggctgctaacaaag',
                               reverse_primer_len=46,
                               max_len=670,
                               overlap_len=18)
        df.to_csv('../predictions_rxnfp_chai/codon_optimized_creep_predictions_thermofiles_splitty.csv', index=False)


    def test_double_oligo_rxnfp_chai(self):
        df = pd.read_csv('../predictions_rxnfp_chai/all-sequences.csv')
        # save the seqs to a fasta file here
        with open('../predictions_rxnfp_chai/optimized_rxnfp_noTAA.fasta', 'w+') as fout:
            for seq_name, seq in df[['Name', 'Optimized Sequence']].values:
                fout.write(f'>{seq_name}\n{seq}\n')
        fasta_file = '../predictions_rxnfp_chai/optimized_rxnfp_noTAA.fasta' #oligopool_random_forest_NT_SwissProt_ActiveSite_prediction_DEHP_25092024.txt'

        df = make_oligo_double(fasta_file, forward_primer='cccctctagaaataattttgtttaactttaagaaggagatatacat',
                               forward_primer_len=38, reverse_primer='CTCGAGCACCACCACCACCACCACTGAgatccggctgctaacaaag',
                               reverse_primer_len=46,
                               max_len=670,
                               overlap_len=18)
        df.to_csv('../predictions_rxnfp_chai/oligo_seq_order_rxnfp_double.csv', index=False)


    def test_long_oligos(self):
        df = pd.read_csv('../scripts/degradeo/predictions_PBP_12042025/pbp_codon_optimized_no_spaces_5.csv')

        min_gc = 0.25
        max_gc = 0.65
        min_tm = 10
        max_tm = 1000
        min_segment_length = 200
        max_segment_length = 250
        max_length = 500
        primer_lower = 'GATCCGGC'.lower()

        output_directory = 'output_long/'
        oligo_df = get_oligos(df, 'Optimized Sequence', 'Name', output_directory,
                              'gaaataattttgtttaactttaagaaggagatatacat',
                              primer_lower, sequence_end='CTCGAGCACCACCACCACCACCACTGA',
                              min_gc=min_gc, max_gc=max_gc, min_tm=min_tm, max_tm=max_tm,
                              min_segment_length=min_segment_length,
                              max_segment_length=max_segment_length,
                              genbank_file="../scripts/degradeo/predictions_PBP_12042025/base-pet22b-base-anm.gb", insert_position=5193, simple=True,
                              codon_optimize=False)

