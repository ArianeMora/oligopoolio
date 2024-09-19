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


