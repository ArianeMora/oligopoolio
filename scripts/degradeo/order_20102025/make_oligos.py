import pandas as pd
import random
import numpy 
import sys
from dnachisel import translate
sys.path.append('../../../oligopoolio')
from oligos import get_oligos
SEED = 128
random.seed(SEED)
numpy.random.seed(SEED)  

min_gc = 0.25
max_gc = 0.65
min_tm = 10
max_tm = 1000
min_segment_length = 20
max_segment_length = 100
max_length = 500

df = pd.read_csv('test_order_twist_codon_optimized.csv')
df['Sequence'] = df['Optimized Sequence'].apply(lambda x: translate(x))
df.to_csv('test_order_twist_codon_optimized_translated.csv', index=False)
output_directory = 'output/DNAWeaver/'
oligo_df = get_oligos(df, 'Sequence', 'Name', output_directory, 'gaaataattttgtttaactttaagaaggagatatacat', 'gatccggctgctaacaaag', sequence_end='TAA',
                     min_gc=min_gc, max_gc=max_gc, min_tm=min_tm, max_tm=max_tm, min_segment_length=min_segment_length, max_segment_length=max_segment_length, 
                     genbank_file="base-pet22b-base-anm.gb", insert_position=5193)
oligo_df.to_csv(f'{output_directory}/oligos_DNAWeaver.csv', index=False)
output_directory = 'output/Simple/'

oligo_df = get_oligos(df, 'Sequence', 'Name', output_directory, 'gaaataattttgtttaactttaagaaggagatatacat', 'gatccggctgctaacaaag', sequence_end='TAA',
                     min_gc=min_gc, max_gc=max_gc, min_tm=min_tm, max_tm=max_tm, min_segment_length=min_segment_length, max_segment_length=max_segment_length, 
                     genbank_file="base-pet22b-base-anm.gb", insert_position=5193, simple=True)
oligo_df.to_csv(f'{output_directory}/oligos_simple.csv', index=False)
