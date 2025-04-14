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
Date: September 2024
"""
import typer
import sys
import pandas as pd
from oligopoolio import *
from oligopoolio.oligos import *
import os
from typing_extensions import Annotated
from os.path import dirname, join as joinpath

app = typer.Typer()

def fasta_to_df(fasta):
    rows = []
    records = list(SeqIO.parse(fasta, "fasta"))
    done_records = []
    # Remove all the ids
    for record in records:
        new_id = re.sub('[^0-9a-zA-Z]+', '', record.id)
        if new_id not in done_records:
            rows.append([new_id, record.seq])
        else:
            u.warn_p(['Had a duplicate record! Only keeping the first entry, duplicate ID:', record.id])
    df = pd.DataFrame(rows, columns=['id', 'seq'])
    return df
@app.command()
def run(input_file: Annotated[str, typer.Argument(help="Full path to query fasta or csv (note have simple IDs "\
                                                     "otherwise we'll remove all funky characters.). If you pass as a csv you need to also supply the ")],
        prime5f_overlap: Annotated[str, typer.Argument(help="Overlap of the 5prime region for the backbone.")],
        prime3f_overlap: Annotated[str, typer.Argument(help="Overlap of the 3prime region for the backbone "
                                                            "(same strand as 5prime).")],
        seq_start: Annotated[str, typer.Argument(help="Codons to add to the start of the sequence")] = "",
        seq_end: Annotated[str, typer.Option(help="End of the sequence")] = "",
        oligo_start: Annotated[str, typer.Argument(help="Codons to add to the start of each oligo e.g. for amplification.")] = "",
        oligo_end: Annotated[str, typer.Argument(help="Codons to add to the end of each oligo e.g. for amplification.")] = "",
        min_oligo_len: Annotated[int, typer.Option(help="Minimum length of the oligo")] = 100,
        max_oligo_len: Annotated[int, typer.Option(help="Maximum length of the oligo")] = 200,
        optimal_tm: Annotated[int, typer.Option(help="Optimal temperature of the overlap between oligos (C).")] = 62,
        genbank_file: Annotated[int, typer.Option(help="Genbank file in which to add the oligos for visualisation")] = 62,
        output_folder: Annotated[str, typer.Option(help="Where to store results (full path!)")] = 'Current Directory',
        run_name: Annotated[str, typer.Option(help="Name of the run, default is oligos")] = 'oligos',
        id_col: Annotated[str, typer.Option(help="id column in df if df passed (csv) rather than fasta")] = 'id',
        seq_col: Annotated[str, typer.Option( help="codon optimized sequence column.)")] = 'seq',
        ):
    """
    Find similar proteins based on sequence or structural identity in order to annotate these using
    BLAST and FoldSeek. Also annotate with ProteInfer and CLEAN.
    """
    output_folder = output_folder if output_folder != 'Current Directory' else os.getcwd()
    if '.fa' in input_file or '.fasta' in input_file and '.csv' not in input_file:
        df = fasta_to_df(input_file)
    else:
        df = pd.read_csv(input_file)

    # Run the pipeline
    get_oligos(run_name, df, id_col, seq_col, output_folder,
             database, clean_dir=clean_dir, proteinfer_dir=proteinfer_dir,
             run_method=run_method, keep_dups=keep_dups, args_blast=args_blast, args_foldseek=args_foldseek,
             args_proteinfer=args_proteinfer, args_clean=args_clean, methods=methods, foldseek_db=foldseek_db)


if __name__ == "__main__":
    app()

# Example command
# oligopoolio input_df.csv Uniprot_reviewed_catalytic_activity_06032025.fasta --methods blast --output-folder output/ --run-name omgprot50