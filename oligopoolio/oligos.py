import dnaweaver as dw
import time
import dnachisel as dnachisel 
from sciutil import SciUtil
from Bio.Seq import Seq
from difflib import SequenceMatcher
import primer3
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import pandas as pd
from primer3 import calc_hairpin, calc_homodimer
import math
from dnachisel import *

u = SciUtil()

import pyswarms as ps
import Levenshtein
import numpy as np
from functools import partial
import random


def create_oligos(seq_df, seq_col, id_col, genbank_file, bsa, min_overlap, max_overlap, optimal_seq_len, min_seq_len, max_seq_len, tm_tolerance, 
                  reverse=False, insert_position=5193):
    rows = []
    for seq, seq_id in seq_df[[seq_col, id_col]].values:
        oligos, best_cost = optimize_fragments_for_gene(seq, min_overlap, max_overlap, optimal_seq_len, min_seq_len, max_seq_len, tm_tolerance)
        translation_label = f"Insert_{seq_id}"
        record = insert_sequence_with_translation(genbank_file, None, insert_position, seq, translation_label, reverse)
        prev_oligo = None
        if len(oligos) > 0:
            for i, oligo in enumerate(oligos):
                seq = oligo
                primer_overlap, primer_tm, primer_len, homodimer_tm, hairpin_tm = None, None, None, None, None
                if prev_oligo:  
                    match = SequenceMatcher(None, prev_oligo, seq).find_longest_match()
                    primer_overlap = prev_oligo[match.a:match.a + match.size]
                    results = check_secondary_structure(primer_overlap)
                    homodimer_tm = results['homodimer']['homodimer_dg']
                    hairpin_tm = results['hairpin']['hairpin_dg']
                    primer_tm = primer3.bindings.calcTm(primer_overlap)
                    primer_len = len(primer_overlap)
                prev_oligo = seq
                orig_seq = seq
                strand = 1
                if bsa:
                    if i % 2 == 1 and i != 0:
                        seq = str(Seq(seq).reverse_complement())
                        strand = -1
                else:
                    if i % 2 == 0:
                        seq = str(Seq(seq).reverse_complement())
                        strand = -1
                oligo_tm = primer3.bindings.calcTm(seq)
                insert_features_from_oligos(record, f"{seq_id}_oligo_{i}", orig_seq, strand, oligo_tm, None)
                rows.append([seq_id, f"{seq_id}_oligo_{i}", seq, len(seq), oligo_tm, strand, primer_overlap, primer_tm, primer_len, homodimer_tm, hairpin_tm, seq, best_cost])
            if genbank_file:
                output_file = genbank_file.replace('.', f'_{seq_id}.')
                record.name = seq_id
                SeqIO.write(record, f'{output_file}', "genbank")
    oligo_df = pd.DataFrame(rows, columns=["id", "oligo_id", "oligo_sequence", "oligo_length", "oligo_tm", "strand", "primer_overlap_with_previous", "overlap_tm_5prime", "overlap_length", "overlap_homodimer_dg", "overlap_hairpin_dg", "original_sequence", "best_cost"])
    return oligo_df

def optimize_fragments_for_gene(seq, min_overlap, max_overlap, optimal_seq_len, min_seq_len, max_seq_len, tm_tolerance):
    seq_len = len(seq)
    num_fragments = int(math.ceil(seq_len / optimal_seq_len))
    if num_fragments % 2 == 1:
        num_fragments += 1
    num_cuts = num_fragments - 1
    num_variables = num_cuts * 2

    cut_lower_bounds = np.arange(1, num_cuts + 1) * 2
    cut_upper_bounds = np.full(num_cuts, seq_len - 2)
    olap_lower_bounds = np.full(num_cuts, min_overlap)
    olap_upper_bounds = np.full(num_cuts, max_overlap)

    lower_bounds = np.concatenate([cut_lower_bounds, olap_lower_bounds])
    upper_bounds = np.concatenate([cut_upper_bounds, olap_upper_bounds])

    bounds = (lower_bounds, upper_bounds)

    optimizer = ps.single.GlobalBestPSO(
        n_particles=150,
        dimensions=num_variables,
        options={'c1': 0.5, 'c2': 0.5, 'w': 0.7},
        bounds=bounds
    )

    best_cost, best_pos = optimizer.optimize(partial(objective_function, seq_len=seq_len, num_fragments=num_fragments, seq=seq,
                                                     max_overlap=max_overlap, min_overlap=min_overlap,
                                                     max_seq_len=max_seq_len, min_seq_len=min_seq_len, tm_tolerance=tm_tolerance), iters=100)

    cuts = np.sort(np.round(best_pos[:num_cuts]).astype(int))
    overlaps = np.round(best_pos[num_cuts:]).astype(int)
    fragments = get_fragments(cuts, seq_len, num_fragments, overlaps)
    fragment_strings = [seq[start:end] for start, end in fragments]

    return fragment_strings, best_cost

def get_fragments(cuts, seq_len, num_fragments, overlap_lens):
    cuts = np.sort(np.round(cuts).astype(int))
    starts = [0] + list(cuts)
    fragments = []

    for i in range(num_fragments):
        start = starts[i]
        if i < num_fragments - 1:
            next_start = starts[i + 1]
            end = next_start + overlap_lens[i]
        else:
            end = seq_len
        end = min(end, seq_len)
        fragments.append((int(start), int(end)))
    return fragments

def objective_function(x, seq_len, num_fragments, seq, tm_optimal=62, min_dist=8, max_overlap=27, min_overlap=18, 
                       min_seq_len=80, max_seq_len=130, tm_tolerance=5):
    penalties = []
    for particle in x:
        cuts = particle[:num_fragments-1]
        overlaps = particle[num_fragments-1:]
        fragments = get_fragments(cuts, seq_len, num_fragments, overlaps)
        frag_lengths = [end - start for start, end in fragments]

        olap_penalty = 0
        overlaps_seqs = []
        temp_penalty = 0
        seq_len_penalty = 0
        homodimer_penalty = 0
        repeat_penalty = 0
        for i in range(len(fragments) - 1):
            start1, end1 = fragments[i]
            start2, end2 = fragments[i + 1]
            overlap_seq = seq[start2: end1]
            results = check_secondary_structure(overlap_seq)
            homodimer_tm = results['homodimer']['homodimer_dg']
            hairpin_tm = results['hairpin']['hairpin_dg']
            primer_tm = primer3.bindings.calcTm(overlap_seq)

            tm_penalty = np.abs(primer_tm - tm_optimal)
            if len(overlap_seq) < min_overlap or len(overlap_seq) > max_overlap:
                olap_penalty += 10 * abs(len(overlap_seq) - (min_overlap + max_overlap)/2)
            temp_penalty += 2*tm_penalty if tm_penalty > tm_tolerance else tm_penalty

            for frag_len in [end1 - start1, end2 - start2]:
                if frag_len > max_seq_len:
                    seq_len_penalty += 10 * abs(max_seq_len - frag_len)
                elif frag_len < min_seq_len:
                    seq_len_penalty += 10 * abs(min_seq_len - frag_len)

            homodimer_penalty += 2*abs(homodimer_tm) if homodimer_tm < -3 else abs(homodimer_tm)
            if overlap_seq[:3] in ['AAA', 'TTT', 'CCC', 'GGG'] or overlap_seq[-3:] in ['AAA', 'TTT', 'CCC', 'GGG']:
                repeat_penalty += 20
            overlaps_seqs.append(overlap_seq)

        min_dist_val = 100000
        for i, x in enumerate(overlaps_seqs):
            for j, y in enumerate(overlaps_seqs):
                if i != j:
                    dist = Levenshtein.distance(x, y)
                    if dist < min_dist_val:
                        min_dist_val = dist
        dist_penalty = max_overlap - min_dist_val
        if min_dist_val < min_dist:
            dist_penalty += max_overlap

        mean_len = np.mean(frag_lengths)
        variance = np.sqrt(np.sum([(l - mean_len)**2 for l in frag_lengths]))

        penalties.append(variance + olap_penalty + homodimer_penalty + temp_penalty + seq_len_penalty + dist_penalty + repeat_penalty)
    return np.array(penalties)


def insert_sequence_with_translation(input_file, output_file, insert_position, new_sequence, translation_label, reverse=False):
    """
    Insert a new sequence at a specific position in a GenBank file, add a translation annotation,
    and adjust existing annotations to avoid overlap.

    Args:
        input_file (str): Path to the original GenBank file.
        output_file (str): Path to save the modified GenBank file.
        insert_position (int): Position to insert the new sequence (0-based).
        new_sequence (str): The DNA sequence to insert.
        translation_label (str): Label for the translation feature.
        reverse (bool): Whether the feature should be on the reverse strand.
    """
    # Read the original GenBank file
    record = SeqIO.read(input_file, "genbank")
    
    # Insert the new sequence at the specified position
    if reverse:
        new_sequence = str(Seq(new_sequence).reverse_complement())  # Reverse complement the sequence if needed
    record.seq = record.seq[:insert_position] + Seq(new_sequence) + record.seq[insert_position:]
    
    # Adjust existing annotations to avoid overlap
    inserted_length = len(new_sequence)
    for feature in record.features:
        if feature.location.start >= insert_position:
            feature.location = FeatureLocation(
                feature.location.start + inserted_length,
                feature.location.end + inserted_length,
                strand=feature.location.strand
            )
    
    # Create the feature label
    strand_label = " (reverse)" if reverse else " (forward)"
    full_label = translation_label + strand_label
    
    # Add a feature for the inserted sequence
    start = insert_position
    end = insert_position + len(new_sequence)
    feature = SeqFeature(
        location=FeatureLocation(start, end, strand=-1 if reverse else 1),
        type="CDS",  # CDS type for coding sequences
        qualifiers={
            "label": full_label,
            "translation": Seq(new_sequence).translate(table=11)  # Translation for the sequence
        }
    )
    record.features.append(feature)
    # Save the modified GenBank file
    if output_file:
        SeqIO.write(record, output_file, "genbank")
        print(f"Updated GenBank file saved as {output_file}")
    
    return record

def insert_features_from_oligos(record, seq_id, seq, strand, tm, output_file):
    """
    Insert features into a GenBank file based on oligo_df.

    Args:
        genbank_file (str): Path to the original GenBank file.
        output_file (str): Path to save the updated GenBank file.
        oligo_df (pd.DataFrame): DataFrame with oligo information. 
                                 Expected columns: ['id', 'oligo_id', 'oligo_sequence', 'oligo_length', ...].
    """
    # Iterate through the oligo DataFrame to add features
    start = record.seq.find(seq.upper())
    if start == -1:
        print(f"Warning: Oligo sequence {seq_id} not found in the GenBank sequence.")
    
    end = start + len(seq)
    feature = SeqFeature(
        location=FeatureLocation(start, end, strand=strand),
        type="misc_feature",
        qualifiers={
            "label": f"{seq_id} {'(reverse)' if strand == -1 else '(forward)'}",
            "note": f"Length: {len(seq)}, TM: {tm}"
        }
    )
    record.features.append(feature)
    
    # Save the updated GenBank file
    if output_file:
        SeqIO.write(record, output_file, "genbank")
        print(f"Updated GenBank file saved as {output_file}")
    return record
    
def check_secondary_structure(sequence, temp=55):
    """
    Check secondary structures like hairpins and homodimers in a given primer sequence.
    
    Args:
        sequence (str): The DNA sequence of the primer to analyze.
        
    Returns:
        dict: Results for hairpin and homodimer properties.
    """
    try:
        # Check for hairpin structure
        hairpin_result = calc_hairpin(sequence, temp_c=temp)
        hairpin_info = {
            "hairpin_found": hairpin_result.structure_found,
            "hairpin_tm": hairpin_result.tm,
            "hairpin_dg": hairpin_result.dg/1000,
            "hairpin_dh": hairpin_result.dh/1000,
            "hairpin_ds": hairpin_result.ds,
        }

        # Check for homodimer structure
        homodimer_result = calc_homodimer(sequence, temp_c=temp)
        homodimer_info = {
            "homodimer_found": homodimer_result.structure_found,
            "homodimer_tm": homodimer_result.tm,
            "homodimer_dg": homodimer_result.dg/1000,
            "homodimer_dh": homodimer_result.dh/1000,
            "homodimer_ds": homodimer_result.ds,
        }
    except Exception as e:
        u.warn_p([f"Warning: issue with secondary structure check. ", sequence, e])
        hairpin_info = {"hairpin_found": False, "hairpin_tm": None, "hairpin_dg": None, "hairpin_dh": None, "hairpin_ds": None}
        homodimer_info = {"homodimer_found": False, "homodimer_tm": None, "homodimer_dg": None, "homodimer_dh": None, "homodimer_ds": None}
    # Combine results
    return {"hairpin": hairpin_info, "homodimer": homodimer_info}

def build_oligos_DNAWeaver(seq_id: str, sequence: str, output_directory: str, min_gc=0.3, max_gc=0.7, min_tm=55, max_tm=70, min_segment_length=40, max_segment_length=100, max_length=1500):
    """ Use DNAweaver to build oligos """
    # Here we use a comercial supplier but don't actually care. 
    cheap_dna_offer = dw.CommercialDnaOffer(
        name="CheapDNA.",
        sequence_constraints=[
            dw.NoPatternConstraint(enzyme="BsaI"),
            dw.SequenceLengthConstraint(max_length=4000),
            dw.GcContentConstraint(min_gc=min_gc, max_gc=max_gc)
        ],
        pricing=dw.PerBasepairPricing(0.10),
    )

    oligo_dna_offer = dw.CommercialDnaOffer(
        name="OliGoo",
        sequence_constraints=[
            dw.GcContentConstraint(min_gc=min_gc, max_gc=max_gc),
            dw.SequenceLengthConstraint(max_length=4000),
        ],
        pricing=dw.PerBasepairPricing(0.07),
        memoize=True
    )

    oligo_assembly_station = dw.DnaAssemblyStation(
        name="Oligo Assembly Station",
        assembly_method=dw.OligoAssemblyMethod(
            overhang_selector=dw.TmSegmentSelector(
                min_size=15, max_size=25, min_tm=min_tm, max_tm=max_tm
            ),
            min_segment_length=min_segment_length,
            max_segment_length=max_segment_length,
            sequence_constraints=[dw.SequenceLengthConstraint(max_length=4000)],
            duration=8,
            cost=30,
        ),
        supplier=oligo_dna_offer,
        coarse_grain=20,
        a_star_factor="auto",
        memoize=True,
    )

    assembly_station = dw.DnaAssemblyStation(
        name="Gibson Assembly Station",
        assembly_method=dw.GibsonAssemblyMethod(
            overhang_selector=dw.TmSegmentSelector(min_tm=min_tm, max_tm=max_tm),
            min_segment_length=min_segment_length,
            max_segment_length=max_segment_length + 20, # add a bit of a buffer
        ),
        supplier=[cheap_dna_offer, oligo_assembly_station],
        logger="bar",
        coarse_grain=100,
        fine_grain=10,
        a_star_factor="auto",
    )
    
    print("Looking for the best assembly plan...")
    t0 = time.time()
    quote = assembly_station.get_quote(sequence, with_assembly_plan=True)
    assembly_plan_report = quote.to_assembly_plan_report()
    assembly_plan_report.write_full_report(f"{output_directory}/oligo_assembly_plan_{seq_id}.zip")
    original_sequence = assembly_plan_report.plan.sequence
    # Then get the sequence 
    rows = []
    for oligo in assembly_plan_report.plan.assembly_plan:
        # If this was chosen then choose it
        if oligo.accepted:
            rows.append([oligo.id, oligo.sequence, original_sequence])
    return rows


def codon_optimize(protein_sequence: str, min_gc=0.3, max_gc=0.7):
    """ Codon optimize the protein sequence using DNA chisel: https://github.com/Edinburgh-Genome-Foundry/DnaChisel"""
    seq = dnachisel.reverse_translate(protein_sequence)
    problem = dnachisel.DnaOptimizationProblem(
        sequence=seq,
        constraints=[
            AvoidPattern("BsaI_site"), # type: ignore
            EnforceGCContent(mini=min_gc, maxi=max_gc, window=50), # type: ignore
            EnforceTranslation(location=(0, len(seq))), # type: ignore
            AvoidStopCodons(location=(0, len(seq)-3)) # type: ignore Let there be stop codons in the last bit
        ],
        objectives=[CodonOptimize(species='e_coli', location=(0, len(seq)))] # type: ignore
    )
    # SOLVE THE CONSTRAINTS, OPTIMIZE WITH RESPECT TO THE OBJECTIVE
    problem.resolve_constraints()
    problem.optimize()

    # PRINT SUMMARIES TO CHECK THAT CONSTRAINTS PASS
    print(problem.constraints_text_summary())
    print(problem.objectives_text_summary())

    # GET THE FINAL SEQUENCE (AS STRING OR ANNOTATED BIOPYTHON RECORDS)
    final_sequence = problem.sequence  # string
    final_record = problem.to_record(with_sequence_edits=True)
    print(protein_sequence)
    print(final_sequence)
    return final_sequence