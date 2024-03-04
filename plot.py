#!/usr/bin/env python
#v1.2.0
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import warnings
warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", category=UserWarning, module="matplotlib.font_manager")

# For automatization
def create_list_runs_from_file(filename: str, linenumber: int):
    """
    Converts line of samples.txt to list_runs
    Example:
        From:   sample1,sample2,igl1_p4,ATCACACGTTCCAGACTACGCTGCTAG,GAACTGACAACTATATGCGAGCTGAGA,301,25
        to:     [['igl1_p4_S1', 'igl1_p4_S2', '23', '30']]
    """
    with open(filename, 'r') as sampe_file:
        samples = sampe_file.read().splitlines()

    list_runs = [
    [fields[-5] + "_S" + str(i)
     for i in range(1, len(fields) - 4)] +
    [str(int(fields[-1]) - 2), str(int(fields[-1]) + 5)]
    for fields in [samples[linenumber].split(',')]
    ]

    return list_runs

#? argv
    # 0: plot.py
    # 1: frameshift
    # 2: reference_sequence
    # 3: sample_file name
    # 4: line number

reference_sequence = str(SeqIO.read(sys.argv[2], 'fasta').seq.translate())
if int(sys.argv[1]) == 1:
    frameshift = None
else:
    frameshift = int(sys.argv[1])-1

list_runs = create_list_runs_from_file(sys.argv[3], int(sys.argv[4]))

for run_num, runs in enumerate([run_only[:-2] for run_only in list_runs]):
    #? Generate mutation dicts
    # for plotting
    mutation_table_full = {f"{ref_aa}{idx}": {aa: 0 for aa in "ACDEFGHIKLMNPQRSTVWXY*"} for idx, ref_aa in enumerate(reference_sequence, start=1)}
    # for rel-calculations
    native_table_full = {f"{ref_aa}{idx}": {aa: 0 for aa in "ACDEFGHIKLMNPQRSTVWXY*"} for idx, ref_aa in enumerate(reference_sequence, start=1)}

    for run in range(len(runs)):
        current_file = [str(record.seq) for record in SeqIO.parse(f"demultiplex/{runs[run]}.fasta", 'fasta')]
        #? Generate mutation dicts
        # for plotting
        mutation_table_run = {f"{ref_aa}{idx}": {aa: 0 for aa in "ACDEFGHIKLMNPQRSTVWXY*"} for idx, ref_aa in enumerate(reference_sequence, start=1)}
        # for rel-calculations
        native_table_run = {f"{ref_aa}{idx}": {aa: 0 for aa in "ACDEFGHIKLMNPQRSTVWXY*"} for idx, ref_aa in enumerate(reference_sequence, start=1)}

        #? Iterate over sequences
        for seq in current_file:
            nucleotide_seq_frame = seq[frameshift:]
            seq_aa = str(Seq(nucleotide_seq_frame).translate())

            for res_num, (ref_aa, seq_aa) in enumerate(zip(reference_sequence, seq_aa), start=1):
                #? count only when native-mutation
                if ref_aa != seq_aa:
                    mutation_table_run[f"{ref_aa}{res_num}"][seq_aa] += 1
                #? count all positions
                native_table_run[f"{ref_aa}{res_num}"][seq_aa] += 1

        #? Merging to table_full
        for key in mutation_table_full:
            mutation_table_full[key] = {k: mutation_table_full[key].get(k, 0) + mutation_table_run[key].get(k, 0) for k in set(mutation_table_full[key]) | set(mutation_table_run[key])}
        for key in native_table_full:
            native_table_full[key] = {k: native_table_full[key].get(k, 0) + native_table_run[key].get(k, 0) for k in set(native_table_full[key]) | set(native_table_run[key])}
  
    #? Normalization
    for key in mutation_table_full:
        for inner_key in mutation_table_full[key]:
            mutation_table_full[key][inner_key] /= len(runs) if len(runs) > 0 else 1
    for key in native_table_full:
        for inner_key in native_table_full[key]:
            native_table_full[key][inner_key] /= len(runs) if len(runs) > 0 else 1

    #? Rel-calculations
    frequenc_rel = {}
    # Total number of all sequences (sum of all counts per position)
    total_seq_count = sum(native_table_full[list(native_table_full)[0]].values())

    for (native_out, native_in), (mut_out, mut_in) in zip(native_table_full.items(), mutation_table_full.items()):
        # Which mutation occurs most frequently; non-observance of the native AA
        max_mut_key = max(mut_in, key=lambda k: mut_in[k])
        # How often does the most common mutation occur
        max_mut_value = mut_in[max_mut_key]
        # How often does the native AA occur
        max_native_value = native_in[native_out[0]]
        # How often does the most frequent mutation occur relative to native AA
        rel_native = round((max_mut_value) / native_in[native_out[0]], 2)
        # How often does the most frequent mutation occur relative to all counts
        rel_all = round((max_mut_value) / total_seq_count, 2)
        # rel_native and rel_all in new dict
        frequenc_rel[native_out] = rel_native, rel_all

    #? Reindex dict for right AA plotting order
    desired_order = [a for a in "ACDEFGHIKLMNPQRSTVWXY*"]
    mutation_df = pd.DataFrame(mutation_table_full).reindex(desired_order)
    freq_rel = pd.DataFrame(frequenc_rel)

    #? Select subset for nnk_positions, set values below a treshold to "-"
    possible_positions = list(mutation_table_full.keys())
    selected_positions = mutation_df.loc[:, possible_positions[int(list_runs[run_num][-2]):int(list_runs[run_num][-1])]]
    best_aa = pd.DataFrame({'Best AA': mutation_df.idxmax(axis=0)[int(list_runs[run_num][-2]):int(list_runs[run_num][-1])],
                            'Occ. best AA rel to native': [value[0] if value[0] >= 0.05 else "-" for key, value in list(frequenc_rel.items())[int(list_runs[run_num][-2]):int(list_runs[run_num][-1])]],
                            'Occ. best AA rel to all': [value[1] if value[1] >= 0.05 else "-" for key, value in list(frequenc_rel.items())[int(list_runs[run_num][-2]):int(list_runs[run_num][-1])]]})

    #? Plotting
    header = f"NGS results {runs}"
    font = 'Hiragino Sans GB'
    
    #? Subplots and ratios
    fig, axs = plt.subplots(3, figsize=(15, 10), gridspec_kw={'height_ratios': [2, 1, 1]})
    ax1 = axs[0]
    ax2 = axs[1]
    ax3 = axs[2]

    #? Large plot for the selected positions
    sns.heatmap(data=selected_positions, cmap='viridis', linewidths=0.5, linecolor='slategrey', cbar_kws={'pad': 0.01}, ax=ax1)
    ax1.set_title(header, fontsize=25, fontname=font, fontweight='bold', loc='left')
    ax1.set_ylabel('AA', fontsize=15, fontname=font, style='italic')
    ax1.set_xlabel('Native', fontsize=15, fontname=font, style='italic')

    #? Table
    table1 = ax2.table(cellText=best_aa.T.values, rowLabels=best_aa.T.index ,colLabels=best_aa.T.columns, loc='center', cellLoc='center', colColours=['green']*len(best_aa), bbox=[0.0, 0.3, 0.84, 0.4])
    ax2.axis('off')

    #? Smaller plot for the entire mutation_df
    sns.heatmap(data=mutation_df, cmap='viridis', linewidths=0.5, linecolor='slategrey', cbar_kws={'pad': 0.01}, ax=ax3)
    ax3.set_title("All Data", fontsize=20, fontname=font, fontweight='bold', loc='left')
    ax3.set_ylabel('AA', fontsize=12, fontname=font, style='italic')
    ax3.set_xlabel('Native', fontsize=12, fontname=font, style='italic')
    
    plt.savefig(f'results/{runs[0][:-3]}.png', dpi=400, bbox_inches="tight", pad_inches=1)
    #print(f'{runs[0][:-3]} plotted')

__author__ = "Tom U. Schlegel"
__contact__ = "tom.schlegel@uni-leipzig.de"
__license__ = "GNU GPLv3"
__version__ = "1.2"
