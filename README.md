## README

I use the FULL fw_barcode and rv_barcode sequence (with the N/C overlap) → ReadingFrame should then always be 1. The length between the barcodes is then exactly the original construct length.
samples.txt should look like this:
sample1,igl1_p4_S1,ATCACACGTTCCAGACTACGCTGCTAG,GAACTGACAACTATATGCGAGCTGAGA,301,25
sample2,igl1_p4_S2,ATCACACGTTCCAGACTACGCTGCTAG,GAACTGACAACTATATGCGAGCTGAGA,301,25

The consecutive lines are understood as duplicates.
The script plots the patch in the zoomed heatmap from 2 positions before and up to 5 positions after the start position. I might adjust that again.
