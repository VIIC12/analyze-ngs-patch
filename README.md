# README

**Usage:** ./run_ngs.sh [-n Path to the reference sequence] File with samples
**Example:** ./run_ngs.sh -n ref_seq.fasta samples.txt

### <u>SampleFile</u>
Your samples file should look like this:
file_name,sample_name,fw_barcode,rv_barcode,length between barcodes,postion of patch

sample1,igl1_p4,ATCACACGTTCCAGACTACGCTGCTAG,GAACTGACAACTATATGCGAGCTGAGA,301,25

If you have several files (duplicates), you can simply add them to the front:
file_name1,file_name2,...,sample_name,fw_barcode,rv_barcode,length between barcodes,postion of patch

sample1,sample2,igl1_p4,ATCACACGTTCCAGACTACGCTGCTAG,GAACTGACAACTATATGCGAGCTGAGA,30125

### <u>Plotting</u>
The script plots the patch in the zoomed heatmap from 2 positions before and up to 5 positions after the start position. I will perhaps make this customizable again.

---
version = 1.1.0