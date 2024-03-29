# README

**Usage:** ./run_ngs.sh [-n <Path to the reference sequence>] [-s <Path to the samples file] [-o overwrite] <br/>
**Example:** ./run_ngs.sh -n ref_seq.fasta -s samples.txt

### <u>SampleFile</u>
Each line is a unique sample, therefore each sample_name has to be unique. <br/>
Your samples file should look like this: <br/>
file_name,sample_name,fw_barcode,rv_barcode,length between barcodes,postion of patch <br/>
file_name,sample_name2,fw_barcode,rv_barcode,length between barcodes,postion of patch <br/>
... <br/>

If you have several files (duplicates), you can simply add them to the front: <br/>
file_name1,file_name2,...,sample_name,fw_barcode,rv_barcode,length between barcodes,postion of patch <br/>

Example:<br/>
`sample1,igl1_p1,ATCACACGTTCCAGACTACGCTGCTAG,GAACTGACAACTATATGCGAGCTGAGA,225,25
sample1,sample2,igl1_p2,ATCACACGTTCCAGACTACGCTGCTAG,GAACTGACAACTATATGCGAGCTGAGA,301,25`

#### <u>Barcodes and length between barcodes</u>
![Explanation length between barcodes](https://docs.google.com/drawings/d/1PCBN5wQSFBwD71xal4ZbwcH0CPPnzmCGk6g9Fs3Qf4k/export/png)
The easiest way is to use the sequence barcode + overlap (red+grey) for fw_barcode and rv_barcode. The _length between barcodes_ then corresponds exactly to your construct length.


### <u>Plotting</u>
The script plots the patch in the zoomed heatmap from 2 positions before and up to 5 positions after the start position. I will probably make this customizable soon.

---
version = 1.2.0