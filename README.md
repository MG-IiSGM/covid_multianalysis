# covid_multianalysis
Quality control, alignment, variant calling, annotation, strain imputation and data analysis for SARS-CoV-2 sequences.

- Based on previously developed variant calling pipeline https://github.com/pedroscampoy/provaryota
- Takes in SARS-CoV2 fastq.[gz] files in a folder, a reference genome in gff3 and annotation files and outputs bam files and bamstats, SNV/INDEL info, sample distances, dendrograms and pangolin's COVID lineage.
