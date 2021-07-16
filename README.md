# TNSeq-Course-Summer-2021

TNSeq analysis materials for the Summer 2021 course

## Analysis directories

The follow along course materials for the hands-on analysis section of the two days can be found in the `analysis` folder.

* **Day 1:** analysis/preprocess.md
* **Day 2:** analysis/analysis.md

For the first day, the course will go through the pipeline for converting the FASTQ files from the sequencer to count files used for statistical analysis. Additionally, QC methods are incorporated to determine quality of the data. On the second day, the course will outline the statistical analysis (in R) for:

* sliding windows and Wilcoxon Rank Sum statistical tests for localized differences in insertion frequencies
* a gene-based, parametric statistical analysis of insertion frequency between conditions

Separate from the Jupyter lab instance for the data processing and R analysis, the course will showcase Artemis and Transit. This will be done by simple screen share. Artemis is used for visualization of reads against the genome, and Transit will be used to show how to perform an HMM analysis. Both software are included in the conda environment (see below).

## Software environment

Finally, there is the `env` directory, which contains the conda environment used in the course. If you wish to use any of the tools, or to replicate the analysis process of the course on your own client, you can use this conda environment YAML file to recreate the environment - i.e. automatically install with conda all the software used in this course.

```bash
# Create a conda environment by the YAML file
conda env create -f tnseq.yml
# Load the conda environment for access to the software
conda activate tnseq
```

For further details on getting conda installed on your system, please see the [conda install documents](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).