# Preprocess your FASTQ data

This day of the course will go through the entire process of converting the FASTQ reads to input files for statistical analysis with QC along the way.

## Set up your work environment

We will organize the work directory as good practice. A few pointers. \* is a wildcard. It will grab anything. If you add it to beginning `.fastq` like so: `\*.fastq`, then any file that ends with `.fastq` will be grabbed. Alternatively, if we add \* to the end, like with `B_sub\*`, then we will grab any files that start with `B_sub`.

```bash
cd ~ # Brings you back to the "home" folder
# List all the files in the directory at the moment
ls | cat

# This should be the output
0h_1.fastq.gz
24h_1.fastq.gz
5h_1.fastq.gz
B_subtilis_subtilis_s168.1.ebwt
B_subtilis_subtilis_s168.2.ebwt
B_subtilis_subtilis_s168.3.ebwt
B_subtilis_subtilis_s168.4.ebwt
B_subtilis_subtilis_s168.TA.bed
B_subtilis_subtilis_s168.fasta
B_subtilis_subtilis_s168.genome
B_subtilis_subtilis_s168.gff3
B_subtilis_subtilis_s168.rev.1.ebwt
B_subtilis_subtilis_s168.rev.2.ebwt
bootstrapped_counts.tsv
```

Then we create the folders.

```bash
mkdir ref # create a folder for your reference files
mkdir fastq # create a folder for your FASTQ files
mkdir aln # create a folder for your SAM/BAM files
mkdir counts # create a folder for your output counts
mkdir logs # create a folder to hold all the output log files
```

Finally, move the files to their respective folders to organize everything.

```bash
# Move all the B. subtilis reference files to /ref/
mv B_subtilis* ref/
# Move the FASTQS to the fastq/ folder.
mv *.fastq.gz  fastq/
# Move the bootstrap counts file to the counts/ folder
mv bootstrapped_counts.tsv counts/

# if you perform the following command which lists all your files
# (that are not special hidden files in the OS; denoted by starting with '.')
find . -not -path '*/\.*' -type f
# This should be your output
./fastq/24h_1.fastq.gz
./fastq/5h_1.fastq.gz
./fastq/0h_1.fastq.gz
./counts/bootstrapped_counts.tsv
./ref/B_subtilis_subtilis_s168.genome
./ref/B_subtilis_subtilis_s168.rev.2.ebwt
./ref/B_subtilis_subtilis_s168.fasta
./ref/B_subtilis_subtilis_s168.1.ebwt
./ref/B_subtilis_subtilis_s168.4.ebwt
./ref/B_subtilis_subtilis_s168.rev.1.ebwt
./ref/B_subtilis_subtilis_s168.2.ebwt
./ref/B_subtilis_subtilis_s168.TA.bed
./ref/B_subtilis_subtilis_s168.3.ebwt
./ref/B_subtilis_subtilis_s168.gff3
```

Modify the .fastq.gz to whatever the ending of your FASTQ files are. Common alternatives are:
* .fastq
* .fq.gz
* .fq

## SRA access

**This is note necessary for the course. It simply shows how I got the data from the publication.** 

SRP: SRP066259

| Sample | SRA accession | Reads    |
| ---    | ---           | ---      |
| 0h     | SRR2918149    | 12235150 |
| 5h     | SRR2918152    | 12201459 |
| 24h    | SRR2918155    | 11217733 |

To download the files, you can use the NCBI toolkit `fastq-dump` (https://ncbi.github.io/sra-tools/fastq-dump.html) to download the files. `fastq-dump` will download the .sra file and extract the FASTQ data from it.

`fastq-dump` parameters:
* --gzip
    * output the FASTQ files as compressed gzipped files
* --split-files
    * if there are multiple reads (as in paired reads), split the paired reads into their own files

```bash
for srr in SRR2918149 SRR2918152 SRR2918155; do # Loop through every SRR name
    # ${srr} means use this variable stand in for SRR names
    fastq-dump --gzip --split-files ${srr};
done;
mv SRR2918149_1.fastq.gz 0h_1.fastq.gz;
mv SRR2918152_1.fastq.gz 5h_1.fastq.gz;
mv SRR2918155_1.fastq.gz 24h_1.fastq.gz;
# * is a wildcard, so .fastq.gz is the only required match
# anything before that will be counted
mv *.fastq.gz ./fastq
```

## Reference genome

Bacillus subtilis subsp. subtilis str. 168 complete genome: [https://www.ncbi.nlm.nih.gov/nuccore/NC_000964.3](https://www.ncbi.nlm.nih.gov/nuccore/NC_000964.3)

## Library quality control

[FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is a popular program for basic QC of a sequencing library. It can be used as a GUI or a command line (CLI) program. If using it as a GUI, simply open the program, and drag-and-drop the FASTQs into the program.

We will use the program, via CLI) to QC all the libraries. FASTQC also comes with a GUI component.

### The programmatic way

```bash
cd ~/fastq # Move to where your FASTQs are
for fq in *.fastq.gz; do # Loop through every FASTQ
    # ${fq} means use this variable that stands for an individual FASTQ from 
    # *.fastq.gz
    # -o output files to this directory
    #    We send the output to the logs folder we created earlier
    fastqc ${fq} -o ../logs;
done; # Done with loop

# Check to see your logs were properly created and placed.
ls ~/logs/
# The results of your logs should look as follows
0h_1_fastqc.html
0h_1_fastqc.zip
24h_1_fastqc.html
24h_1_fastqc.zip
5h_1_fastqc.html
5h_1_fastqc.zip
```

### Example of the manual way

```bash
fastqc 0h_1.fastq.gz -o ~/logs;
```


## Adapter trimming

Here we will use cutadapt to remove the Illumina adapters and to trim the reads.

In cutadapt, we provide the 3' and paired 3' adapter with -a and -A respectively. For Illumina adapters, this is set as: `-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT`. We also define the output files for the first and second reads with -o and -p respectively. For the input and output, cutadapt will automatically detect if they are compressed by the inclusion of .gz at the end of the file. (Note: this only works for gzip compression.) When the `-o` and `-p` are used, the cutadapt report will be sent to stdout / the terminal. We can redirect this output to a file with `>`.

For specific restriction digest fragmentation library construction, you may want to trim the read length. For example, MmeI restriction digestion for library fragmentation should, in theory, produce reads of exactly 20bp. Therefore, it makes sense to trim reads to 20bp. This can improve specificity of the reads. We will be giving a littel buffer around the read size to account for small indels created in read sequencing, so we will allow reads of size 19-21bp. We have a choice to make:
* To trim reads only, thereby keeping all reads, but trimming all reads to 21bp. It will also keep and do nothing to reads shorter than 21bp.
* Throw away reads that do not fall within the intended range.
In theory, the latter would be the most specific option, as reads should be within the indicated range.

* Trim
    * `--length 21`
* Throw away
    * `--minimum-length 19 --maximum-length 21`

We will proceed with the latter method. It will create a set of reads with high specifity for the targets we are desire.

**However for this dataset, we will forego trimming. The reads are too short - 16bp.**

```bash
gzip -dc 0h_1.fastq.gz | head -4
>@SRR2918149.1 DGL9ZZQ1:641:C5A53ACXX:8:1101:1704:2185 length=16
>CGCATCAGATAATGAT
>+SRR2918149.1 DGL9ZZQ1:641:C5A53ACXX:8:1101:1704:2185 length=16
>FBDE@:<EEF?CGF<?
```

### The programmatic way

Parameters for cutadapt
* `--minimum-length`
    * Minimum allowed length for the read (filters reads)
* `--maximum-length`
    * Maximum allowed length for the read (filters reads)
* `--length`
    * Trim reads to this length 
* `-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA`
    * for indicating the Illumina Truseq adapter sequence for the first read / single end read
* `-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT` 
    * for indicating the Illumina adapter sequence for the paired read
* `-p ${second_trimmed}`
    * output name for second read

#### Single end sequencing (no trimming)

For this workshop, use the following code. This is because the reads are 16bp, so trimming to a defined range 19bp <= lib <= 21bp would remove **ALL** reads. For longer read sequencing, and MmeI mariner based TnSeq protocols, you may want to include the size trimming.

```bash
for first in *_1.fastq.gz; do # Loop through every first read FQ
    # Define the output trimmed read file names.
    first_trimmed=${first%_1.fastq.gz}.trimmed_1.fastq.gz;
    # Define the name of the cutadapt report.
    report=${first%_1.fastq.gz}.cutadapt.log;
    # Run cutadapt
    cutadapt \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -o ${first_trimmed} \
        ${first} > ../logs/${report};
done;

# Look at the files in the FASTQ directory now
ls 
# You should have the following:
0h_1.fastq.gz
0h.trimmed_1.fastq.gz
24h_1.fastq.gz
24h.trimmed_1.fastq.gz
5h_1.fastq.gz
5h.trimmed_1.fastq.gz
```

The following code snippets are for alternative experimental setups on TnSeq. If you have any questions about which 

#### Single end sequencing with size filtering

```bash
for first in *_1.fastq.gz; do # Loop through every first read FQ
    # Define the output trimmed read file names.
    first_trimmed=${first%_1.fastq.gz}.trimmed_1.fastq.gz;
    # Define the name of the cutadapt report.
    report=${first%_1.fastq.gz}.cutadapt.log;
    # Run cutadapt
    cutadapt \
        --minimum-length 19 \
        --length 21 \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -o ${first_trimmed} \
        ${first} > ../logs/${report};
done;
```

#### Paired end sequencing with size trimming

```bash
for first in *_1.fastq.gz; do # Loop through every first read FQ
    # Define the paired read file based on the first's name
    second=${first%_1.fastq.gz}_2.fastq.gz;
    # Define the output trimmed read file names.
    first_trimmed=${first%_1.fastq.gz}.trimmed_1.fastq.gz;
    second_trimmed=${second%_1.fastq.gz}.trimmed_1.fastq.gz;
    # Define the name of the cutadapt report.
    report=${first%_1.fastq.gz}.cutadapt.log;
    # Run cutadapt
    cutadapt \
        --minimum-length 19 \
        --maximum-length 21 \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -o ${first_trimmed} \
        -p ${second_trimmed} \
        ${first} ${second}> ../logs/${report};
done;
```

## Alignment

There are a large number of short read aligners out there. Bowtie is one choice that is very good at reads 35bp or less. For longer reads, you may consider other read aligners if you want to consider speed improvements and improvements in handling mismatches and indels, which become more likely with the longer reads. Other suggested aligners than are:
* bowtie2
* bwa

### bowtie alignment

Before you can align to your genome of interest, you will need to index genome's FASTA file.

`-f` indicates that the input reference is in FASTA format.

```bash
cd ../ref/
bowtie-build \
    -f \
    B_subtilis_subtilis_s168.fasta \
    B_subtilis_subtilis_s168
```

Once you have the genome indexed, you can proceed with aligning the FASTQs.

The bowtie parameters we will be sending is as follows:
* `-v 3`
    * report end-to-end hits w/ <=v mismatches; ignore qualities
    * no more than 3 mismatches
* `-a`
    * report all alignments per read (much slower than low -k)
    * report all 
    * This can significantly slow down the runtime
* `--best`
    * hits guaranteed best stratum; ties broken by quality
    * Alignments are separated into strata by their quality of alignment. Thus creating a best strata.
* `--strata`
    * hits in sub-optimal strata aren't reported (requires --best)
    * Only show reads alignments in the best strata.
    * In combination with `-a` and `--best`, we choose to show all read alignments that are "best" alignments.
* `-m 1`
    * suppress all alignments if > <int> exist (def: no limit)
    * Hide any reads that have morethan `1` alignment
    * It may seem contradictory to the above `-a --best --strata`, but what we have done is force bowtie to not just choose a "best" alignment, but show all equivalent alignments. If we have many ties for "best" alignment, we effectively throw that read away as we cannot determine it has a unique placement.
* `-q`
    * query input files are FASTQ .fq/.fastq (default)
    * Tells bowtie to expect a FASTQ file.
        * While today this may seem odd, bowtie was written in a time when Illumina was not the ubiquitious, gold standard, and heavily standardized sequencing technology that it is today.
* `--sam`
    * outputs the alignment in SAM format
    * SAM has become the universal output alignment format, for which many downstream tools expect. (As opposed to the old, original bowtie output.)

Note that bowtie outputs the data to the screen, so you must send it to a file with `>`. Additionally, bowtie sends its statistics to the screen as an "error," so we send that to a file with `2>`.

```bash
cd ../fastq;
for first in *.trimmed_1.fastq.gz; do
    base=${first%_1.fastq.gz};
    bowtie \
        -v 3 \
        -a \
        --best \
        --strata \
        -m 1 \
        -q \
        --sam \
        -x ../ref/B_subtilis_subtilis_s168 \
        ${first} > ${base}.sam 2> ${base}.bowtie_stats.txt;
done;
```

We've finished aligning, but there is a little more pre-processing to do. We need to convert the text-based SAM file to a smaller, faster binary blob BAM (that we cannot read without special programs). This will also be coordinate-sorted.

```bash
# Organize your files
mv *.bowtie_stats.txt ../logs
mv *.sam ../aln
```

### Alternative aligner

`bwa` is another short read aligner. If your read lengths surpass 35bp, I would suggest using one of `bwa` or `bowtie2`.

```bash
# Index the genome
bwa index B_subtilis_subtilis_s168.fasta
bwa mem \
    B_subtilis_subtilis_s168.fasta \
    required_first_read.fastq \
    optional_second_read.fastq \
> sample.sam 
```

## Alignment QC and statistics

`samtools view` will convert the SAM file to the binary BAM file, but will output is to "text" on the screen. We must point it to the file to which it would be read. However, we can instead send it to another program, `samtools`, to have it sort it first. `samtools sort` recognizes that it is being sent the file with `-`. This is output to a file with `-o`. We need to index this binary with `samtools index` for rapid search in the file.

```bash
cd ../aln
for sam in *.sam; do
    bam=${sam%.sam}.bam;
    # Convert SAM to BAM and sort
    samtools view \
        -bht ../ref/B_subtilis_subtilis_s168.fasta \
        ${sam}\
    | samtools sort -o ${bam} -;
    # Index the files
    samtools index ${bam};
    # Get Alignment statistics with samtools
    samtools idxstats ${bam} > ../logs/${bam%.bam}.idx_stats.txt
done;
```

## Counting

#### File formats

The course aims to make use of common standards for output files, specifically BED and WIG files. These are standards that are accepted by a large swathe of bioinformatics software. Therefore, it is good practice to try to fit your sequencing data into these formats, so that you may use the wide range of tools available to further process your data.

The BED format is used to define features in the genome. It is not required to be continuous nor regular in feature size. The BED format is minimally defined as (example):

| contig | start | end |
| ---    | ---   | --- |
| chrom1 | 1     | 100 |
| chrom1 | 150   | 200 |
| chrom2 | 1     | 100 |

A header is not included. You can find the full specification at: [https://m.ensembl.org/info/website/upload/bed.html](https://m.ensembl.org/info/website/upload/bed.html).

The WIG format is designed to display continuous data, such as probability scores or depth. An example would be:

```
variableStep chrom=chr2
300701  12.5
300702  12.5
300704  22.0
300704  12.5
```

A header is usually included to define the step and contig. You can find the full specification at: [https://m.ensembl.org/info/website/upload/wig.html](https://m.ensembl.org/info/website/upload/wig.html). The wig will be used as the import data file for Transit.

### By TA sites

`bedtools` is a powerful and open source tool to manipulate BED files with BAMs, GTFs, etc. This means less use of custom scripts, and instead a consistent program with standard file formats.
It is worth going through the [bedtools documents](https://bedtools.readthedocs.io/en/latest/index.html) to get an understanding of the breadth of tools within the bedtools suite.

Also, here I have included the only custom script in this process `dinucleotide_loci_from_FASTA.py`. I created this script to catalogue the TA sites in the genome as a BED file. It accepts FASTA genome files. It has an option to change the dinucleotide sequence.

```bash
# First create the TA loci BED file
python dinucleotide_loci_from_FASTA.py \
    ../ref/B_subtilis_subtilis_s168.fasta \
> ../ref/B_subtilis_subtilis_s168.TA.bed;

# Then create a .genome file
# This lists the contig sizes for the genome
# This is only necessary for the optional counting step later
samtools view -H 0h.trimmed.bam
>@HD	VN:1.0	SO:coordinate
>@SQ	SN:NC_000964.3	LN:4215606
>@PG	ID:Bowtie	VN:1.3.0	CL:"/home/drdeconti/bin/miniconda3/envs/tnseq/bin/bowtie-align-s --wrapper basic-0 -v 3 -a --best --strata -m 1 -q --sam -x ../ref/B_subtilis_subtilis_s168 0h.trimmed_1.fastq.gz"
>@PG	ID:samtools	PN:samtools	PP:Bowtie	VN:1.13	CL:samtools view -bht ../ref/B_subtilis_subtilis_s168.fasta 0h.trimmed.sam
>@PG	ID:samtools.1	PN:samtools	PP:samtools	VN:1.13	CL:samtools sort -o 0h.trimmed.bam -
>@PG	ID:samtools.2	PN:samtools	PP:samtools.1	VN:1.13	CL:samtools view -H 0h.trimmed.bam

# Write a tab-delimited file outlining the contig sizes
# @SG will list out all contigs and their respective lengths
echo -e "NC_000964.3\t4215606" > ../ref/B_subtilis_subtilis_s168.genome

# bedtools will find all the reads that overlap the TA BED we created
for bam in *.bam; do
    bedtools intersect \
        -c \
        -sorted \
        -a ../ref/B_subtilis_subtilis_s168.TA.bed \
        -b ${bam} > ${bam%.bam}.TA_counts.bed;
done;

# Optionally, we may want a more unique read counting
# We know that the TA insertion locus is found near the 
# beginning of the read.
# So, let us only consider the central part of the read.
# To do this, we will trim the read from the right by 8 with bedtools.
# 1. Convert BAM to BED
# 2. Trim the BAM's BED file by 2
#    -s = trim length by strand orientation
#         so reads on the - strand will reverse -l/-r
#    -l = add this length off the left side of read (add -8)
#    -r = add this length off the right side of read (0)
# 3. Sort BED (just to be safe)
# 4. Send into the intersect
for bam in *.bam; do
    bedtools bamtobed -i ${bam} \
    | bedtools slop \
        -i - \
        -g ../ref/B_subtilis_subtilis_s168.genome \
        -l 0 -r -12 -s \
    | bedtools sort -i - \
    | bedtools intersect \
        -c \
        -sorted \
        -a ../ref/B_subtilis_subtilis_s168.TA.bed \
        -b - > ${bam%.bam}.trimmed_TA_counts.bed;
done;
```


### By genes

```bash
# Use featureCounts to count the gene coverage
# -t use 'gene' feature annotation
# -g use 'ID' in feature column of GFF3 to get gene name
#       This may be different per GFF/GTF, so read the file first
featureCounts \
    -t gene \
    -g ID \
    -a ../ref/B_subtilis_subtilis_s168.gff3 \
    -o gene_counts.txt \
    *.bam;
# Remove the useless header
tail -n +2 gene_counts.txt > gene_counts.no_header.txt;
mv gene_counts.txt gene_counts.no_header.txt ../counts/
mv gene_counts.txt.summary ../logs/
```

### Conversion to wig for TRANSIT

```bash
cd ../aln/
for bam in *.bam; do
    # 1. Convert bam to BED
    # 2. Reduce read size
    # 3. Sort
    echo "variableStep chrom=NC_000964.3" > ${bam%.bam}.wig;
    bedtools bamtobed -i ${bam} \
    | bedtools slop \
        -i - \
        -g ../ref/B_subtilis_subtilis_s168.genome \
        -l 0 -r -12 -s \
    | bedtools sort -i - \
    | bedtools intersect \
        -c \
        -sorted \
        -a - \
        -b ../ref/B_subtilis_subtilis_s168.TA.bed \
    | bedtools genomecov \
        -i - \
        -g ../ref/B_subtilis_subtilis_s168.genome \
        -bga -5 \
    | cut -f2,4 >> ${bam%.bam}.wig
done;
```

## QC report

Finally, we take all the processing steps's logs, and combine them together into a single QC report. MultiQC is a python program that is designed to read a multitude of bioinformatics software's output logs to summarize in a single HTML report.

```bash
cd ../logs/
multiqc .
```