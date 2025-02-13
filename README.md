# EXTRA-seq
Material and source code for the EXTRA-seq manuscript

## Citation
Currently in preprint at BioRXiv:<br/>
**EXTRA-seq: a genome-integrated extended massively parallel reporter assay to quantify enhancer-promoter communication**<br/>
Judith F. Kribelbauer-Swietek,  Vincent Gardeux,  Gerard Llimos-Aubach, Katerina Faltejskova, Julie Russeil, Nadia Grenningloh, Lucas Levassor, Clémence Steiner,  Jiri Vondrasek,  Bart Deplancke<br/>
DOI: [https://doi.org/10.1101/2024.12.08.627402](https://doi.org/10.1101/2024.12.08.627402)

## 0. Raw sequencing data
All sequencing files were deposited on ArrayExpress under the accession number [E-MTAB-14800](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14800). It comprises raw **.fastq** files of the different sequencing protocols (gDNA, mRNA, Nanopore, and STARR-seq), as well as **.fasta** files of the construct sequences, and processed count matrices as **.csv** files for 9 different libraries (3.5, 3.8, 3.9, 4.8, 4.9, 5.1, 5.6, 5.7, 5.8).<br/>

## 1. Preprocessing of the Nanopore long-read sequencing data

### 1.1 Basecalling using guppy

The protocol used Oxford Nanopore Technologies (ONT) barcodes (SQK-NBD114-24) and adaptor sequences for ONT sequencing. Libraries were sequenced on a Minion or a PromethION 2 Solo device using MINFlow114 and FLO-PRO114M flowcells respectively. The raw fast5 and/or pod5 data files were processed using the guppy basecaller v6.5.7 provided by ONT on a GPU server.

```bash
guppy_basecaller -i ${root_dir}/pod5 -s ${root_dir}/guppy_out --flowcell FLO-PRO114M --kit SQK-NBD114-24 -x "cuda:0" --barcode_kits "SQK-NBD114-24" --compress_fastq --gpu_runners_per_device 100 --chunks_per_runner 128
```

Then, for each barcode, we merged the sequences into a unique fastq.gz file:
```bash
mkdir ${root_dir}/fastq/
for bc in barcode01 barcode02 barcode03 barcode04 barcode05 barcode06 barcode07 barcode08 barcode09 barcode10 barcode11 barcode12 barcode13 barcode14 barcode15 barcode16 barcode17 barcode18 barcode19 barcode20 barcode21 barcode22 barcode23 barcode24 unclassified
do
	cat ${root_dir}/guppy_out/pass/${bc}/*.fastq.gz > ${root_dir}/fastq/${bc}.fastq.gz
done
```

### 1.2 Alignment using minimap2

Once we had our .fastq files. We kept the ones matching to our barcodes of interest. We then aligned these on the construct sequences using minimap2, and generated additional QCs:

```bash
# Parameters
barcode=${1}
nthreads=12

# Folders
mkdir -p bam/
mkdir -p qc/
mkdir -p bigwigs/
	
echo "Processing barcode: $barcode"

# Aligning to the construct
minimap2 -ax map-ont --secondary=no fasta/lib6p0_plasmid_seqeunces_for_nanopore_alignment_2024_12_02.fasta fastq/${barcode}.fastq.gz -o bam/${barcode}.sam

# Change sam to bam 
samtools view -S -b bam/${barcode}.sam > bam/${barcode}.bam
rm bam/${barcode}.sam

# Sort bam file 
samtools sort bam/${barcode}.bam -o bam/${barcode}.sorted.bam
mv bam/${barcode}.sorted.bam bam/${barcode}.bam

# Generate bam index
samtools index bam/${barcode}.bam

# QC
fastqc -t ${nthreads} --memory 10000 --outdir=qc/ fastq/${barcode}.fastq.gz
fastqc -t ${nthreads} --memory 10000 --outdir=qc/ bam/${barcode}.bam
samtools flagstat bam/${barcode}.bam > qc/${barcode}.flagstat.txt
samtools idxstats bam/${barcode}.bam > qc/${barcode}.idxstats.txt

# BigWig
bamCoverage -p ${nthreads} --bam bam/${barcode}.bam -o bigwigs/${barcode}.bw --binSize 10 --normalizeUsing RPKM --outFileFormat "bigwig"

# Remove multiple mapped reads and repeat QC
samtools view -b -F0x900 bam/${barcode}.bam > bam/${barcode}.nomult.bam
samtools index bam/${barcode}.nomult.bam
samtools idxstats bam/${barcode}.nomult.bam  > qc/${barcode}.nomult.idxstats.txt 
samtools flagstat bam/${barcode}.nomult.bam > qc/${barcode}.nomult.flagstat.txt
bamCoverage -p ${nthreads} --bam bam/${barcode}.nomult.bam -o bigwigs/${barcode}.nomult.bw --binSize 10 --normalizeUsing RPKM --outFileFormat "bigwig"
```

### 1.3 Counting barcode occurrences in long-reads

Then, we used homemade scripts ([./software/MinionBarcodes-1.0.jar](./software/MinionBarcodes-1.0.jar)) to extract the barcode sequences at a given position (see the .csv files deposited on ArrayExpress which describes, for each .fasta sequence construct, the exact location of the barcode sequences).
```bash
java -jar ./software/MinionBarcodes-1.0.jar Counter --bam bam/barcode23.nomult.bam --config csv/6p0_library_enh_seq_BC_extraction_coordinates_hublib4p75_for_NanoP_alignment_2024_12_02_updated.csv --mapq 0 -o output/barcode23.nomult.tsv >output/barcode23.nomult.log
```

### 1.4 running Enformer predictions

We provide the .yml file to set up a conda environment including all dependencies used to run Enformer on custom seqeunces. Within the ([./software/Enformer/] directory, there are two markdowns that provide example code on 

i) how to run Enformer on genomic seqeunces modified to reflect EXTRA-seq inserted seqeunces &

ii) how we processed downstream predictions to derive summarized predictions scores

