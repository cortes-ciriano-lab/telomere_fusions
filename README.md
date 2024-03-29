# TelFusDetector: detection of telomere fusions using sequencing data

**IMPORTANT NOTICE:** This GitHub repository presents the old version of *TelFusDetector*. For an upgraded and much faster version of this tool, entirely coded in *Python* and executable through a single command, please look at the latest release of *TelFusDetector* available in the [new version.](https://github.com/cortes-ciriano-lab/TelFusDetector) 

TelFusDetector provides functionalities for the detection of telomere fusions using whole-genome sequencing data.
The pipeline, which is split into the 6 steps described below, requires an aligned bam file as input and a sample/run ID, which is used as prefix in the output files. For illustration purposes, we use below "test" as sample ID. 
We provide a test bam file ("test.bam") with the repository for testing purposes.<br>
The code can be run serially. However, for large files we recommend to run steps 2 and 3 below in parallel using subsets of the data as input (e.g. performing the analyses on a per chromosome basis).<br>

# License
**TelFusDetector is free for academic use only.** If you are not a member of a public funded academic and/or education and/or research institution you must obtain a commercial license from EMBL Enterprise Management GmbH (EMBLEM); please email EMBLEM (info@embl-em.de).

# Installation and requirements
The code requires samtools (https://github.com/samtools/samtools), python version >=3.7.0 and R version >=3.5.0.

Please install the dependencies required by running:<br>
```python
pip install -r requirements.txt
```
```R
Rscript r_requirements_install.R
```
<br>

## Step 1: extracting sequencing coverage information

The first step consists of extracting read alignment information using samtools flagstat.<br>
```
sampleid=test
./scripts/coverage_info.sh test.bam ${sampleid}.cov
```

Output file:<br>
Read alignemnt information, which will be used in Step 4. 

## Step 2: extracting reads with telomeric repeats (TTAGGG and CCCTAA) allowing for one mismatch

From a bam file input via stdin<br>
```python
samtools view test.bam | python scripts/fusion_caller.py --mode callfusions  --outprefix ${sampleid}
```
Output files:<br>
test_fusions: file containing reads with candidate telomere fusions<br>
test_readIDs: file containing the IDs for the reads with candidate telomere fusions<br>
<br>
If running in an HPC cluster, one might want to split the computation by e.g., chromosome. One can do that by running, for example:<br>

```bash
for chr in {1..22} X Y
do
bsub -o o.log -e e.log -J ${sampleid}job -M 2G "samtools view test.bam ${chr} | python scripts/fusion_caller.py --mode callfusions  --outprefix ${sampleid}_${chr}"
done
# we also process unmapped reads:
bsub -o o.log -e e.log -J ${sampleid}job -M 2G "samtools view -f 4 test.bam | python scripts/fusion_caller.py --mode callfusions  --outprefix ${sampleid}_unmapped"
```

When running in parallel, merge the read IDs into a single file:<br>
```bash
touch ${sampleid}_readIDs
for chr in {1..22} X Y unmapped
do
cat ${sampleid}_${chr}_readIDs >> ${sampleid}_readIDs
done
```
<br>

## Step 3: extracting the mates for the reads containing candidate telomere fusions
```python
samtools view test.bam | python scripts/fusion_caller.py --mode extractmates --outprefix ${sampleid}
```
When splitting the calculation by chromosome in step 1, pass as input the file containing the IDs for all reads across all chromsomes with potential telomere fusions.<br>
```python
samtools view test.bam | python scripts/fusion_caller.py --mode extractmates --outprefix ${sampleid} --readIDs ${sampleid}_readIDs
```
Output file:<br>
test_mates: file containing the mate reads for the reads with candidate fusions identified in step 1.<br>
<br>

Step 3 can also be run on a per chromosome basis as follows:
```bash
for chr in {1..22} X Y
do
bsub -o o.log -e e.log -J ${sampleid}job -M 2G "samtools view test.bam ${chr} | python scripts/fusion_caller.py --mode extractmates  --outprefix ${sampleid}_${chr}"
done
# we also process unmapped reads:
bsub -o o.log -e e.log -J ${sampleid}job -M 2G "samtools view -f 4 test.bam | python scripts/fusion_caller.py --mode extractmates  --outprefix ${sampleid}_unmapped"
```

When running in parallel, merge the read IDs and fusions into a single file:<br>
```bash
touch ${sampleid}_mates
for chr in {1..22} X Y unmapped
do
cat ${sampleid}_${chr}_mates >> ${sampleid}_mates
done

touch ${sampleid}_fusions
for chr in {1..22} X Y unmapped
do
cat ${sampleid}_${chr}_fusions >> ${sampleid}_fusions
done

```
<br>

## Step 4: generating a summary file containing the fusions and additional information
```python
python scripts/fusion_caller.py --mode summarise --outprefix ${sampleid} --alignmentinfo ${sampleid}.cov
```
Output file:<br>
test_fusions_summary: file reporting candidate fusions, and information about both reads the pair.<br>
<br>
Note that when calling the summarise mode, the programme expects to find in the running directory two files, namely:<br>
test_fusions<br>
test_mates<br>
Note that these two files are expected to start with the prefix passed to the option "--outprefix"<br>
<br>
Alternative, one can pass the two files as input as follows:<br>
```python
python scripts/fusion_caller.py --mode summarise --outprefix ${sampleid} --matesfile ${sampleid}_mates --fusionsfile ${sampleid}_fusions
```
<br>

## Step 5: quality control (QC) step 
This script provides functionalities to filter potential false positive telomere fusions, flag reads mapping to regions of the human genome containing telomere fusion-like patterns (such as the relic of an ancestral telomere fusion in chr2), and further classify the fusions detected into categories according to the patterns of repeats detected.<br>

```R
Rscript scripts/FusionReadsQC.R --summary_file ${sampleid}_fusions_summary --ref_genome Hg38 --outprefix QC/${sampleid} --read_length 150 
```
Output file:<br>
- QC/test.fusions.unfiltered.tsv: Updated summary file, with extra columns with the information used to perform the QC computation. It provides the QC decision for each read, as well as the reason underlying each decision.<br>
- QC/test.fusions.pass.tsv: List of reads supporting telomere fusions that passed all QC filters (PASS reads from QC/test.fusions.unfiltered.tsv file).<br>
- QC/test.fusions.false_positives.tsv: List of reads supporting telomere fusions that did not pass the QC filters (Filtered reads from QC/test.fusions.unfiltered.tsv file).<br>
- QC/test.fusions.pass.collapsed.tsv: List of reads supporting telomere fusions (only PASS reads) collapsed by chromosome, breakpoint sequence and fusion subtype.<br>
- QC/test.fusions.sample_stats.tsv: Summary table reporting the number of reads supporting fusions found for each sample, as well as the reads that were filtered out.<br>
- QCtest.fusions.QC.Rdata: R environment used in the computation (to be ignored by most users).<br>


## Step 6: breakpoint sequence correction
This script corrects the breakpoint sequences of the telomere fusions detected. The breakpoint sequence of a telomere fusion is the sequence flanked by the forward (TTAGGG) and reverse (CCCTAA) repeats.<br>

```R
Rscript scripts/CollapseCorrectFusions.R --summary_file_collapsed QC/${sampleid}.fusions.pass.collapsed.tsv --outprefix Collapsed_results/${sampleid}
```
Output file:<br>
- Collapsed_results/Possible_breakpoint_sequences.pure.tsv: All possible breakpoint sequences that can be originated from the canonical fusion of two telomeres.<br>
- Collapsed_results/test.breakpoint_correction_steps.tsv: Different steps of the error correction component showing how the breakpoint sequence of each fusion has been curated.<br>
- Collapsed_results/test.corrected.tsv: Telomere fusions obtained after breakpoint sequence correction. All fusions are collapsed by chromosome, breakpoint sequence, and the fusion subtype.<br>
- Collapsed_results/test.proportion_correct_endo9.tsv: Table reporting the proportion of reads annotated as endogenous_9 showing the expected TTAA breakpoint sequence.<br>

# Contact
If you have any comments or suggestions please raise an issue or contact us:
Francesc Muyas: fmuyas@ebi.ac.uk
Isidro Cortes-Ciriano: icortes@ebi.ac.uk
