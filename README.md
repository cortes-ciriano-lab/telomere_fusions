## Step 1: extracting reads with telomere repeats (TTAGGG) and inverted telomere repeats (CCCTAA) allowing for one mismatch

From a bam file passed via stdin<br>
```python
samtools view test.bam | python fusion_caller.py --mode callfusions  --outprefix test
```
Output files:<br>
test_fusions: file containing reads with candidate telomere fusions<br>
test_readIDs: file containing the IDs for the reads with candidate telomere fusions<br>
<br>
If running in an HPC cluster, one might want to split the computation by e.g., chromosome. One can do that by running, for example:<br>

```bash
for chr in {1..22} X Y
do
bsub -o o.log -e e.log -J testjob -M 2G "samtools view test.bam ${chr} | python fusion_caller.py --mode callfusions  --outprefix test_${chr}"
done
# we also process unmapped reads:
bsub -o o.log -e e.log -J testjob -M 2G "samtools view -f 4 test.bam | python fusion_caller.py --mode callfusions  --outprefix test_unmapped"
```

When running in parallel, merge the read IDs into a single file:<br>
```bash
touch test_readIDs
for chr in {1..22} X Y unmapped
do
cat test_${chr}_readIDs >> test_readIDs
done
```
<br>

## Step 2: extract the mates for the reads containing candidate telomere fusions
```python
samtools view test.bam | python fusion_caller.py --mode extractmates --outprefix test
```
When splitting the calculation by chromosome in step 1, pass as input the file containing the IDs for all reads across all chromsomes with potential telomere fusions.<br>
```python
samtools view test.bam | python fusion_caller.py --mode extractmates --outprefix test --readIDs test_readIDs
```
Output file:<br>
test_mates: file containing the mate reads for the reads with candidate fusions identified in step 1.<br>
<br>

Step 2 can also be run on a per chromosome basis as follows:
```bash
for chr in {1..22} X Y
do
bsub -o o.log -e e.log -J testjob -M 2G "samtools view test.bam ${chr} | python fusion_caller.py --mode extractmates  --outprefix test_${chr}"
done
# we also process unmapped reads:
bsub -o o.log -e e.log -J testjob -M 2G "samtools view -f 4 test.bam | python fusion_caller.py --mode extractmates  --outprefix test_unmapped"
```

When running in parallel, merge the read IDs into a single file:<br>
```bash
touch test_mates
for chr in {1..22} X Y unmapped
do
cat test_${chr}_mates >> test_mates
done
```
<br>

## Step 3: generate a summary file containing the fusions and additional information
```python
python fusion_caller.py --mode summarise --outprefix test
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
python fusion_caller.py --mode summarise --outprefix test --matesfile test_mates --fusionsfile test_fusions
```
<br>

## Step 4: quality control (QC) step 
This script aims to filter potential false positives telomere fusions, mark endogenous regions of human genome prone to be confused by telomere fusions and check divide the fusions in different subtypes.<br>

```R
Rscript FusionReadsQC.R --summary_file test_fusions_summary --ref_genome Hg38 --project test --prefix QC/test
```
Output file:<br>
QC/test.fusions.unfiltered.tsv: Updated summary file, with extra columns with the information used to perform the QC computation. It provides the QC decision for each read, as well as the reason because each read was filtered or not.<br>
QC/test.fusions.pass.tsv: List of reads supporting telomere fusions that passed all QC filters (PASS reads from QC/test.fusions.unfiltered.tsv file).<br>
QC/test.fusions.false_positives.tsv: List of reads supporting that did not pass the QC filters (Filtered reads from QC/test.fusions.unfiltered.tsv file).<br>
QC/test.fusions.pass.collapsed.tsv: List of reads supporting telomere fusions (only PASS reads) collapsed by chromosome, middle sequence, and the different fusion subtype criteria for each sample separately.<br>
QC/test.fusions.sample_stats.tsv: Summary table showing the number of reads supporting fusions found for each samples, as well as the reads of them that were filtered out.<br>
QC/test.fusions.project_stats.tsv: Summary table showing the number of reads supporting fusions found for each project group, as well as the number of reads that were filtered out.<br>
QCtest.fusions.QC.Rdata: R environment used in the computation (to be ignored by most of users).<br>


## Step 5: middle sequence correction
This script aims to correct the middle sequences of the telomere fusions to better distinguish the juntion point between the forward (TTAGGG) and reverse (CCCTAA) repeats of the telomeres fused.<br>

```R
Rscript CollapseCorrectFusions.R --summary_file_collapsed QC/test.fusions.pass.collapsed.tsv --prefix Collapsed_results/test
```
Output file:<br>
Collapsed_results/Possible_middle_sequences.pure.tsv: All possible middle sequences that can be originated from the canonical fusion of two telomeres.<br>
Collapsed_results/test.middle_correction_steps.tsv: Different steps of the error correction showing how the middle sequence of each fusion has been curated.<br>
Collapsed_results/test.corrected.tsv: Final version of the telomere fusions obtained after middle correction. All fusions are collapsed by chromosome, middle sequence, and the different fusion subtype criteria (for each sample separately).<br>
Collapsed_results/test.blacklisted_samples.tsv: Sample IDs listed in the blacklist based on the proportion of endogenous 9 fusion-like events showing the middle sequence TTAA.<br>
Collapsed_results/test.proportion_correct_endo9.all_samples.pdf: Distribution of reads annoated in endogenous_9 showing the expected TTAA middle sequence (all samples).<br>
Collapsed_results/test.proportion_correct_endo9.top100.pdf: Distribution of reads annoated in endogenous_9 showing the expected TTAA middle sequence (top 100 samples).<br>

## Requirements

Please install the dependencies required by running:<br>
```python
pip install -r requirements.txt
```
```R
Rscript r_requirements_install.R
```

