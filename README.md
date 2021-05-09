## Step 1: extracting reads with telomere repeats (TTAGG) and inverted telomere repeats (CCCTAA) allowing for one mismatch

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

## Requirements

Please install the dependencies required by running:<br>
```python
pip install -r requirements.txt
```


