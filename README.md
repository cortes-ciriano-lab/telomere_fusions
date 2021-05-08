# Step 1: extracting reads with telomere repeats (TTAGG) and inverted telomere repeats (CCCTAA) allowing for one mismatch

# From a bam file passed via stdin
samtools view test.bam | python fusion_caller.py --mode callfusions  --outprefix test

# Step 2: extract the mates for the reads containing candidate telomere fusions
samtools view test.bam | python fusion_caller.py --mode extractmates --outprefix test

# Step 3: generate a summary file containing the fusions and additional information
python fusion_caller.py --mode summarise --outprefix test
