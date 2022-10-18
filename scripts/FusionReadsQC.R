#---------------------
# Script to perform the QC in the summary file
#---------------------

#-------
# Packages required
#-------
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(stringdist))
options(warn=-1)
#-------
# Function to run the QC fusion detection step
#-------

# Function to detect the number of breakpoints in the middle sequence
NUM_BREAKPOINTS <- function(SEQ, orientation){
  
  # Replace known repeats
  SEQ <- gsub(pattern = 'TTAGGG', replacement = '__fw__',x = SEQ)
  SEQ <- gsub(pattern = 'CCCTAA', replacement = '__rv__',x = SEQ)
  
  # Replace old substitutions
  SEQ <- gsub(pattern = '_forward_', replacement = '__fw__',x = SEQ)
  SEQ <- gsub(pattern = '_reverse_', replacement = '__rv__',x = SEQ)
  
  ALL_ISLANDS <- unlist(strsplit(SEQ, split = '__'))
  ALL_ISLANDS <- ALL_ISLANDS[!ALL_ISLANDS %in%  c('','fw','rv')]
  
  
  if (orientation == 'outward'){
    SEQ2 <- paste('__rv__',SEQ,'__fw__',sep = '')
    ISLAND <- unlist(str_extract_all(string = SEQ2, pattern = "__rv__[A-Z]*__fw__"))
  } else {
    SEQ2 <- paste('__fw__',SEQ,'__rv__',sep = '')
    ISLAND <- unlist(str_extract_all(string = SEQ2, pattern = "__fw__[A-Z]*__rv__"))
  }
  
  return(length(ISLAND))
}

# Function to analyse mate read
# It returns the number of times the middle in read 1 is found in the mate
search_middle_in_mate_one_mismatch <- function(middle, seq, seq_mate, orientation){
  
  # Get original middle sequence
  middle <- gsub(pattern = "_reverse_", replacement = 'CCCTAA', x = middle)
  middle <- gsub(pattern = "_forward_", replacement = 'TTAGGG', x = middle)
  
  # Extend middle sequence
  # 6 bp up and downstream
  if (orientation == 'outward'){
    middle_pat <- paste0("CCCTAACCCTAA",middle,"TTAGGGTTAGGG")
  } else {
    middle_pat <- paste0("TTAGGGTTAGGG",middle,"CCCTAACCCTAA")
  }
  
  # Search for the original middle sequence
  middle2 <- matchPattern(pattern = middle_pat, subject = seq,
                          max.mismatch=8, min.mismatch=0,
                          with.indels=FALSE, fixed=TRUE,
                          algorithm="auto")
  
  # In case the match is found several times in the read
  # We filter and clean them
  if (length(middle2) > 1) {
    middle2 <- middle2[nchar(middle2) == nchar(middle_pat)]
    middle2 <- middle2[substr(middle2, 13, nchar(middle_pat)-12) == middle]
    # Check each one of the two hexamers up- and down-stream of the middle sequence
    if (orientation == 'outward'){
      middle2 <- middle2[stringdist::stringdist(a = 'CCCTAA', b =substr(middle2, 1, 6),method = "hamming") <= 2 &
                           stringdist::stringdist(a = 'CCCTAA', b = substr(middle2, 7, 12),method = "hamming") <= 2 &
                           stringdist::stringdist(a = 'TTAGGG', b = substr(middle2, nchar(middle_pat)-12+1, nchar(middle_pat)-6),method = "hamming") <= 2 &
                           stringdist::stringdist(a = 'TTAGGG', b = substr(middle2, nchar(middle_pat)-6+1, nchar(middle_pat)),method = "hamming") <= 2]
    } else {
      middle2 <- middle2[stringdist::stringdist(a = 'TTAGGG', b =substr(middle2, 1, 6),method = "hamming") <= 2 &
                           stringdist::stringdist(a = 'TTAGGG', b = substr(middle2, 7, 12),method = "hamming") <= 2 &
                           stringdist::stringdist(a = 'CCCTAA', b = substr(middle2, nchar(middle_pat)-12+1, nchar(middle_pat)-6),method = "hamming") <= 2 &
                           stringdist::stringdist(a = 'CCCTAA', b = substr(middle2, nchar(middle_pat)-6+1, nchar(middle_pat)),method = "hamming") <= 2]
    }
  }
  
  # Let's search in the mate if we found the sequence of the middle properly
  if (length(middle2) == 1){
    SEQ <- as.character(middle2[[1]])
    SEQ <- substr(SEQ, 7, nchar(middle_pat)-6)
    
    # Reverse complement the middle sequence
    REV <- as.character(reverseComplement(DNAString(SEQ)))
    
    # Search the middle sequence in the mate using + and - strand
    x <- DNAString(seq_mate)
    Count1 <- countPattern(pattern = SEQ, subject = x, max.mismatch = 1, with.indels = T)
    Count2 <- countPattern(pattern = REV, subject = x, max.mismatch = 1, with.indels = T)
    
    # Sum up the counts in one value
    Count <- max(Count1,Count2)
  } else if (length(middle2)  == 0){
    Count <- -1
  } else {
    Count <- length(middle2) * (-1)
  }
  return(Count)
}

# Function to look for inversion-like fusions
# It considers the order of the read (R1 or R2), as well as if it has been reverse complemented or not
FUSION_INVERSION <- function(orientation,FORWARD_mate,REVERSE_mate,FLAG){
  CHECK1 <- 0
  CHECK2 <- 0
  if (orientation == 'outward'){
    if (str_count(FLAG, pattern = 'Minus') == 1){
      if (str_count(FLAG, pattern = 'isSecondMateRead') == 1 & grepl(pattern = "isMateMinusStrand", x = FLAG)){
        CHECK1 <- FORWARD_mate
        CHECK2 <- REVERSE_mate
      } else if (str_count(FLAG, pattern = 'isSecondMateRead') == 1 & !grepl(pattern = "isMateMinusStrand", x = FLAG)){
        CHECK1 <- REVERSE_mate
        CHECK2 <- FORWARD_mate
      } else if (str_count(FLAG, pattern = 'isFirstMateRead') == 1 & grepl(pattern = "isMateMinusStrand", x = FLAG)){
        CHECK1 <- FORWARD_mate
        CHECK2 <- REVERSE_mate
      } else if (str_count(FLAG, pattern = 'isFirstMateRead') == 1 & !grepl(pattern = "isMateMinusStrand", x = FLAG)){
        CHECK1 <- REVERSE_mate
        CHECK2 <- FORWARD_mate
      }
    } else if (str_count(FLAG, pattern = 'Minus') > 0){
      if (str_count(FLAG, pattern = 'isSecondMateRead') == 1){
        CHECK1 <- FORWARD_mate
        CHECK2 <- REVERSE_mate
      } else if (str_count(FLAG, pattern = 'isFirstMateRead') == 1){
        CHECK1 <- FORWARD_mate
        CHECK2 <- REVERSE_mate
      }
    } else {
      if (str_count(FLAG, pattern = 'isSecondMateRead') == 1){
        CHECK1 <- REVERSE_mate
        CHECK2 <- FORWARD_mate
      } else if (str_count(FLAG, pattern = 'isFirstMateRead') == 1){
        CHECK1 <- REVERSE_mate
        CHECK2 <- FORWARD_mate
      }
    }
  } else{
    # End-to-head 
    if (str_count(FLAG, pattern = 'Minus') == 1){
      if (str_count(FLAG, pattern = 'isSecondMateRead') == 1 & grepl(pattern = "isMateMinusStrand", x = FLAG)){
        CHECK1 <- REVERSE_mate
        CHECK2 <- FORWARD_mate
      } else if (str_count(FLAG, pattern = 'isSecondMateRead') == 1 & !grepl(pattern = "isMateMinusStrand", x = FLAG)){
        CHECK1 <- FORWARD_mate
        CHECK2 <- REVERSE_mate
      } else if (str_count(FLAG, pattern = 'isFirstMateRead') == 1 & grepl(pattern = "isMateMinusStrand", x = FLAG)){
        CHECK1 <- REVERSE_mate
        CHECK2 <- FORWARD_mate
      } else if (str_count(FLAG, pattern = 'isFirstMateRead') == 1 & !grepl(pattern = "isMateMinusStrand", x = FLAG)){
        CHECK1 <- FORWARD_mate
        CHECK2 <- REVERSE_mate
      }
    } else if (str_count(FLAG, pattern = 'Minus') > 0){
      if (str_count(FLAG, pattern = 'isSecondMateRead') == 1){
        CHECK1 <- REVERSE_mate
        CHECK2 <- FORWARD_mate
      } else if (str_count(FLAG, pattern = 'isFirstMateRead') == 1){
        CHECK1 <- REVERSE_mate
        CHECK2 <- FORWARD_mate
      }
    } else {
      if (str_count(FLAG, pattern = 'isSecondMateRead') == 1){
        CHECK1 <- FORWARD_mate
        CHECK2 <- REVERSE_mate
      } else if (str_count(FLAG, pattern = 'isFirstMateRead') == 1){
        CHECK1 <- FORWARD_mate
        CHECK2 <- REVERSE_mate
      }
    }
  }
  # Classifier
  if (CHECK1 > 0 & CHECK2 < 1){
    TYPE <- 'Fusion-like'
  } else if (CHECK1 < 1 & CHECK2 > 0){
    TYPE <- 'Inversion-like'
  } else {
    TYPE <- 'Other'
  }
  return(TYPE)
}

# Function to mark endogenous fusion-like regions of the human genome
# Expected coordinates are passed in the for loop where this function is called
MARK_ENDO_REGION <- function(DATA,CHROM,MIN_COORDINATE,MAX_COORDINATE,Flag_endo,orientation){
  # Flag the endogenous sequence in read of fusion
  DATA$chr =gsub("chr","",as.vector(DATA$chr))
  DATA$chr <- as.character(DATA$chr)
  if (Flag_endo == '2_endogenous'){
    idx1 = which(DATA$chr == CHROM & 
                   DATA$pos >= MIN_COORDINATE & DATA$pos <= MAX_COORDINATE & 
                   DATA$MAPQ >= 8 & DATA$orientation %in% orientation)
    idx2 = which((grepl(pattern = 'TTGGGGTTGGGG', x = DATA$SEQ) | grepl(pattern = 'CCCCAACCCCAA', x = DATA$SEQ)) &
                   grepl(pattern = 'TTAGCTAA', x = DATA$SEQ))
    
    idx = unique(c(idx1,idx2))
  } else {
    idx = which(DATA$chr == CHROM & 
                  DATA$pos >= MIN_COORDINATE & DATA$pos <= MAX_COORDINATE & 
                  DATA$MAPQ >= 8 & DATA$orientation %in% orientation)
  }
  DATA$chr[idx] = Flag_endo
  
  # Flag the endogenous sequence in read mate
  DATA$chr2 =gsub("chr","",as.vector(DATA$chr2))
  DATA$chr2 <- as.character(DATA$chr2)
  idx = which(DATA$chr2 == CHROM & 
                DATA$PNEXT >= MIN_COORDINATE & DATA$PNEXT <= MAX_COORDINATE & 
                DATA$MAPQ_mate >= 8 & DATA$orientation %in% orientation)
  
  DATA$chr2[idx] = Flag_endo
  
  return(DATA)
}



# Main function to run Quality Control (QC) step with summary files 
MAIN_FUNCTION_FUSION_QC <- function(file, project,REFERENCE, read_length){
  # Initial print
  cat("\n")
  print("___")
  print(paste0('Project: ', project))
  
  # Load project file and add header
  d <- tryCatch(read.csv(file, sep = '\t', header = F), error=function(e) NULL)

  # Check if empty files
  if (is.null(d)){
    print(paste0('Warning: ',project, ' files are empty or not readable for the tool'))
    return(list(fusions = NULL, stats = NULL))
  }
  
  # Check if empty files
  if (nrow(d) < 1){
    print(paste0('Warning: ',project, ' files are empty or not readable for the tool'))
    return(list(fusions = NULL, stats = NULL))
  }

  # Add colnames                
  if (nrow(d) > 0){
    names(d) <- HEADER
  }
                
  d <- unique(d)
  
  # Remove wrong reads
  d$MAPQ <- as.numeric(d$MAPQ)
  d$MAPQ_mate <- as.numeric(d$MAPQ_mate)
  d$pos <- as.numeric(d$pos)
  d$pos2 <- as.numeric(d$pos2)
  d$PNEXT <- as.numeric(d$PNEXT)
  
  d <- d[!is.na(d$MAPQ) & !is.na(d$MAPQ_mate),]
  
  if (nrow(d) < 1){
    print(paste0(project, ' files are empty or not readable for the tool'))
    return(list(fusions = NULL, stats = NULL))
  }
  
  # Clean sample name
  d$sample = gsub( "_rawFus_full_chr[0-9]_python.fq", "", d$sample)
  d$sample = gsub( "_rawFus_full_chrX_python.fq", "", d$sample)
  d$sample = gsub( "_rawFus_full_chrY_python.fq", "", d$sample)
  d$sample = gsub( "_rawFus_full_chr[(0-9)][(0-9)]_python.fq", "", d$sample)
  d$sample = gsub( "_rawFus_unmapped_python.fq", "", d$sample)
  d$sample = gsub( ".sorted", "", d$sample)
  d$sample = gsub( ".sam", "", d$sample)
  d$sample <- gsub(pattern = '.final.cram',replacement = '', x = d$sample)
  d$sample <- gsub(pattern = '.recal.md',replacement = '', x = d$sample)
  d$sample <- gsub(pattern = '.MD',replacement = '', x = d$sample)
  d$sample <- gsub(pattern = '.cram',replacement = '', x = d$sample)
  
  d$project <- project
  
  # Fix NA in forward and reverse counts
  if (sum(is.na(d$forward_mate)) > 0){
    d[is.na(d$forward_mate),]$forward_mate <- 0
  }
  
  if (sum(is.na(d$reverse_mate)) > 0){
    d[is.na(d$reverse_mate),]$reverse_mate <- 0
  }
  
  ##-----------------------
  # Clean middle sequences from the beginning and the end
  ##-----------------------
  d$middle <- gsub(pattern = "^[_forward_|_reverse_]+", replacement = '', x = d$middle, perl = F)
  d$middle <- gsub(pattern = "([_forward_|_reverse_]+)$", replacement = '', x = d$middle, perl = F)
  
  ##-----------------------
  # Checking mismatches in fusions
  # Fusions with mismatches
  ##-----------------------
  
  idx= which( (d$forward0  == 0 | d$reverse0 == 0) & ( (d$forward0 > 0 & d$reverse1>0) | (d$forward1 > 0 & d$reverse0>0) ) )
  
  # Fusions without mismatches
  idx2= which( (d$forward0 > 0 & d$reverse0>0))
  
  # Change labels
  dd1 = d[idx2,]
  dd1$label = "No mismatches"
  
  if(length(idx) > 0){
    dd2 = d[idx,]
    dd2$label = "One mismatch"
    dd=rbind(dd1,dd2)
  }else{
    dd=dd1
  }
  
  ##-----------------------
  # Mark endogenous regions
  # endo chr2
  # endo chr4
  # endo chr9
  ##-----------------------
  
  # Endogenous 2
  if (REFERENCE == 'Hg19'){
    # Hg19 coordinates
    CHROM <- "2"
    MIN_COORDINATE <- 114355250
    MAX_COORDINATE <- 114365750
    Flag_endo <- "2_endogenous"
    orientation <- c("end-to-head")
  } else {
    # Hg38 coordinates
    CHROM <- "2"
    MIN_COORDINATE <- 113597750
    MAX_COORDINATE <- 113608250
    Flag_endo <- "2_endogenous"
    orientation <- c("end-to-head")
  }
  
  dd <- MARK_ENDO_REGION(dd,CHROM,MIN_COORDINATE,MAX_COORDINATE,Flag_endo,orientation)
  
  # Endogenous 4
  if (REFERENCE == 'Hg19'){
    # Hg19 coordinates
    CHROM <- "j"
    MIN_COORDINATE <- 0
    MAX_COORDINATE <- 0
    Flag_endo <- "4_endogenous"
    orientation <- c("end-to-head","outward")
  } else {
    # Hg38 coordinates
    CHROM <- "4"
    MIN_COORDINATE <- 190177250
    MAX_COORDINATE <- 190188000
    Flag_endo <- "4_endogenous"
    orientation <- c("end-to-head","outward")
  }
  
  dd <- MARK_ENDO_REGION(dd,CHROM,MIN_COORDINATE,MAX_COORDINATE,Flag_endo,orientation)
  
  # Endogenous 9
  if (REFERENCE == 'Hg19'){
    # Hg19 coordinates
    CHROM <- "9"
    MIN_COORDINATE <- 130911250
    MAX_COORDINATE <- 130922000
    Flag_endo <- "9_endogenous"
    orientation <- c("end-to-head")
  } else {
    # Hg38 coordinates
    CHROM <- "9"
    MIN_COORDINATE <- 128149000
    MAX_COORDINATE <- 128159750
    Flag_endo <- "9_endogenous"
    orientation <- c("end-to-head")
  }
  
  dd <- MARK_ENDO_REGION(dd,CHROM,MIN_COORDINATE,MAX_COORDINATE,Flag_endo,orientation)
  
  
  #-------------------
  # Substitute or mark chromosome, step 1
  #-------------------
  dd$chr_final <- apply(dd, 1, function(x) if(x['chr'] %in% c('2_endogenous','4_endogenous','9_endogenous')){
    x['chr']
  } else if (x['chr2'] %in% c('2_endogenous','4_endogenous','9_endogenous')){
    x['chr2']
  } else {
    x['chr']
  }
  )
  
  ##-----------------------
  # Translate flag code
  ##-----------------------
  
  dd$FLAG1 <- apply(dd, 1, function(x) paste(FLAG_BITNAMES[ as.logical( bamFlagAsBitMatrix(as.integer(x['flag'])))],sep = ',', collapse = ','))
  dd$CODE <- with(dd, paste(sample,chr, pos, nchar(SEQ),RNEXT, PNEXT,nchar(read2), CIGAR, sep = '-'))
  # Mark if read 1 (1 if yes)
  dd$R1 <- (grepl(pattern = 'isFirstMate', x = dd$FLAG1))*1
  
  #------------------------
  # Fusion-like or inversion-like
  #------------------------
  dd$Subtype <- apply(dd,1,function(x) if (!x['chr_final'] %in% c('2_endogenous','4_endogenous','9_endogenous')){
    if (as.numeric(x['reverse_mate']) == 0 & as.numeric(x['forward_mate']) == 0){
      'No_repeats_mate'
    } else {
      FUSION_INVERSION(as.character(x['orientation']),
                       as.numeric(x['forward_mate']),
                       as.numeric(x['reverse_mate']),
                       as.character(x['FLAG1']))
    }
  } else {
    x['chr_final']
  })
  
  ##-----------------------
  # Mark duplicates extracting the information from flag code
  # Only remove duplicate reads if there are more than one read mapping at the same coordinates
  ##-----------------------
  
  DUPLICATE_TABLE <- as.data.frame(with(dd,table(CODE)))
  DUPLICATE_TABLE <- DUPLICATE_TABLE[DUPLICATE_TABLE$Freq > 1,]
  
  # Mark if it is a duplicate read
  dd$DUPLICATE <- (grepl(pattern = 'isDuplicate',x = dd$FLAG1) & 
                     dd$CODE %in% DUPLICATE_TABLE$CODE)
  
  if (TRUE %in% dd$DUPLICATE){
    dd[dd$DUPLICATE == T,]$DUPLICATE <- 'isDuplicate'
  }
  
  if (FALSE %in% dd$DUPLICATE){
    dd[dd$DUPLICATE == F,]$DUPLICATE <- 'Correct'
  }
  
  CHECK <- TRUE %in% grepl(pattern = 'isDuplicate',x = dd$FLAG1) 
  
  
  ##-----------------------
  # Mark supplementary alignments
  # Mark secondary
  ##-----------------------
  if(nrow(dd) > 0){
    dd$Supplementary <- (grepl(pattern = 'isSupplementaryAlignment',x = dd$FLAG1))
    if (TRUE %in% dd$Supplementary){
      dd[dd$Supplementary == T,]$Supplementary <- 'isSupplementaryAlignment'
    }
    
    if (FALSE %in% dd$Supplementary){
      dd[dd$Supplementary == F,]$Supplementary <- 'Correct'
    }
  }
  
  if(nrow(dd) > 0){
    dd$Secondary <- (grepl(pattern = 'isSecondaryAlignment',x = dd$FLAG1))
    if (TRUE %in% dd$Secondary){
      dd[dd$Secondary == T,]$Secondary <- 'isSecondaryAlignment'
    }
    
    if (FALSE %in% dd$Secondary){
      dd[dd$Secondary == F,]$Secondary <- 'Correct'
    }
  }
  
  ##-----------------------
  # Mark reads which mate coordinates are wrong
  ##-----------------------
  if(nrow(dd) > 0){
    dd$proper_mate <- !is.na(dd$pos2) & !is.na(dd$pos) & !is.na(dd$PNEXT) & dd$pos2 == dd$PNEXT
    if (TRUE %in% dd$proper_mate){
      dd[dd$proper_mate == TRUE,]$proper_mate <- 'Correct_mate_found'
    }
    
    if (FALSE %in% dd$proper_mate){
      dd[dd$proper_mate == F,]$proper_mate <- 'Wrong_mate_found'
    }
  }
  
  ##-----------------------
  # Number of breakpoints per middle
  ##----------------------- 
  dd$Num_breakpoints <- apply(dd,1,function(x) NUM_BREAKPOINTS(x['middle'], x['orientation']))
  dd$Num_breakpoints <- paste0('Num_breakpoints_middle=',dd$Num_breakpoints)
  
  ##-----------------------
  # Remove false positives
  ##-----------------------
  
  # Check percentage of reads mapped per read
  dd$Bases_mapped_read_from_CIGAR <- cigarWidthAlongQuerySpace(dd$CIGAR,before.hard.clipping = T)
  dd$Bases_mapped_mate_from_CIGAR <- cigarWidthAlongQuerySpace(dd$CIGAR_mate,before.hard.clipping = T)
  dd$Prop_read_mapped <- nchar(dd$SEQ)/dd$Bases_mapped_read_from_CIGAR
  dd$Prop_mate_mapped <- nchar(dd$read2)/dd$Bases_mapped_mate_from_CIGAR
  
  # Get potential false positives
  # Not in endogenous regions, with good mapping into non-telomeric regions (middle of chromosomes)
  # We check afterwards if read wiith fusion maps to sub-telomeric regions
  idx = which(dd$MAPQ >= 8 & nchar(as.vector(dd$SEQ)) >= 60 & !dd$chr_final %in% c('2_endogenous','4_endogenous','9_endogenous'))
  
  # Type of fusions
  dd$type_fus1 = "fusion"
  
  if(length(idx) > 0){
    
    # Mark them as sub-telomeric 
    fp_now =  dd[idx,]
    fp_now$type_fus1="subtelomeric read"
    
    dd = dd[-idx,]
    
    # Now from the false positives, identify potential insertions of fusions, and fusions with the mate mapping to sub-telomeric regions
    # Which ones map close to the beginning of the chr?
    subtelo <- NULL
    for (ff in 1:nrow(fp_now)){
      if (REFERENCE == 'Hg38'){
        endtelonow = telo_hg38$start[match(fp_now$chr[ff], telo_hg38$chr)]
      } else{
        endtelonow = telo_hg19$start[match(fp_now$chr[ff], telo_hg19$chr)]
      }
      
      if( ((fp_now$pos[ff] < 100000) |  (fp_now$pos[ff] > endtelonow) ) & !is.na(endtelonow)){
        subtelo = rbind.fill(subtelo, fp_now[ff,])
        dd = rbind(dd, fp_now[ff,])
        
      } else{
        fp_now[ff,]$type_fus1 <- 'chromosomic read'
        dd = rbind(dd, fp_now[ff,])
      }
    }
  } 
  
  # Now in mate
  # Not in endogenous regions, with good mapping into non-telomeric regions (middle of chromosomes)
  # We check afterwards if mates map to sub-telomeric regions
  idx = which(dd$MAPQ_mate >= 8 & nchar(as.vector(dd$read2)) >= 60 & !dd$chr_final %in% c('2_endogenous','4_endogenous','9_endogenous'))
  
  # Type of fusions
  dd$type_fus2 = "fusion"
  
  if(length(idx) > 0){
    
    # Mark them as subtelomeric 
    fp_now =  dd[idx,]
    fp_now$type_fus2="subtelomeric mate"
    dd = dd[-idx,]
    
    # Now from the false positives, identify potential insertions of fusions, and fusions with the mate mapping to sub-telomeric regions
    # Which ones map close to the beginning of the chr?
    subtelo <- NULL
    for (ff in 1:nrow(fp_now)){
      if (REFERENCE == 'Hg38'){
        endtelonow = telo_hg38$start[match(fp_now$chr2[ff], telo_hg38$chr)]
      } else{
        endtelonow = telo_hg19$start[match(fp_now$chr2[ff], telo_hg19$chr)]
      }
      
      if( ((fp_now$pos2[ff] < 100000) |  (fp_now$pos2[ff] > endtelonow) ) & !is.na(endtelonow)){
        subtelo = rbind.fill(subtelo, fp_now[ff,])
        dd = rbind(dd, fp_now[ff,])
        
      } else{
        fp_now[ff,]$type_fus2 <- 'chromosomic mate'
        dd = rbind(dd, fp_now[ff,])
      }
    }
  }
  
  ##-----------------------
  # Remove reads with inconsistent mate mapping read
  # With low mapping quality
  # Long ones and without lots of hard-clip bases (more than 80 bps)
  # Without telomeric repeats
  # Not in endogenous
  # To remove fusions with unknown mate sequence without telomeric repeats
  ##-----------------------
  idx = which(dd$telo_repeats_mate == "No" & dd$MAPQ_mate < 8 & nchar(as.vector(dd$read2)) >= 60 & !dd$chr_final %in% c('2_endogenous','4_endogenous','9_endogenous') & dd$type_fus1 == 'fusion' & dd$type_fus2 == 'fusion')
  
  # Type of fusions
  dd$mate_mapping = "expected_mate_mapping"
  
  if (length(idx) > 0){
    fp_now =  dd[idx,]
    fp_now$mate_mapping="inconsistent_mate_mapping"
    dd = dd[-idx,]
    dd <- rbind(dd, fp_now)
  }
  
  ##-----------------------
  # Analyse mate reads
  # Check number of breakpoints in mate
  # Check if middle sequence in read 1 is found in mate
  ##-----------------------
  
  # Substitute Fwd and Rev Telomere Repeats in mate sequence
  dd$read2.2 <- gsub(pattern = 'TTAGGG', replacement = '__fw__',x = dd$read2)
  dd$read2.2 <- gsub(pattern = 'CCCTAA', replacement = '__rv__',x = dd$read2.2)
  
  # Get orientation of mate
  dd$orientation_mate <- apply(dd,1,function(x) if (grepl(pattern = '__fw__[ACTGN]*__rv__', x = x['read2.2'])){
    'end-to-head'
  } else if (grepl(pattern = '__rv__[ACTGN]*__fw__', x = x['read2.2'])) {
    'outward'
  } else {
    'no_orientation'
  })
  
  # Get number of break points in mate
  # For inward (end-to-head) and outward separately
  dd$Inward_breakpoint_count_mate <- str_count(string = dd$read2.2, pattern = '__fw__[ACTGN]*__rv__')
  dd$Outward_breakpoint_count_mate <- str_count(string = dd$read2.2, pattern = '__rv__[ACTGN]*__fw__')
  dd$Total_breakpoint_count_mate <- dd$Inward_breakpoint_count_mate + dd$Outward_breakpoint_count_mate
  
  # Search middle sequence in read 1 in mate
  # It is useful to check overlapping reads (it considers rev-comp)
  # middle_in_mate_count > 0 : Number of time that the middle was found in mate
  # middle_in_mate_count = 0 : Middle not found in mate
  # middle_in_mate_count = -1 : Middle in read one wrong
  # middle_in_mate_count < -1 : Number of time middle sequence is found in Read1
  dd$middle_in_mate_count <- apply(dd,1,function(x) search_middle_in_mate_one_mismatch(as.character(x['middle']),as.character(x['SEQ']),as.character(x['read2']),as.character(x['orientation'])))
  dd$middle_check <- apply(dd, 1, function(x) if (as.numeric(x['middle_in_mate_count']) < -1){
    'Multiple_middle_in_read1'
  } else if (as.numeric(x['middle_in_mate_count']) == -1) {
    'Middle_problem_in_read1'
  } else if (as.numeric(x['middle_in_mate_count']) == 0) {
    'No_middle_in_read2'
  } else if (as.numeric(x['middle_in_mate_count']) == 1) {
    'Middle_in_read2'
  } else if (as.numeric(x['middle_in_mate_count']) > 1) {
    'Multiple_middle_in_read2'
  } else {
    'Other'
  }
  )
  
  dd$Pair_comp <- apply(dd,1,function(x) if (x['middle_check'] %in% c('Other','Multiple_middle_in_read1','Multiple_middle_in_read2','Middle_problem_in_read1')){
    x['middle_check']
  } else if (x['Total_breakpoint_count_mate'] > 1){
    'Multiple_breakpoints_in_read2'
  } else if (x['Total_breakpoint_count_mate'] < 1) {
    'No_breakpoints_in_read2'
  } else if (x['Total_breakpoint_count_mate'] == 1) {
    if (x['middle_check'] == 'No_middle_in_read2' & x['orientation'] == x['orientation_mate']){
      'No_middle_same_TF_orient'
    } else if (x['middle_check'] == 'Middle_in_read2' & x['orientation'] == x['orientation_mate']){
      'Same_middle_same_TF_orient'
    } else if (x['orientation'] != x['orientation_mate']) {
      'In-out'
    } else {
      'Other'
    }
  } else {
    'Other'
  }
  )
  
  
  ##-----------------------
  # Analyse mate reads
  # Check number of breakpoints in mate
  # Check if middle sequence in read 1 is found in mate
  ##-----------------------
  
  # Substitute Fwd and Rev Telomere Repeats in mate sequence
  dd$read2.2 <- gsub(pattern = 'TTAGGG', replacement = '__fw__',x = dd$read2)
  dd$read2.2 <- gsub(pattern = 'CCCTAA', replacement = '__rv__',x = dd$read2.2)
  
  # Get orientation of mate
  dd$orientation_mate <- apply(dd,1,function(x) if (grepl(pattern = '__fw__[ACTGN]*__rv__', x = x['read2.2'])){
    'end-to-head'
  } else if (grepl(pattern = '__rv__[ACTGN]*__fw__', x = x['read2.2'])) {
    'outward'
  } else {
    'no_orientation'
  })
  
  # Get number of break points in mate
  # For inward (end-to-head) and outward separately
  dd$Inward_breakpoint_count_mate <- str_count(string = dd$read2.2, pattern = '__fw__[ACTGN]*__rv__')
  dd$Outward_breakpoint_count_mate <- str_count(string = dd$read2.2, pattern = '__rv__[ACTGN]*__fw__')
  dd$Total_breakpoint_count_mate <- dd$Inward_breakpoint_count_mate + dd$Outward_breakpoint_count_mate
  
  # Search middle sequence in read 1 in mate
  # It is useful to check overlapping reads (it considers rev-comp)
  # middle_in_mate_count > 0 : Number of time that the middle was found in mate
  # middle_in_mate_count = 0 : Middle not found in mate
  # middle_in_mate_count = -1 : Middle in read one wrong
  # middle_in_mate_count < -1 : Number of time middle sequence is found in Read1
  dd$middle_in_mate_count <- apply(dd,1,function(x) search_middle_in_mate_one_mismatch(as.character(x['middle']),as.character(x['SEQ']),as.character(x['read2']),as.character(x['orientation'])))
  dd$middle_check <- apply(dd, 1, function(x) if (as.numeric(x['middle_in_mate_count']) < -1){
    'Multiple_middle_in_read1'
  } else if (as.numeric(x['middle_in_mate_count']) == -1) {
    'Middle_problem_in_read1'
  } else if (as.numeric(x['middle_in_mate_count']) == 0) {
    'No_middle_in_read2'
  } else if (as.numeric(x['middle_in_mate_count']) == 1) {
    'Middle_in_read2'
  } else if (as.numeric(x['middle_in_mate_count']) > 1) {
    'Multiple_middle_in_read2'
  } else {
    'Other'
  }
  )
  
  dd$Pair_comp <- apply(dd,1,function(x) if (x['middle_check'] %in% c('Other','Multiple_middle_in_read1','Multiple_middle_in_read2','Middle_problem_in_read1')){
    x['middle_check']
  } else if (x['Total_breakpoint_count_mate'] > 1){
    'Multiple_breakpoints_in_read2'
  } else if (x['Total_breakpoint_count_mate'] < 1) {
    'No_breakpoints_in_read2'
  } else if (x['Total_breakpoint_count_mate'] == 1) {
    if (x['middle_check'] == 'No_middle_in_read2' & x['orientation'] == x['orientation_mate']){
      'No_middle_same_TF_orient'
    } else if (x['middle_check'] == 'Middle_in_read2' & x['orientation'] == x['orientation_mate']){
      'Same_middle_same_TF_orient'
    } else if (x['orientation'] != x['orientation_mate']) {
      'In-out'
    } else {
      'Other'
    }
  } else {
    'Other'
  }
  )
  
  
  ##-----------------------
  # Remove one of the mates when both show the same middle
  # Possibly caused because of small insert size
  ##-----------------------
  
  # Count only one of the mates if they overlap and represent the same dna fragment (same)
  fp_now <- NULL
  dd$READ_NAME_CODE <- with(dd, paste(sample,read_name,DUPLICATE,Supplementary,Secondary, sep = '-'))
  dd2 <- NULL
  
  # Let's check for overlapping reads
  # with two entries in the file
  REPEATS <- dd %>%
    dplyr::group_by(READ_NAME_CODE) %>%
    dplyr::summarise(N = length(orientation),
                     READ1 = length(unique(R1)),
                     INFO = paste(unique(Pair_comp), collapse = ','),
                     INFO_l = length(unique(Pair_comp)),
                     Error_like = sum(Pair_comp %in% c('Middle_problem_in_read1','Multiple_breakpoints_in_read2','Multiple_middle_in_read1','Multiple_middle_in_read2','No_middle_same_TF_orient')))
  
  
  # Take a decision for this filter
  REPEATS$Decision <- apply(REPEATS, 1, function(x) if (as.numeric(x['N']) < 2){
    'Correct1'
  } else if (as.numeric(x['N']) == 2){
    if (as.numeric(x['Error_like']) > 0){
      'Mate_problems'
    } else {
      'Correct2'
    }
  } else {
    'Other'
  }
  )
  
  # Split them in groups
  # Supported by one mate
  ONE_MATE <- REPEATS[REPEATS$Decision == 'Correct1',]
  one_mate <- dd[dd$READ_NAME_CODE %in% ONE_MATE$READ_NAME_CODE,]
  if (nrow(one_mate) > 0){
    one_mate$Paired_type <- 'Correct'
  }
  
  # Supported by both mates (only one must be counted)
  BOTH_MATES <- REPEATS[REPEATS$Decision == 'Correct2' & !grepl(pattern = 'In-out',x = REPEATS$INFO),]
  both_mate <- dd[dd$READ_NAME_CODE %in% BOTH_MATES$READ_NAME_CODE,]
  if (nrow(both_mate) > 0){
    # The one of both mates to be considered for downstream analysis
    both_mate_included <- both_mate %>%
      dplyr::group_by(READ_NAME_CODE) %>%
      dplyr::slice(1)
    both_mate_included$Paired_type <- 'Supported by both mates (overlapping)'
    
    # The one of both mates to be considered for downstream analysis
    both_mate_ignored <- both_mate %>%
      dplyr::group_by(READ_NAME_CODE) %>%
      dplyr::slice(n())
    both_mate_ignored$Paired_type <- 'Supported by both mates (fragment already counted)'
    
    # Merge results in one variable
    both_mate <- as.data.frame(rbind(both_mate_included,both_mate_ignored))
  }
  
  # Fusion found in both mates (in different orientation)
  # Cirular-like
  CIRCULAR_MATES <- REPEATS[REPEATS$Decision == 'Correct2' & grepl(pattern = 'In-out',x = REPEATS$INFO),]
  circular_mate <- dd[dd$READ_NAME_CODE %in% CIRCULAR_MATES$READ_NAME_CODE,]
  if (nrow(circular_mate) > 0){
    # Mark circular mates (1 or 2)
    circular_mate1 <- circular_mate %>%
      dplyr::group_by(READ_NAME_CODE) %>%
      dplyr::slice(1)
    circular_mate1$Paired_type <- 'Circular-like (mate 1)'
    circular_mate1$Subtype <- 'Circular-like (mate 1)'
    
    # Mark circular mates (1 or 2)
    circular_mate2 <- circular_mate %>%
      dplyr::group_by(READ_NAME_CODE) %>%
      dplyr::slice(n())
    circular_mate2$Paired_type <- 'Circular-like (mate 2)'
    circular_mate2$Subtype <- 'Circular-like (mate 2)'
    
    # Merge results in one variable
    circular_mate <- as.data.frame(rbind(circular_mate1,circular_mate2))
  }
  
  # Mate problems
  # Middle different between mates
  # More than two mates
  # Both read1 (or read2)
  MATE_PROBLEMS <- REPEATS[REPEATS$Decision %in% c('Mate_problems','Other'),]
  mate_problems <- dd[dd$READ_NAME_CODE %in% MATE_PROBLEMS$READ_NAME_CODE,]
  if (nrow(mate_problems) > 0){
    mate_problems$Paired_type <- 'Mate_problems'
  }
  
  # Merge all in one
  dd2 <- rbind(one_mate, both_mate,circular_mate,mate_problems)
  
  # Re-check In-out
  idx = which(dd2$Pair_comp == 'In-out' & dd2$Paired_type == 'Correct')
  in_out =  dd2[idx,]
  if (nrow(in_out) > 0){
    in_out$Paired_type <- 'Circular-like (mate 1)'
  }
  if(length(idx) > 0){
    dd2 = dd2[-idx,]
    dd2 <- rbind(dd2,in_out)
  } else {
    dd2 <- dd2
  }
  
  ##-----------------------
  # Change chromosome final 
  ##-----------------------
  # Trusting mapping of the mate first, and read with fusion the second
  if (nrow(dd2) > 0){
    dd2$chr_final <- apply(dd2, 1, function(x) if(x['chr_final'] %in% c('2_endogenous','4_endogenous','9_endogenous')){
      x['chr_final']
    } else {
      if (x['type_fus2'] %in% c('chromosomic mate')){
        x['chr2']
      } else if (x['type_fus1'] %in% c('chromosomic read')){
        x['chr']
      } else if (x['type_fus2'] %in% c('subtelomeric mate')){
        x['chr2']
      } else if (x['type_fus1'] %in% c('subtelomeric read')){
        x['chr']
      } else {
        'Unknown'
      }
    }
    )
  }
  
  # General subtype
  if (nrow(dd2) > 0){
    dd2$type_fus <- apply(dd2, 1, function(x) if(x['chr_final'] %in% c('2_endogenous','4_endogenous','9_endogenous')){
      x['chr_final']
    } else {
      if (x['type_fus2'] %in% c('chromosomic mate')){
        gsub(pattern = ' mate', replacement = '', x = x['type_fus2'])
      } else if (x['type_fus1'] %in% c('chromosomic read')){
        gsub(pattern = ' read', replacement = '', x = x['type_fus1'])
      } else if (x['type_fus2'] %in% c('subtelomeric mate')){
        gsub(pattern = ' mate', replacement = '', x = x['type_fus2'])
      } else if (x['type_fus1'] %in% c('subtelomeric read')){
        gsub(pattern = ' read', replacement = '', x = x['type_fus1'])
      } else {
        'fusion'
      }
    }
    )
  }
  
  ##-----------------------
  # Remove fusions where we find repeats forward and reverse in the mates. These are likely FP
  ##-----------------------
  dd <- dd2
  if(nrow(dd)>0){
    ##-----------------
    # Filter column
    ##-----------------
    if (nrow(dd) > 0){
      dd$Filter <- with(dd, paste(DUPLICATE,Supplementary,Secondary,Num_breakpoints,Subtype,proper_mate,type_fus,Pair_comp,Paired_type,type_fus1,type_fus2,mate_mapping, sep = ','))
      
      dd$Filter <- gsub(pattern = ",Supported by one mate$",replacement = '', x = dd$Filter)
      dd$Filter <- gsub(pattern = ",Supported by both mates \\(overlapping\\)$",replacement = '', x = dd$Filter)
      dd$Filter <- gsub(pattern = "Num_breakpoints_middle=1,",replacement = '', x = dd$Filter)
      dd$Filter <- gsub(pattern = "9_endogenous,",replacement = '', x = dd$Filter)
      dd$Filter <- gsub(pattern = "2_endogenous,",replacement = '', x = dd$Filter)
      dd$Filter <- gsub(pattern = "4_endogenous,",replacement = '', x = dd$Filter)
      dd$Filter <- gsub(pattern = "Other,",replacement = '', x = dd$Filter)
      dd$Filter <- gsub(pattern = "Circular-like \\(mate 1\\)",replacement = '', x = dd$Filter)
      dd$Filter <- gsub(pattern = "Circular-like \\(mate 2\\)",replacement = '', x = dd$Filter)
      dd$Filter <- gsub(pattern = "Fusion-like,",replacement = '', x = dd$Filter)
      dd$Filter <- gsub(pattern = "Inversion-like,",replacement = '', x = dd$Filter)
      dd$Filter <- gsub(pattern = "Correct,",replacement = '', x = dd$Filter)
      dd$Filter <- gsub(pattern = "In-out,",replacement = '', x = dd$Filter)
      dd$Filter <- gsub(pattern = "No_breakpoints_in_read2,",replacement = '', x = dd$Filter)
      dd$Filter <- gsub(pattern = "Correct_mate_found,",replacement = '', x = dd$Filter)
      dd$Filter <- gsub(pattern = "fusion",replacement = '', x = dd$Filter)
      dd$Filter <- gsub(pattern = "expected_mate_mapping",replacement = '', x = dd$Filter)
      dd$Filter <- gsub(pattern = "subtelomeric read",replacement = '', x = dd$Filter)
      dd$Filter <- gsub(pattern = "subtelomeric mate",replacement = '', x = dd$Filter)
      dd$Filter <- gsub(pattern = "subtelomeric",replacement = '', x = dd$Filter)
      dd$Filter <- gsub(pattern = ",+",replacement = ',', x = dd$Filter)
      dd$Filter <- gsub(pattern = "^,",replacement = '', x = dd$Filter)
      dd$Filter <- gsub(pattern = ",$",replacement = '', x = dd$Filter)
      
      
      dd$PASS <- apply(dd, 1, function(x) if (x['DUPLICATE'] == 'Correct' & 
                                              x['Supplementary'] == 'Correct' &
                                              x['Secondary'] == 'Correct' &
                                              x['Num_breakpoints'] == 'Num_breakpoints_middle=1' &
                                              x['proper_mate'] == 'Correct_mate_found' &
                                              x['type_fus'] %in% c('2_endogenous', '4_endogenous', '9_endogenous', 'fusion','subtelomeric') & 
                                              x['mate_mapping'] == "expected_mate_mapping" &
                                              !(x['Subtype'] == 'No_repeats_mate' & x['type_fus'] == 'fusion') & 
                                              !x['Pair_comp'] %in% c('Middle_problem_in_read1','Multiple_breakpoints_in_read2','Multiple_middle_in_read1','Multiple_middle_in_read2','No_middle_same_TF_orient') &
                                              x['Paired_type'] %in% c('Correct',
                                                                      'Supported by one mate',
                                                                      'Supported by both mates (overlapping)',
                                                                      'Circular-like (mate 1)',
                                                                      'Circular-like (mate 2)')) {
        "PASS"
      } else {
        'False_Positive'
      }
      )
      
      # Final print
      Read_lengths <- c(nchar(dd[!is.na(dd$SEQ),]$SEQ),nchar(dd[!is.na(dd$SEQ) & !is.na(dd$MAPQ_mate),]$SEQ))
      READ_length_mean <- round(mean(Read_lengths, na.rm = T),2)
      READ_length_approx <- max(nchar(dd$SEQ),nchar(dd$read2), na.rm = T)
      print (paste0('Number of samples with fusions: ', length(unique(dd$sample))))
      print (paste0('Number of original fusion calls: ', nrow(dd)))
      print (paste0('Number of pass fusion calls: ', sum(dd$PASS == 'PASS')))
      print (paste0('Number of potential false positives: ', sum(dd$PASS != 'PASS')))
      
      # Add read_length and calculate coverage
      if (is.null(read_length)){
        dd$read_length <- READ_length_approx
      } else {
        dd$read_length <- read_length
      }
      
      # Calculate coverage
      dd$Total_reads <- as.numeric(dd$Total_reads)
      dd$Supplementary_reads <- as.numeric(dd$Supplementary_reads)
      dd$Duplicate_reads <- as.numeric(dd$Duplicate_reads)
      dd$Total_reads_used <- dd$Total_reads - dd$Supplementary_reads - dd$Duplicate_reads
      
      dd$coverage <- dd$read_length*dd$Total_reads_used/3000000000
      
      # Filtering data.frame
      STATS_temp <- data.frame(project = project, N_samples_with_fusions = length(unique(dd$sample)), Mean_read_length = READ_length_mean, Max_read_length = READ_length_approx, raw_fusion_reads = nrow(dd),
                               pass_fusion_reads = sum(dd$PASS == 'PASS'), fp_fusion_reads = sum(dd$PASS != 'PASS'), markduplicates = CHECK)
    }
  } else {
    STATS_temp <- NULL
  }
  
  return(list(fusions = dd, stats = STATS_temp))
}

#-------
# Telomeric regions
#-------

# https://blog.gene-test.com/telomeric-regions-of-the-human-genome/
# Hg19 reference genome
telo_hg19 = data.frame(c("1","10","11","12","13","14","15","16","18","19","2","20","21","22","3","4","5","6","7","8","9","X","Y"),
                       c(249240621,135524747,134996516,133841895,115159878,107339540,102521392,90344753,78067248,59118983,243189373,63015520,48119895,51294566,198012430,191144276,180905260,171105067,159128663,146354022,141203431,155260560,59363566),
                       c(249250621,135534747,135006516,133851895,115169878,107349540,102531392,90354753,78077248,59128983,243199373,63025520,48129895,51304566,198022430,191154276,180915260,171115067,159138663,146364022,141213431,155270560,59373566),
                       c("telomere","telomere","telomere","telomere","telomere","telomere","telomere","telomere","telomere","telomere","telomere","telomere","telomere","telomere","telomere","telomere","telomere","telomere","telomere","telomere","telomere","telomere","telomere"))
colnames(telo_hg19) <- c("chr", "start", "end", "region")
telo_hg19$start = telo_hg19$start - 100000
# at the start we consider the telo region as 10Kb

# Hg38 reference genome
telo_hg38 <- data.frame(c("1","10","11","12","13","14","15","16","17","18","19","2","20","21","22","3","4","5","6","7","8","9","X","Y"),
                       c(248946422,133787422,135076622,133265309,114354328,107033718,101981189,90328345,83247441,80363285,58607616,242183529,64434167,46699983,50808468,198285559,190204555,181528259,170795979,159335973,145128636,138384717,156030895,57217415),
                       c(248956422,133797422,135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,242193529,64444167,46709983,50818468,198295559,190214555,181538259,170805979,159345973,145138636,138394717,156040895,57227415),
                       c("telomere","telomere","telomere","telomere","telomere","telomere","telomere","telomere","telomere","telomere","telomere","telomere","telomere","telomere","telomere","telomere","telomere","telomere","telomere","telomere","telomere","telomere","telomere","telomere"))
colnames(telo_hg38) = c("chr", "start", "end", "region")
telo_hg38$start = telo_hg38$start - 100000
# at the start we consider the telo region as 10Kb

#-------
# Header of the summaries
#-------
HEADER = c(
  "read_name","flag","chr",
  "pos",
  "MAPQ",
  "CIGAR",
  "RNEXT",
  "PNEXT",
  "TLEN",
  "SEQ",
  "forward0",
  "reverse0",
  "forward1",
  "reverse1",
  "forward2",
  "reverse2",
  "orientation",
  "left",
  "middle",
  "right",
  "flag2",
  "chr2",
  "pos2",
  "MAPQ_mate",
  "CIGAR_mate",
  "read2",
  "telo_repeats_mate",
  "forward_mate",
  "reverse_mate",
  "total_count",
  "sample",
  "Path_file",
  "File","Total_reads",
  "Supplementary_reads",
  "Duplicate_reads",
  "Paired_reads"
)


#---------------------------------------
#---------------------------------------
# Main computation
# Computing the QC functions
#---------------------------------------
#---------------------------------------

# 1. Arguments
parser <- ArgumentParser()

# setting parameters
parser$add_argument("-summary_file", "--summary_file", type="character", help="Summary file with fusions", metavar="file", required=T)
parser$add_argument("-ref_genome", "--ref_genome", choices=c('Hg19', 'Hg38'), type="character", help="Reference genome used for mapping", nargs=1, required=T)
parser$add_argument("-outprefix", "--outprefix", help="Prefix id of the out files. It can contain the folder path", nargs=1, required=T)
parser$add_argument("-read_length", "--read_length", type="integer",help="Read length used for sequencing. If not provided, it will be estimated from summary file", nargs=1, required=F)

# 1. Reading parameters
args <- parser$parse_args()

file <- args$summary_file
REFERENCE <- args$ref_genome
prefix <- args$outprefix
read_length <- args$read_length

# Get project ID - obtained from prefix
project <- basename(prefix)

# 2. Running main function
if (file.exists(file)){
  # Run tool
  RESULTS <- MAIN_FUNCTION_FUSION_QC(file, project,REFERENCE,read_length)
  fusions <- RESULTS$fusions
  stats <- RESULTS$stats
  
  } else {
  print (c(file, ' does not exist'))
}

# 3. Creating output files and folder
outf <- dirname(prefix)
dir.create(outf, showWarnings = FALSE)

if (!is.null(fusions)){
  # 3.1: All read-fusions
  out1 <- paste(prefix,'.fusions.unfiltered.tsv',sep = '')
  write.table(x = fusions, file = out1, sep = '\t', quote = F, row.names = F, col.names = T)
  
  # 3.2: All pass reads
  out2 <- paste(prefix,'.fusions.pass.tsv',sep = '')
  PASS <- fusions[fusions$PASS == 'PASS',]
  write.table(x = PASS, file = out2, sep = '\t', quote = F, row.names = F, col.names = T)
  
  # 3.3: False positive reads
  out3 <- paste(prefix,'.fusions.false_positives.tsv',sep = '')
  FP <- fusions[fusions$PASS != 'PASS',]
  write.table(x = FP, file = out3, sep = '\t', quote = F, row.names = F, col.names = T)
  
  # 3.5: Stats per sample
  out5 <- paste(prefix,'.fusions.sample_stats.tsv',sep = '')
  fusions$PASS <- (fusions$PASS == 'PASS')*1
  fusions$Total <- 1
  COUNT <- aggregate(cbind(PASS,Total) ~ sample + project,data = fusions,sum)
  COUNT$Prop_pass <- COUNT$PASS / COUNT$Total
  colnames(COUNT) <- c('sample','project','pass_fusion_reads','total_fusion_reads','Prop_pass')
    
  write.table(x = COUNT, file = out5, sep = '\t', quote = F, row.names = F, col.names = T)
  
  # 3.6: Collapsed events
  out6 <- paste(prefix,'.fusions.pass.collapsed.tsv',sep = '')
  all = ddply(PASS, .(chr_final,label,sample,project,Subtype,type_fus,
                      orientation,middle,Total_reads,Supplementary_reads,
                      Duplicate_reads,Paired_reads,Total_reads_used,
                      read_length, coverage), summarise, n=length(SEQ))
  colnames(all)[1] <- 'chr'
  write.table(x = all, file = out6, sep = '\t', quote = F, row.names = F, col.names = T)
  
  # 3.7: Save image
  out7 <- paste(prefix,'.fusions.QC.Rdata',sep = '')
  save.image(out7)

} else {
  print (paste0('No results obtained. Check if input file is empty'))
}









