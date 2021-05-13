#---------------------
# Script to perform the QC in the summary file
#---------------------

#-------
# Packages required
#-------
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(argparse))
options(warn=-1)
#-------
# Function to run the QC fusion detection step
#-------

# Function to detect the number of breakpoints in the middle sequence
NUM_BREAKPOINTS <- function(SEQ, orientation){
  
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
  idx = which(DATA$chr == CHROM & 
                DATA$pos >= MIN_COORDINATE & DATA$pos <= MAX_COORDINATE & 
                DATA$MAPQ > 0 & DATA$orientation %in% orientation)
  DATA$chr[idx] = Flag_endo
  
  # Flag the endogenous sequence in read mate
  DATA$chr2 =gsub("chr","",as.vector(DATA$chr2))
  DATA$chr2 <- as.character(DATA$chr2)
  idx = which(DATA$chr2 == CHROM & 
                DATA$PNEXT >= MIN_COORDINATE & DATA$PNEXT <= MAX_COORDINATE & 
                DATA$MAPQ_mate > 0 & DATA$orientation %in% orientation)
  
  DATA$chr2[idx] = Flag_endo
  
  return(DATA)
}

# Main function to run Quality Control (QC) step with summary files 
MAIN_FUNCTION_FUSION_QC <- function(file, project,REFERENCE){
  # Initial print
  cat("\n")
  print("___")
  print(paste0('Project: ', project))
  
  
  # Load project file and add header
  d0 = as.data.frame(fread(cmd = paste0("awk -F'\t' '{if (NF == 29) print $0}' ", file), sep = '\t', header = F))
  if (nrow(d0) > 0){
    d0 <- add_column(d0, extra2 = 0, .after = "V27")
    d0 <- add_column(d0, extra1 = 0, .after = "V27")
    names(d0) <- HEADER
  }
  
  # Load project file and add header
  d1 = as.data.frame(fread(cmd = paste0("awk -F'\t' '{if (NF == 31) print $0}' ", file), sep = '\t', header = F))
  if (nrow(d1) > 0){
    names(d1) <- HEADER
  }
  
  # Remove extra columns
  d2 = as.data.frame(fread(cmd = paste0("awk -F'\t' '{if (NF == 33) print $0}' ", file), sep = '\t', header = F))
  if (nrow(d2 > 0)){
    d2$V28 <- NULL
    d2$V29 <- NULL
    colnames(d2) <- colnames(d1)
  }
  
  # Check if empty files
  d <- rbind(d0,d1, d2)
  if (nrow(d) < 1){
    print(paste0('Warning: ',project, ' files are empty or not readable for the tool'))
    return(list(fusions = NULL, stats = NULL))
  }
  names(d) <- HEADER
  
  # Clean sample name
  d$sample = gsub( "_rawFus_full_chr[0-9]_python.fq", "", d$sample)
  d$sample = gsub( "_rawFus_full_chrX_python.fq", "", d$sample)
  d$sample = gsub( "_rawFus_full_chrY_python.fq", "", d$sample)
  d$sample = gsub( "_rawFus_full_chr[(0-9)][(0-9)]_python.fq", "", d$sample)
  d$sample = gsub( "_rawFus_unmapped_python.fq", "", d$sample)
  d$sample = gsub( ".sorted", "", d$sample)
  
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
    dd$proper_mate <- !is.na(dd$pos2) & dd$pos2 == dd$PNEXT 
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
  dd$Bases_mapped_read_from_CIGAR <- cigarWidthAlongQuerySpace(dd$CIGAR)
  dd$Bases_mapped_mate_from_CIGAR <- cigarWidthAlongQuerySpace(dd$CIGAR_mate)
  dd$Prop_read_mapped <- dd$Bases_mapped_read_from_CIGAR/nchar(dd$SEQ)
  dd$Prop_mate_mapped <- dd$Bases_mapped_mate_from_CIGAR/nchar(dd$read2)
  
  # Get potential false positives
  # Not in endogenous regions, mapped perfectly in non-telomeric regions (middle of chromosomes)
  # We check afterwards if mates map to sub-telomeric regions
  idx = which(dd$telo_repeats_mate == "No" & dd$MAPQ_mate >= 20 & nchar(as.vector(dd$read2)) >= 80 & !dd$chr_final %in% c('2_endogenous','4_endogenous','9_endogenous') & dd$Prop_mate_mapped > 0.8)
  
  # Type of fusions
  dd$type_fus = "fusion"
  
  if(length(idx) > 0){
    
    # Mark them as subtelomeric 
    fp_now =  dd[idx,]
    fp_now$type_fus="subtelomeric"
    dd = dd[-idx,]
    dd$type_fus="fusion"
    
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
        fp_now[ff,]$type_fus <- 'chromosomic mate'
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
  idx = which(dd$telo_repeats_mate == "No" & dd$MAPQ_mate < 20 & nchar(as.vector(dd$read2)) >= 80 & !dd$chr_final %in% c('2_endogenous','4_endogenous','9_endogenous') & dd$type_fus == 'fusion')
  
  # Type of fusions
  dd$mate_mapping = "expected_mate_mapping"
  
  if (length(idx) > 0){
    fp_now =  dd[idx,]
    fp_now$mate_mapping="inconsistent_mate_mapping"
    dd = dd[-idx,]
    dd <- rbind(dd, fp_now)
  }
  
  ##-----------------------
  # Remove one of the mates when both show the same middle
  # Possibly caused because of small insert size
  ##-----------------------
  
  # Count only one of the mates if they overlap and represent the same dna fragment (same)
  fp_now <- NULL
  dd$READ_NAME_CODE <- with(dd, paste(sample,read_name,DUPLICATE,Supplementary,Secondary, sep = '-'))
  dd2 <- NULL
  
  TABLE <- as.data.frame(table(dd$READ_NAME_CODE))
  TABLE1 <- TABLE[TABLE$Freq < 2,]
  TABLE2 <- TABLE[TABLE$Freq > 1,]
  
  dd2 <- dd[dd$READ_NAME_CODE %in% TABLE1$Var1,]
  dd2$Paired_type <- 'Supported by one mate'
  
  dd_temp <- dd[dd$READ_NAME_CODE %in% TABLE2$Var1,]
  
  for (read_code in unique(TABLE2$Var1)){
    d_temp <- dd_temp[dd_temp$READ_NAME_CODE == read_code,]
    
    # If more than one mate counted
    if (nrow(d_temp) > 1){
      MIDDLE <- unique(d_temp$middle)
      
      if (length(MIDDLE) > 1){
        d_temp$Paired_type <- 'Mate problems'
        dd2 <- rbind(dd2,d_temp)
      } else {
        d_temp$Paired_type <- 'Supported by both mates (overlapping)'
        
        # We only keep one of them
        d_temp1 <- d_temp[1,]
        dd2 <- rbind(dd2,d_temp1)
        
        # Move the rest as FPs
        d_temp2 <- d_temp[seq(2,nrow(d_temp)),]
        d_temp2$Paired_type <- 'Supported by both mates (fragment already counted)'
        dd2 <- rbind(dd2,d_temp2)
      }
    } else {
      d_temp$Paired_type <- 'Supported by one mate'
      dd2 <- rbind(dd2,d_temp)
    }
  }
  
  ##-----------------------
  # Change chromosome final 
  ##-----------------------
  if (nrow(dd2) > 0){
    dd2$chr_final <- apply(dd2, 1, function(x) if(!x['chr_final'] %in% c('2_endogenous','4_endogenous','9_endogenous') & x['type_fus'] %in% c('chromosomic mate','subtelomeric')) {
      x['chr2']
    } else if (x['chr_final'] %in% c('2_endogenous','4_endogenous','9_endogenous')){
      x['chr_final']
    } else {
      'Unknown'
    }
    )
  }
  
  ##-----------------------
  # Remove fusions where we find repeats forward and reverse in the mates. These are likely FP
  ##-----------------------
  dd <- dd2
  if(nrow(dd)>0){
    #print("Check mates with repeats forward and reverse. Ignoring overlapping mates showing the same middle")
    idx = which(!dd$chr_final %in% c('2_endogenous','4_endogenous','9_endogenous') & 
                  dd$Paired_type == "Supported by one mate" & 
                  dd$forward_mate > 0 & dd$reverse_mate > 0 &
                  !(dd$RNEXT == '=' & abs(dd$pos - dd$PNEXT) < 75))
    
    # Update false positives
    fp_now =  dd[idx,]
    if (nrow(fp_now) > 0){
      fp_now$Paired_type <- 'Mate problems'
    }
    
    # Keep good ones
    if(length(idx) > 0){
      dd = dd[-idx,]
      dd <- rbind(dd,fp_now)
    } else {
      dd <- dd
    }
    
    ##-----------------
    # Filter column
    ##-----------------
    if (nrow(dd) > 0){
      dd$Filter <- with(dd, paste(DUPLICATE,Supplementary,Secondary,Num_breakpoints,Subtype,proper_mate,type_fus,mate_mapping,Paired_type, sep = ','))
      
      dd$Filter <- gsub(pattern = ",Supported by one mate$",replacement = '', x = dd$Filter)
      dd$Filter <- gsub(pattern = ",Supported by both mates \\(overlapping\\)$",replacement = '', x = dd$Filter)
      dd$Filter <- gsub(pattern = "Num_breakpoints_middle=1,",replacement = '', x = dd$Filter)
      dd$Filter <- gsub(pattern = "9_endogenous,",replacement = '', x = dd$Filter)
      dd$Filter <- gsub(pattern = "2_endogenous,",replacement = '', x = dd$Filter)
      dd$Filter <- gsub(pattern = "4_endogenous,",replacement = '', x = dd$Filter)
      dd$Filter <- gsub(pattern = "Fusion-like,",replacement = '', x = dd$Filter)
      dd$Filter <- gsub(pattern = "Inversion-like,",replacement = '', x = dd$Filter)
      dd$Filter <- gsub(pattern = "Correct,",replacement = '', x = dd$Filter)
      dd$Filter <- gsub(pattern = "Correct_mate_found,",replacement = '', x = dd$Filter)
      dd$Filter <- gsub(pattern = "fusion",replacement = '', x = dd$Filter)
      dd$Filter <- gsub(pattern = "expected_mate_mapping",replacement = '', x = dd$Filter)
      dd$Filter <- gsub(pattern = "subtelomeric",replacement = '', x = dd$Filter)
      dd$Filter <- gsub(pattern = ",+",replacement = ',', x = dd$Filter)
      dd$Filter <- gsub(pattern = "^,",replacement = '', x = dd$Filter)
      dd$Filter <- gsub(pattern = ",$",replacement = '', x = dd$Filter)
      
      
      dd$PASS <- apply(dd, 1, function(x) if (x['DUPLICATE'] == 'Correct' & 
                                              x['Supplementary'] == 'Correct' &
                                              x['Secondary'] == 'Correct' &
                                              x['Num_breakpoints'] == 'Num_breakpoints_middle=1' &
                                              x['proper_mate'] == 'Correct_mate_found' &
                                              x['type_fus'] %in% c('fusion','subtelomeric') & 
                                              x['mate_mapping'] == "expected_mate_mapping" &
                                              x['Subtype'] != 'Other' &
                                              x['Paired_type'] %in% c('Supported by one mate','Supported by both mates (overlapping)')) {
        "PASS"
      } else {
        'False_Positive'
      }
      )
      
      # Final print
      READ_length_mean <- round(mean(nchar(dd$SEQ), na.rm = T),2)
      print (paste0('Number of samples with fusions: ', length(unique(dd$sample))))
      print (paste0('Number of original fusion calls: ', nrow(dd)))
      print (paste0('Number of pass fusion calls: ', sum(dd$PASS == 'PASS')))
      print (paste0('Number of potential false positives: ', sum(dd$PASS != 'PASS')))
      
      # Filtering data.frame
      STATS_temp <- data.frame(project = project, N_samples_with_fusions = length(unique(dd$sample)), Mean_read_length = READ_length_mean, raw_fusion_reads = nrow(dd),
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
  "sample"
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
parser$add_argument("-project", "--project", default='Project',help="Project Id to be used", nargs=1, required=F)
parser$add_argument("-prefix", "--prefix", default = NULL, help="Prefix id of the out files. It can contain the folder path", nargs=1, required=F)

# 1. Reading parameters
args <- parser$parse_args()

file <- args$summary_file
REFERENCE <- args$ref_genome
project <- args$project
prefix <- args$prefix

if (is.null(prefix)){
  prefix <- paste('.',project, sep = '/')
}

# 2. Running main function
if (file.exists(file)){
  # Run tool
  RESULTS <- MAIN_FUNCTION_FUSION_QC(file, project,REFERENCE)
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
  
  # 3.4: Write stats
  out4 <- paste(prefix,'.fusions.project_stats.tsv',sep = '')
  stats$Prop_pass <- stats$pass_fusion_reads/stats$raw_fusion_reads
  write.table(x = stats, file = out4, sep = '\t', quote = F, row.names = F, col.names = T)
  
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
  all = ddply(PASS, .(chr_final,label,sample,project,Subtype,type_fus,orientation,middle), summarise, n=length(SEQ))
  colnames(all)[1] <- 'chr'
  write.table(x = all, file = out6, sep = '\t', quote = F, row.names = F, col.names = T)
  
  # 3.7: Save image
  out7 <- paste(prefix,'.fusions.QC.Rdata',sep = '')
  save.image(out7)

} else {
  print (paste0('No results obtained. Check if input file is empty'))
}









