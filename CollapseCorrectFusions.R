#---------------------
# Script to perform the fusion error correction and final collapse
#---------------------

#-------
# Packages required
#-------
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(argparse))
options(warn=-1)

#-------
# Functions
#-------
# Creates pure middle sequence table
CreatePureMiddleTable <- function(){
  # Forward repeat
  fwd_template <- 'TTAGGG'
  fwd_split <- strsplit(fwd_template, '', fixed=FALSE)[[1]]
  
  # Reverse repeat
  rev_template <- 'CCCTAA'
  rev_split <- strsplit(rev_template, '', fixed=FALSE)[[1]]
  
  # Creating end-to-head (also known as inward) possible middle sequences
  INWARD_MIDDLE <- NULL
  for (i in seq(0,5)){
    # Forward
    seq_fwd <- seq(0, i)
    TEMP_fwd <- fwd_split[seq_fwd]
    A_fwd <- paste(TEMP_fwd, collapse = '')
    
    for (i2 in seq(1,6)){
      # Reverse
      if (i2 == 1){
        A_rev <- ''
      } else {
        seq_rev <- seq(i2, 6)
        TEMP_rev <- rev_split[seq_rev]
        A_rev <- paste(TEMP_rev, collapse = '')
      }
      
      PASTED <- paste(A_fwd, A_rev, sep = '')
      DF_temp <- data.frame('end-to-head', PASTED, A_fwd, A_rev, i, i2-1)
      colnames(DF_temp) <- c('Type','Sequence', 'Seq_fwd','Seq_rev','Cut_in_fwd','Cut_in_rev')
      
      INWARD_MIDDLE <- rbind(INWARD_MIDDLE, DF_temp)
    }
  }
  INWARD_MIDDLE <- unique(INWARD_MIDDLE)
  
  # Creating outward possible middle sequences
  fwd_template <- 'TTAGGG'
  fwd_split <- strsplit(fwd_template, '', fixed=FALSE)[[1]]
  
  rev_template <- 'CCCTAA'
  rev_split <- strsplit(rev_template, '', fixed=FALSE)[[1]]
  
  OUTWARD_MIDDLE <- NULL
  for (i in seq(0,5)){
    
    # Reverse
    seq_rev <- seq(0, i)
    TEMP_rev <- rev_split[seq_rev]
    A_rev <- paste(TEMP_rev, collapse = '')
    
    for (i2 in seq(1,6)){
      # Forward
      if (i2 == 1){
        A_fwd <- ''
      } else {
        seq_fwd <- seq(i2, 6)
        TEMP_fwd <- fwd_split[seq_fwd]
        A_fwd <- paste(TEMP_fwd, collapse = '')
      }
      
      PASTED <- paste(A_rev, A_fwd, sep = '')
      DF_temp2 <- data.frame('outward', PASTED, A_rev, A_fwd, i, i2-1)
      colnames(DF_temp2) <- c('Type','Sequence', 'Seq_rev','Seq_fwd','Cut_in_rev','Cut_in_fwd')
      
      OUTWARD_MIDDLE <- rbind(OUTWARD_MIDDLE, DF_temp2)
    }
  }
  
  # Final table with all possible pure middle sequences
  PURE_MIDDLE <- rbind(INWARD_MIDDLE, OUTWARD_MIDDLE)
  PURE_MIDDLE$Sequence <- as.character(PURE_MIDDLE$Sequence)
  PURE_MIDDLE$LEN <- nchar(PURE_MIDDLE$Sequence)
  
  return(PURE_MIDDLE)
}

# Middle Correction starting from left
MiddleCorrectLeft <- function(SEQ){
  fwd_seq <- 'TTAGGG'
  rev_seq <- 'CCCTAA'
  
  SEQ <- gsub(pattern = '_reverse_', replacement = rev_seq, x = SEQ)
  SEQ <- gsub(pattern = '_forward_', replacement = fwd_seq, x = SEQ)
  
  SEQ2 <- unlist(strsplit(SEQ, '', fixed=FALSE))
  N <- nchar(SEQ)
  
  i_max <- N - 5
  new_seq <- ''
  i <- 1
  while (i <= N){
    start <- i
    BASE <- SEQ2[start]
    end <- start + 5
    end <- min(N,end)
    kmer_list <- SEQ2[start:end]
    kmer <- paste(kmer_list, collapse = '')
    
    FWD_dist <- adist(kmer,fwd_seq)[[1]]
    REV_dist <- adist(kmer,rev_seq)[[1]]
    
    if (FWD_dist < 2 & nchar(kmer) == 6){
      new_seq <- paste(new_seq,'__fw__',sep = '')
      i <- i + 6
    } else if(REV_dist < 2 & nchar(kmer) == 6){
      new_seq <- paste(new_seq,'__rv__',sep = '')
      i <- i + 6
    } else {
      new_seq <- paste(new_seq, BASE ,sep = '')
      i <- i + 1
    }
  }
  
  return(new_seq)
}

# Middle Correction starting from right
MiddleCorrectRight <- function(SEQ){
  fwd_seq <- 'TTAGGG'
  rev_seq <- 'CCCTAA'
  
  SEQ <- gsub(pattern = '_reverse_', replacement = rev_seq, x = SEQ)
  SEQ <- gsub(pattern = '_forward_', replacement = fwd_seq, x = SEQ)
  
  SEQ2 <- unlist(strsplit(SEQ, '', fixed=FALSE))
  N <- nchar(SEQ)
  
  i_min <- 5
  new_seq <- ''
  i <- N
  while (i > 0){
    end <- i
    start <- end - 5
    start <- max(1,start)
    BASE <- SEQ2[end]
    
    kmer_list <- SEQ2[start:end]
    kmer <- paste(kmer_list, collapse = '')
    
    FWD_dist <- adist(kmer,fwd_seq)[[1]]
    REV_dist <- adist(kmer,rev_seq)[[1]]
    
    if (FWD_dist < 2 & nchar(kmer) == 6){
      new_seq <- paste('__fw__',new_seq, sep = '')
      i <- i - 6
    } else if(REV_dist < 2 & nchar(kmer) == 6){
      new_seq <- paste('__rv__',new_seq,sep = '')
      i <- i - 6
    } else {
      new_seq <- paste(BASE,new_seq,sep = '')
      i <- i - 1
    }
    #print (c(start, end,kmer,new_seq))
  }
  
  return(new_seq)
}

# Clean corrected middle sequence
CleanMiddle <- function(SEQ){
  
  SEQ <- gsub(pattern = '_reverse_', replacement = '__rv__', x = SEQ)
  SEQ <- gsub(pattern = '_forward_', replacement = '__fw__', x = SEQ)
  
  SEQ2 <- gsub(pattern = "^([A|C|T|G|N]?[__fw__|__rv__]+)+", replacement = '', x = SEQ, perl = F)
  SEQ3 <- gsub(pattern = "([__fw__|__rv__]+[A|C|T|G|N]?)+$", replacement = '', x = SEQ2, perl = F)
  return(SEQ3)
}

# Correct endo/pure sequences sequences allowing 1 mismatch or 1 indel
MiddleCorrectionEndoPure <- function(SEQ,orientation,chr,PURE_MIDDLE){
  if (chr != '9_endogenous'){
    # Step 1. Separate the middle sequences in Pure/Classical and alternative
    PURE2 <- PURE_MIDDLE[PURE_MIDDLE$Type == orientation & PURE_MIDDLE$LEN == nchar(SEQ),]
    seqs_exact <- PURE_MIDDLE[PURE_MIDDLE$Type == orientation & PURE_MIDDLE$LEN == nchar(SEQ),]$Sequence
    seqs_all <- PURE_MIDDLE[PURE_MIDDLE$Type == orientation,]$Sequence
    
    # Step 2. Check if the alternative is 1 mismatch distance to the expected pure sequences
    if (nchar(SEQ) >= 4 & nchar(SEQ) <= 10 & !SEQ %in% seqs_all){
      DISTANCES <- seqs_exact[adist(SEQ,seqs_exact) < 2]
      if (length(DISTANCES) < 1){
        DISTANCES <- seqs_all[adist(SEQ,seqs_all) < 2]
      }
    } else{
      DISTANCES <- character(0)
    }
    
    if (length(DISTANCES) > 0){
      new_seq <- DISTANCES[1]
    } else {
      new_seq <- SEQ
    }
    
  } else {
    endo_seq <- 'TTAA'
    
    N <- nchar(SEQ)
    
    if (chr == '9_endogenous'){
      dist <- adist(SEQ,endo_seq)[[1]]
    } else {
      dist <- 10
    }
    
    if (dist < 2){
      new_seq <- endo_seq
    } else {
      new_seq <- SEQ
    }
  }
  return (new_seq)
}

# Function to correct islands of sequences in the resulting middle (2 mismatches allowed)
IslandCorrection <- function(SEQ){
  fwd_seq <- 'TTAGGG'
  rev_seq <- 'CCCTAA'
  
  seqs <- c(fwd_seq,rev_seq)
  seqs_labels <- c('__fw__','__rv__')
  
  if (grepl(pattern = '[__fw__|__rv__]', SEQ)){
    SEQ2 <- unlist(strsplit(SEQ, split = '__'))
    SEQ2 <- SEQ2[SEQ2 != '']
    SEQ2 <- gsub(pattern = 'rv', replacement = '__rv__', x = SEQ2)
    SEQ2 <- gsub(pattern = 'fw', replacement = '__fw__', x = SEQ2)
    
    
    new_seq <- ''
    for (i in SEQ2){
      # One mismatch or 1 indel
      DIST <- seqs_labels[adist(i,seqs) < 2]
      if (length(DIST) > 0){
        seq_temp <- DIST[1]
      } else{
        # Two mismatches or 1 indel + 1 missmatch
        DIST <- seqs_labels[adist(i,seqs) < 3]
        if (length(DIST) > 0 & nchar(i) %in% seq(5,7)){
          seq_temp <- DIST[1]
        } else {
          seq_temp <- i
        }
      }
      new_seq <- paste(new_seq,seq_temp,sep = '')
    }
    return(new_seq)
  } else {
    return(SEQ)
  }
}

# Function to simplify middle sequences when we have different breakpoints
SimplifyMiddle <- function(SEQ, orientation){
  
  SEQ <- gsub(pattern = '_forward_', replacement = '__fw__',x = SEQ)
  SEQ <- gsub(pattern = '_reverse_', replacement = '__rv__',x = SEQ)
  
  # Split different groups between fw and rv groups
  ALL_ISLANDS <- unlist(strsplit(SEQ, split = '__'))
  ALL_ISLANDS <- ALL_ISLANDS[!ALL_ISLANDS %in%  c('','fw','rv')]

  # For each orientation
  if (orientation == 'outward'){
    SEQ2 <- paste('__rv__',SEQ,'__fw__',sep = '')
    ISLAND <- unlist(str_extract_all(string = SEQ2, pattern = "__rv__[A-Z]*__fw__"))
  } else {
    SEQ2 <- paste('__fw__',SEQ,'__rv__',sep = '')
    ISLAND <- unlist(str_extract_all(string = SEQ2, pattern = "__fw__[A-Z]*__rv__"))
  }
  
  # Check the number of different islands
  if (length(ISLAND) > 1) {
    return(SEQ)
  } else {
    ISLAND <- gsub(pattern = "[__fw__|__rv__]", replacement = '', x = ISLAND)
    OTHER_ISLANDS <- ALL_ISLANDS[ALL_ISLANDS != ISLAND]
    
    if (length (OTHER_ISLANDS) > 0){
      max_len <- max(nchar(OTHER_ISLANDS))
    } else {
      max_len <- 0
    }
    if (max_len <= 250){
      return(ISLAND)
    } else {
      return(SEQ)
    }
  }
}

# Calculate the number of different possible solutions after correction
PossibleChoicesPureCorrection <- function(SEQ,orientation,chr,PURE_MIDDLE){
  if (!chr %in% c('2_endogenous','4_endogenous','9_endogenous')){
    # Step 1. Separate the middle sequences in Pure/Classical and alternative
    PURE2 <- PURE_MIDDLE[PURE_MIDDLE$Type == orientation & PURE_MIDDLE$LEN == nchar(SEQ),]
    seqs_exact <- unique(PURE_MIDDLE[PURE_MIDDLE$Type == orientation & PURE_MIDDLE$LEN == nchar(SEQ),]$Sequence)
    seqs_all <- unique(PURE_MIDDLE[PURE_MIDDLE$Type == orientation,]$Sequence)
    
    # Step 2. Check if the alternative is 1 mismatch distance to the expected pure sequences
    if (nchar(SEQ) >= 4 & nchar(SEQ) <= 10 & !SEQ %in% seqs_all){
      DISTANCES <- seqs_exact[adist(SEQ,seqs_exact) < 2]
      if (length(DISTANCES) < 1){
        DISTANCES <- seqs_all[adist(SEQ,seqs_all) < 2]
      }
    } else if (SEQ %in% seqs_all){
      new_seq <- 1
      return(new_seq)
    } else{
      DISTANCES <- character(0)
    }
    
    # Check the number of possibilities
    if (length(DISTANCES) > 0){
      new_seq <- length(DISTANCES)
    } else {
      new_seq <- NA
    }
    
  } else {
    new_seq <- 1
  }
  return (new_seq)
}

# Re-check pure sequences (final step)
CorrectFinalPure <- function(SEQ, orientation,middle_0,middle_left,A){
  
  # Find the possible matches for the specific pure
  MIDDLE_LEFT <- A[A$Type == orientation & (A$middle_left_clean == SEQ | A$middle_right_clean == SEQ),]
  MIDDLE_LEFT <- MIDDLE_LEFT[!MIDDLE_LEFT$Sequence %in% c('CCCTAAGGG','CCCTTAGGG'),]
  
  final_seq <- NULL
  DF <- NULL
  if (nrow(MIDDLE_LEFT) > 1 & length(unique(MIDDLE_LEFT$Sequence)) > 1 & SEQ != middle_left){
    
    if (orientation == 'outward'){
      SEQ2 <- paste('__rv__',middle_left,'__fw__',sep = '')
      ISLAND <- unlist(str_extract_all(string = SEQ2, pattern = "__rv__[A-Z]*__fw__"))
    } else {
      SEQ2 <- paste('__fw__',middle_left,'__rv__',sep = '')
      ISLAND <- unlist(str_extract_all(string = SEQ2, pattern = "__fw__[A-Z]*__rv__"))
    }
    
    if (middle_left == SEQ){
      final_seq <- paste(SEQ,1,sep = ";")
    }
    
    middle_left0 <- middle_left
    middle0_copy <- middle_0
    if (orientation == 'outward'){
      if (!grepl(pattern = '^__rv__', x = middle_left0)){
        middle_left <- paste0('__rv__',middle_left)
        middle_0 <- paste0('CCCTAA',middle_0)
      }
      
      if (!grepl(pattern = '__fw__$', x = middle_left0)){
        middle_left <- paste0(middle_left,'__fw__')
        middle_0 <- paste0(middle_0,'TTAGGG')
      }
    } else {
      if (!grepl(pattern = '^__fw__', x = middle_left0)){
        middle_left <- paste0('__fw__',middle_left)
        middle_0 <- paste0('TTAGGG',middle_0)
      }
      
      if (!grepl(pattern = '__rv__$', x = middle_left0)){
        middle_left <- paste0(middle_left,'__rv__')
        middle_0 <- paste0(middle_0,'CCCTAA')
      }
    }
    
    
    if (length(ISLAND) < 2){
      
      # Get original sequence
      POS <- gregexpr(pattern = ISLAND, middle_left)[[1]][1]
      LEN <- nchar(ISLAND)
      ORIGINAL_SEQ <- substr(middle_0, POS, POS+LEN-1)
      
      BEST_DIST <- NULL
      BEST_SEQ <- NULL
      
      for (ROW in seq(1,nrow(MIDDLE_LEFT))){
        TEMP <- MIDDLE_LEFT[ROW,]
        ORIGINAL_PURE <- TEMP$Sequence
        LEFT <- TEMP$middle_left
        
        ORIGINAL_PURE2 <- ORIGINAL_PURE
        dist_subtract <- 0
        if (orientation == 'outward'){
          if (!grepl(pattern = '^__rv__', x = LEFT)){
            ORIGINAL_PURE2 <- paste0('CCCTAA',ORIGINAL_PURE2)
          }
          
          if (!grepl(pattern = '__fw__$', x = LEFT)){
            ORIGINAL_PURE2 <- paste0(ORIGINAL_PURE2,'TTAGGG')
          }
        } else {
          if (!grepl(pattern = '^__fw__', x = LEFT)){
            ORIGINAL_PURE2 <- paste0('TTAGGG',ORIGINAL_PURE2)
          }
          
          if (!grepl(pattern = '__rv__$', x = LEFT)){
            ORIGINAL_PURE2 <- paste0(ORIGINAL_PURE2,'CCCTAA')
          }
        }
        
        DIST <- adist(ORIGINAL_SEQ,ORIGINAL_PURE2)
        PENALTY_GAP <- abs(nchar(ORIGINAL_SEQ)-nchar(ORIGINAL_PURE2))
        DIST <- DIST + PENALTY_GAP
        
        DF_temp <- data.frame(SEQ,orientation,ORIGINAL_SEQ, middle_left, middle_left0,ORIGINAL_PURE2, ORIGINAL_PURE, LEFT, DIST)
        
        DF <- rbind(DF,DF_temp)
        
      }
      MIN <- min(DF$DIST)
      DF <- DF[DF$DIST <= min(MIN,2),]
      
      if (nrow(DF) > 0){
        PURES <- unique(DF$ORIGINAL_PURE)
        N <- length(PURES)
        if (N > 1){
          if (SEQ %in% PURES){
            final_seq <- paste(SEQ,N,sep = ";")
          } else{
            final_seq <- paste(PURES[1],N,sep = ";")
          }
        } else{
          final_seq <- paste(PURES[1],N,sep = ";")
        }
        
      } else {
        final_seq <- paste(SEQ,1,sep = ";")
      }
    } else {
      final_seq <- paste(SEQ,1,sep = ";")
    }
  } else{
    final_seq <- paste(SEQ,1,sep = ";")
  }
  return(final_seq)
}

# Main function
# It integrates most of previous functions
MainCorrection <- function(telo, PURE_MIDDLE){
  
  # Clean sample names
  telo$sample = gsub(".sorted.recal.md", "", telo$sample)
  telo$sample = gsub(".bam.sorted.recal.md", "", telo$sample)
  telo$sample = gsub("_rawFus_unmapped_python.fq", "", telo$sample)
  telo$sample = gsub(".bam", "", telo$sample)
  telo$sample = gsub(".sam", "", telo$sample)
  telo$sample = gsub(".cram", "", telo$sample)
  telo$sample = gsub(".final", "", telo$sample)
  telo$sample = gsub(".recal.md", "", telo$sample)
  
  # Intersect with telo data
  PURE_MIDDLE_SEQUENCES <- unique(paste(PURE_MIDDLE$Type, PURE_MIDDLE$Sequence, sep = ':'))
  
  #--------
  # Step 0
  #--------
  # Determine original fusion type
  telo$type <- apply(telo, 1, function(x) if (x['chr'] %in% c('2_endogenous','4_endogenous','9_endogenous')){
    x['chr']
  } else if (paste(x['orientation'],x['middle'],sep = ':') %in% PURE_MIDDLE_SEQUENCES | nchar(x['middle']) == 0){
    'Pure'
  } else{
    'Alternative'
  })
  
  #-------
  # Step1
  #-------
  # Clean middle sequences if breakpoints with perfect match 
  # with forward and reverse sequences are found at the ends of the sequence
  telo$middle0 <- apply(telo,1,function(x) if (x['type'] != 'Pure') {
    SimplifyMiddle(x['middle'], x['orientation'])
  } else {
    x['middle']
  }
  )
  
  # Find the closes pure to a middle sequence 
  # Allowing for one mismatch or 1 indel
  # Save the number of possible solutions for correction 
  # For example, when 2 pure sequences are equally possible
  telo$middle1 <- apply(telo,1,function(x) MiddleCorrectionEndoPure(x['middle0'],x['orientation'],x['chr'],PURE_MIDDLE))
  telo$possible_pure_solutions1 <- apply(telo,1,function(x) PossibleChoicesPureCorrection(x['middle0'],x['orientation'],x['chr'],PURE_MIDDLE))
  
  telo$type1 <- apply(telo, 1, function(x) if (x['chr'] %in% c('2_endogenous','4_endogenous','9_endogenous')){
    x['chr']
  } else if (paste(x['orientation'],x['middle1'],sep = ':') %in% PURE_MIDDLE_SEQUENCES | nchar(x['middle1']) == 0){
    'Pure'
  } else{
    'Alternative'
  })
  
  #-------
  # Step 2
  #-------
  # Find fwd|rev patters allowing one error starting sequence from left
  telo$middle_left <- apply(telo,1,function(x) if (x['type1'] != 'Pure') {
    MiddleCorrectLeft(x['middle1'])
  } else {
    x['middle1']
  })
  
  # Find fwd|rev patters allowing one error starting sequence from right
  telo$middle_right <- apply(telo,1,function(x) if (x['type1'] != 'Pure') {
    MiddleCorrectRight(x['middle1'])
  } else {
    x['middle1']
  })
  
  # Clean middle sequences (remove fwd and rev repeats at the ends of the fragments (plus small insertions)
  telo$middle_raw_clean <- CleanMiddle(telo$middle)
  telo$middle_left_clean <- CleanMiddle(telo$middle_left)
  telo$middle_right_clean <- CleanMiddle(telo$middle_right)
  
  # Select consensus sequence from the cleaned solution
  # We select the shortest
  telo$middle2 <- apply(telo, 1, function(x) if (x['chr'] != '9_endogenous'){
    if (x['type1'] != 'Pure'){
      a <- c(x['middle_raw_clean'], x['middle_left_clean'], x['middle_right_clean'])
      a[nchar(a)==min(nchar(a))][1]
    } else {
      x['middle1']
    }
  } else{
    a <- c(x['middle_raw_clean'], x['middle_left_clean'], x['middle_right_clean'])
    if (4 %in% nchar(a)){
      a[nchar(a)==4][1]
    } else {
      a[nchar(a)==min(nchar(a))][1]
    }
  }
  )
  
  # Re-call fusion types
  telo$type2 <- apply(telo, 1, function(x) if (x['chr'] %in% c('2_endogenous','4_endogenous','9_endogenous')){
    x['chr']
  } else if (paste(x['orientation'],x['middle2'],sep = ':') %in% PURE_MIDDLE_SEQUENCES | nchar(x['middle2']) == 0){
    'Pure'
  } else{
    'Alternative'
  })
  
  #--------
  # Step 3
  #--------
  # Correct endo and pure sequences after the correction step 2
  telo$middle3 <- apply(telo,1,function(x) MiddleCorrectionEndoPure(x['middle2'],x['orientation'],x['chr'],PURE_MIDDLE))
  telo$possible_pure_solutions3 <- apply(telo,1,function(x) PossibleChoicesPureCorrection(x['middle2'],x['orientation'],x['chr'],PURE_MIDDLE))
  
  # Re-call fusion types
  telo$type3 <- apply(telo, 1, function(x) if (x['chr'] %in% c('2_endogenous','4_endogenous','9_endogenous')){
    x['chr']
  } else if (paste(x['orientation'],x['middle3'],sep = ':') %in% PURE_MIDDLE_SEQUENCES | nchar(x['middle3']) == 0){
    'Pure'
  } else{
    'Alternative'
  })
  
  #--------
  # Step 4
  #--------
  # Correct middle islands between interspersed __fw__ISLAND__rv__ patterns
  telo$middle_island <- apply(telo,1,function(x) IslandCorrection(x['middle3']))
  telo$middle_island_clean <-  CleanMiddle(telo$middle_island)
  telo$middle4 <- apply(telo,1,function(x) MiddleCorrectionEndoPure(x['middle_island_clean'],x['orientation'],x['chr'],PURE_MIDDLE))
  telo$possible_pure_solutions4 <- apply(telo,1,function(x) PossibleChoicesPureCorrection(x['middle_island_clean'],x['orientation'],x['chr'],PURE_MIDDLE))
  
  # Re-call fusion type
  telo$type4 <- apply(telo, 1, function(x) if (x['chr'] %in% c('2_endogenous','4_endogenous','9_endogenous')){
    x['chr']
  } else if (paste(x['orientation'],x['middle4'],sep = ':') %in% PURE_MIDDLE_SEQUENCES | nchar(x['middle4']) == 0){
    'Pure'
  } else{
    'Alternative'
  })
  
  #--------
  # Step 5
  #--------
  # Simplify middle sequence when we find a breakpoint (fw-SOMETHING-rv or rv-SOMETHING-fw in the middle sequence)
  telo$middle_simplified <- apply(telo,1,function(x) if (x['type4'] != 'Pure') {
    SimplifyMiddle(x['middle4'], x['orientation'])
  } else {
    x['middle4']
  }
  )
  
  # Correct endo-pure sequences
  telo$middle5 <- apply(telo,1,function(x) MiddleCorrectionEndoPure(x['middle_simplified'],x['orientation'],x['chr'],PURE_MIDDLE))
  telo$possible_pure_solutions5 <- apply(telo,1,function(x) PossibleChoicesPureCorrection(x['middle_simplified'],x['orientation'],x['chr'],PURE_MIDDLE))
  
  # Re-call fusion types
  telo$type5 <- apply(telo, 1, function(x) if (x['chr'] %in% c('2_endogenous','4_endogenous','9_endogenous')){
    x['chr']
  } else if (paste(x['orientation'],x['middle5'],sep = ':') %in% PURE_MIDDLE_SEQUENCES | nchar(x['middle5']) == 0){
    'Pure'
  } else{
    'Alternative'
  })
  
  #--------
  # Step 6
  #--------
  # Re-check pure sequences
  A <- PURE_MIDDLE
  A$middle_left <- apply(A,1,function(x) MiddleCorrectLeft(x['Sequence']))
  A$middle_right <- apply(A,1,function(x) MiddleCorrectRight(x['Sequence']))
  A$middle_left_clean <- CleanMiddle(A$middle_left)
  A$middle_right_clean <- CleanMiddle(A$middle_right)
  A <- unique(A[,c('Type','Sequence',"middle_left","middle_right","middle_left_clean","middle_right_clean")])
  
  # Pure fusions
  telo$middle_final_pure_correction <- NA
  Pure <- telo[telo$type5 == 'Pure',]
  if (nrow(Pure) > 0){
    for (Row in seq(1, nrow(Pure))){
      Temp <- Pure[Row,]
      SEQ <- Temp$middle5
      orientation <- Temp$orientation
      middle_0 <- Temp$middle0
      middle_left <- Temp$middle_left
      TYPE <- Temp$type5
      
      RESULTS <- unique(CorrectFinalPure(SEQ, orientation,middle_0,middle_left,A))
      Pure[Row,]$middle_final_pure_correction <- RESULTS
    }
  }
  
  # Not pure fusions
  Not_pure <- telo[telo$type5 != 'Pure',]
  if (nrow(Not_pure) > 0){
    Not_pure$middle_final_pure_correction <- with(Not_pure,paste(middle5,1,sep = ';'))
  }
  telo <- rbind(Pure,Not_pure)
  telo <- telo %>% separate(middle_final_pure_correction, c("middle6", "possible_pure_solutions6"),sep = ';')
  
  
  # Change type of short alternative middles to Pure_type2
  telo$type6 <- telo$type5
  
  # Add pure_type2 category to the fusion types
  # New type
  telo$type6 <- apply(telo, 1, function(x) if (x['type6'] == 'Alternative' & nchar(x['middle6']) <= 12) {
    'Pure_type2'
  } else{
    x['type6']
  })
  # Original type
  telo$type <- apply(telo, 1, function(x) if (x['type'] == 'Alternative' & nchar(x['middle']) <= 12) {
    'Pure_type2'
  } else{
    x['type']
  })
  
  # Collapse possible solution info
  telo$possible_pure_solutions6 <- apply(telo,1,function(x) max(as.numeric(x['possible_pure_solutions1']),
                                                                as.numeric(x['possible_pure_solutions2']), 
                                                                as.numeric(x['possible_pure_solutions3']),
                                                                as.numeric(x['possible_pure_solutions4']),
                                                                as.numeric(x['possible_pure_solutions5']),
                                                                as.numeric(x['possible_pure_solutions6']), na.rm = T))
  
  if (min(telo$possible_pure_solutions6) < 0){
    telo[telo$possible_pure_solutions6 < 0,]$possible_pure_solutions6 <- NA
  }
  
  # Check if the obtained sequence is corrected or not
  telo$n_with_correction <- apply(telo, 1, function(x) if (x['middle'] != x['middle6']){
    as.numeric(x['n'])
  } else {
    0
  })
  
  # Check if it was possible more than one solution
  telo$n_with_more_than_1_solution <- apply(telo, 1, function(x) if (!is.na(x['possible_pure_solutions6']) & as.numeric(x['possible_pure_solutions6']) > 1){
    as.numeric(x['n'])
  } else {
    0
  })
  
  return(telo)
}

# Search for blacklist samples
SearchBlacklistSamples <- function(telo_corrected, cutoff,min_endo9,prefix){
  # Get unique samples
  SAMPLES <- unique(telo_corrected[,c('sample','project')])
  SAMPLES$Count <- 1
  
  # Check only the endo9 sequences, where TTAA sequence is expected
  ENDO <- telo_corrected[telo_corrected$chr ==  '9_endogenous',]
  if (nrow(ENDO) > 0){
    # Check TTAA sequence
    ENDO$n_TTAA_ENDO <- apply(ENDO, 1, function(x) if (x['chr'] == '9_endogenous' & x['orientation'] == 'end-to-head' & x['middle6'] == 'TTAA'){
      as.numeric(x['n'])
    } else {
      0
    })
    
    # Count the number of endo9 counts per sample
    ENDO_AGG <- aggregate(cbind(n,n_TTAA_ENDO) ~ sample + project, data = ENDO, sum)
    ENDO_AGG$ENDO_CORRECT <-  ENDO_AGG$n_TTAA_ENDO/ENDO_AGG$n
    ENDO_AGG <- merge(ENDO_AGG, SAMPLES, by = c('sample','project'), all = T)
    ENDO_AGG[is.na(ENDO_AGG)] <- 0
    if (min(ENDO_AGG$n) < 1){
      ENDO_AGG[ENDO_AGG$n == 0,]$ENDO_CORRECT <- NA 
    }
    ENDO_AGG$sample <- factor(ENDO_AGG$sample, levels=ENDO_AGG[order(ENDO_AGG$ENDO_CORRECT, decreasing = F),]$sample)
    
    # Mark blacklist samples
    BLACKLIST <- ENDO_AGG[(ENDO_AGG$ENDO_CORRECT < cutoff | ENDO_AGG$n < min_endo9) & !ENDO_AGG$project %in% c('Maciejowski_2020_NatGen','Nat2020','Cleal2019','Cell2015'),]
    telo_corrected$blacklist <- (telo_corrected$sample %in% BLACKLIST$sample)*1
    
    # General view of endo9 TTAA proportion distribution
    Proportion_ENDO_CORRECT_PLOT <- ggplot(ENDO_AGG, aes (x = sample, y = ENDO_CORRECT)) + 
      geom_bar(stat = 'identity', alpha = 0.80, size = 0.1) + 
      ylim(0,1)+
      geom_hline(yintercept=cutoff, linetype="dashed", color = "red") +
      common_ggplot2 +
      labs (y = 'Proportion of TTAA (after correction) middle sequences', x = 'Samples') +
      theme(axis.text.x=element_blank())
    
    # General view of endo9 TTAA proportion distribution (top 100 samples)
    Proportion_ENDO_CORRECT_TOP100_PLOT <- ggplot(ENDO_AGG[ENDO_AGG$sample %in% head(ENDO_AGG$sample, 100) & ENDO_AGG$n >= 5,], aes (x = sample, y =ENDO_CORRECT, fill = project)) + 
      geom_bar(stat = 'identity', color = 'white', width = 0.80) + 
      common_ggplot2 +
      ylim(0,1)+
      geom_hline(yintercept=cutoff, linetype="dashed", color = "red") +
      labs (y = 'Proportion of TTAA (after correction) middle sequences', x = 'Samples') +
      theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5),
            plot.margin = margin(10, 10, 10, 10)) 
    
    # Out plot 1
    outplot1 <- paste(prefix, '.proportion_correct_endo9.all_samples.pdf', sep = '')
    ggplot2::ggsave(plot = Proportion_ENDO_CORRECT_PLOT, filename = outplot1, device = 'pdf', useDingbats=FALSE, height = 5, width = 10)
    
    # Out plot 2
    outplot2 <- paste(prefix, '.proportion_correct_endo9.top100.pdf', sep = '')
    ggplot2::ggsave(plot = Proportion_ENDO_CORRECT_TOP100_PLOT, filename = outplot2, device = 'pdf', useDingbats=FALSE, height = 10, width = 15)
  }
  
  return(telo_corrected)
}

# Common ggplot style
common_ggplot2 <- theme_bw() + theme(axis.text.x=element_text(size=8,angle=0,colour="black"),
                                     axis.text.y=element_text(size=8,colour="black"),
                                     axis.title.y=element_text(size=8,colour="black"),
                                     # axis.ticks.y = element_blank(),
                                     axis.title.x=element_text(size=8,colour="black"),
                                     legend.title = element_blank(),#text(size=16),
                                     legend.text = element_text(size=8,colour="black"),
                                     legend.position = "right",
                                     legend.key = element_blank(),
                                     axis.ticks= element_line(size=0.2), #0.5 is defect size
                                     plot.margin=unit(c(0,0.1,0.1,.1),"cm"),
                                     plot.title=element_text(size=10,colour="black")) +
  theme(    legend.position="bottom",
            legend.justification = c(0,0.5),
            text=element_text(colour="black"),
            legend.title.align = 0,
            legend.title=element_text(size=7,colour="black"),
            legend.key.size = unit(0.5,"line")) +
  theme(#panel.border = element_blank(),
    plot.background = element_blank()
  )

#---------------------------------------
#---------------------------------------
# Main computation
# Collapse and error correction fusion script
#---------------------------------------
#---------------------------------------

# 1. Arguments
parser <- ArgumentParser()

# setting parameters
parser$add_argument("-summary_file_collapsed", "--summary_file_collapsed", type="character", help="Summary file with fusions", metavar="file", required=T)
parser$add_argument("-project", "--project", default='Project',help="Project Id to be used", nargs=1, required=F)
parser$add_argument("-prefix", "--prefix", default = NULL, help="Prefix id of the out files. It can contain the folder path", nargs=1, required=F)
parser$add_argument("-endo9_prop", "--endo9_prop", default = 0.9, help="Minimum proportion of TTAA reads found in endo 9 region [Default: 0.9]", nargs=1, required=F)
parser$add_argument("-min_endo9", "--min_endo9", default = 5, help="Minimum number of reads mapping in endo9 region [Default: 5]", nargs=1, required=F)

# 0. Reading parameters
args <- parser$parse_args()

file <- args$summary_file
project <- args$project
prefix <- args$prefix
cutoff <- args$endo9_prop
min_endo9 <- args$min_endo9

if (is.null(prefix)){
  prefix <- paste('.',project, sep = '/')
}

# 1. Create out folder
outf <- dirname(prefix)
dir.create(outf, showWarnings = FALSE)

# 2. Create possible pure middle sequences
PURE_MIDDLE <- CreatePureMiddleTable()
# Save table with possible pure middle sequences and where the fusion-cut occurred
outfile0 <- paste(outf, 'Possible_middle_sequences.pure.tsv', sep = '/')
write.table(x = PURE_MIDDLE, file = outfile0,sep = '\t', quote = F, row.names = F, col.names = T)

# 3. Perform fusion correction
telo <- as.data.frame(fread(file))
if (nrow(telo) > 0){
  # 4. Main correction
  telo_corrected <- MainCorrection(telo,PURE_MIDDLE)
  
  
  # 5. Search for blacklist samples
  telo_corrected2 <- SearchBlacklistSamples(telo_corrected = telo_corrected, 
                                            cutoff,
                                            min_endo9,
                                            prefix)
  
  
  # 6. Save blacklist samples
  if ("blacklist" %in% colnames(telo_corrected2)){
    BLACKLIST <- unique(telo_corrected2[telo_corrected2$blacklist == 1, c('sample','project')])
    outfile1 <- paste(prefix, '.blacklisted_samples.tsv', sep = '')
    write.table(BLACKLIST, file = outfile1, sep = '\t', col.names = T, row.names = F, quote = F)
  } else {
    telo_corrected2$blacklist <- 0
  }
  
  # 7. Save step corrections
  outfile2 <- paste(prefix, '.middle_correction_steps.tsv', sep = '')
  write.table(x = telo_corrected2, file = outfile2,sep = '\t', quote = F, row.names = F, col.names = T)
  
  # 8. Save collapsed fusions
  collapsed_fusion <- aggregate(cbind(n,n_with_correction,n_with_more_than_1_solution) ~ chr + sample + project + type_fus + orientation + type6 + Subtype + middle6  + blacklist, data = telo_corrected2, sum) 
  colnames(collapsed_fusion) <- c("chr", "sample","project", "region_type", "orientation","type", "Subtype","middle","blacklist","n","n_with_correction","n_with_more_than_1_solution")
  collapsed_fusion$middle_len <- nchar(collapsed_fusion$middle)
  
  outfile3 <- paste(prefix, '.corrected.tsv', sep = '')
  write.table(x = collapsed_fusion, file = outfile3,sep = '\t', quote = F, row.names = F, col.names = T)
} else {
  print(paste0(file, ' is empty'))
}
