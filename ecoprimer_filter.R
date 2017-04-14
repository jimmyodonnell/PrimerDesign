#!/usr/bin/env Rscript

#-------------------------------------------------------------------------------
# Author: Jimmy O'Donnell <jodonnellbio@gmail.com>
#-------------------------------------------------------------------------------

# test function:
test <- function(x) ifelse(x, "pass", stop("test failed!"))

# some test strings
lengths <- 1:200

seqs <- vector()
for(i in 1:length(lengths)){
  tmp <- sample(LETTERS, size = lengths[i], replace = TRUE)
  seqs[i] <- paste0(tmp, collapse = "")
}

primer_length <- function(x, min = 18, max = 30) {
  ( nchar(x) > min & nchar(x) < max )
}

test(identical(which(primer_length(seqs)), 19:29))

amplicon_length <- function(x, min = 70, max = 150){
  ( nchar(x) > min & nchar(x) < max )
}

test(identical(which(amplicon_length(seqs)), 71:149))

gc_content <- function(gc_count, primer_len, min = 0.4, max = 0.6){
  gc_prop <- gc_count/primer_len
  ( gc_prop >= min & gc_prop <= max )
}

test(identical(which(sapply(1:20, gc_content, primer_len = 20)), 8:12))

repeats <- function(x, max = 4){
  the_regex <- paste0('(\\w)\\1{', max, ',}')
  !grepl(the_regex, x = x)
}

mystrings <- c("atcg", "aychaydgg", "diehkaaa", "djhfossss", "fjrndttttt", "ggggggsjfiejcndjskeid")

test(identical(which(repeats(mystrings)), 1:4))

melt_temp <- function(x, min = 50, max = 65){
  ( x >= min & x <= max )
}

test(identical(which(melt_temp(1:100)), 50:65))

melt_temp_diff <- function(x, y, max = 2){
  abs(x - y) <= 2
}

test(identical(which(melt_temp_diff(rep(5, 10), 1:10)), 3:7))

write.fasta <- function(vec, outfile)
# write a fasta file from a named character vector
{
  if(file.exists(outfile)){
    stop("fasta file already exists!")
  }
  for(i in 1:length(vec)){
    write(paste0(">", names(vec)[i]), file = outfile, append = TRUE)
    write(vec[i], file = outfile, append = TRUE)  
  }
}

# it turns out this is unnecessary because orientation is not considered anyway
revcomp <- function(str)
{
  orig <- toupper(str)
  orig_v <- unlist(strsplit(orig, ""))
  # source: http://droog.gs.washington.edu/parc/images/iupac.html
  comp <- c(
    "A" = "T", "T" = "A", 
    "G" = "C", "C" = "G", 
    "Y" = "R", "R" = "Y", 
    "S" = "S", "W" = "W", 
    "K" = "M", "M" = "K", 
    "D" = "H", "H" = "D", 
    "B" = "V", "V" = "B", 
    "N" = "N" 
    )
  out <- paste0(rev(unname(comp[orig_v])), collapse = "")
  return(out)
}

# ------------------------------------------
# Table result description : 
# column  1 : serial number
# column  2 : primer1
# column  3 : primer2
# column  4 : primer1 Tm without mismatch
# column  5 : primer1 lowest Tm against exemple sequences
# column  6 : primer2 Tm without mismatch
# column  7 : primer2 lowest Tm against exemple sequences
# column  8 : primer1 G+C count
# column  9 : primer2 G+C count
# column 10 : good/bad
# column 11 : amplified example sequence count
# column 12 : amplified counterexample sequence count
# column 13 : yule
# column 14 : amplified example taxa count
# column 15 : amplified counterexample taxa count
# column 16 : ratio of amplified example taxa versus all example taxa (Bc index)
# column 17 : unambiguously identified example taxa count
# column 18 : ratio of specificity unambiguously identified example taxa versus all example taxa (Bs index)
# column 19 : minimum amplified length
# column 20 : maximum amplified length
# column 21 : average amplified length
# ------------------------------------------

library(data.table)
library(magrittr)
library(tools) # file_path_sans_ext()

the_file <- commandArgs(trailingOnly = TRUE)[1]

dt <- fread(the_file)

primer_len <- nchar(dt[1,2])

dt[
  V19 >= 70 & V19 <= 150 &
  V20 >= 70 & V20 <= 150 & 
  V21 >= 70 & V21 <= 150 & 
  gc_content(V8, primer_len) & gc_content(V9, primer_len) & 
  repeats(V2, max = 4) & repeats(V3, max = 4) & 
  melt_temp(V4) & melt_temp(V5) & melt_temp(V6) & melt_temp(V7) & 
  melt_temp_diff(V4, V6) & melt_temp_diff(V5, V7) & 
  V10 %in% "GG"
] -> primers_filtered 

primer_seqs <- unique(c(primers_filtered$V2, primers_filtered$V3))
names(primer_seqs) <- paste0("primer-", sprintf("%03d", 1:length(primer_seqs)))

# write fasta file
primer_file <- paste0(file_path_sans_ext(the_file), "_primers.fasta")
write.fasta(primer_seqs, primer_file)

# write filtered table
outfile <- paste0(file_path_sans_ext(the_file), "_filt.csv")
fwrite(primers_filtered, file = outfile)
