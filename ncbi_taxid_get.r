#!/usr/bin/env Rscript

library(taxize)

# get ncbi taxon id (taxid) for a taxon name
args <- commandArgs(trailingOnly = TRUE)
# args <- "Boraras brigittae"

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("You must supply a quoted taxon name as an argument!", call. = FALSE)
}

record <- get_ids(args, db = "ncbi", verbose = FALSE)

ncbi_taxid <- record$ncbi[[1]]

cat(paste0(ncbi_taxid, "\n"))
