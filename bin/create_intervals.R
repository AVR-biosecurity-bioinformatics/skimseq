#!/usr/bin/env Rscript

tryCatch({

args <- commandArgs(trailingOnly = TRUE)

projectDir              <- args[1]
params.rdata            <- args[2]
params.interval_size    <- args[3]

sys.source(paste0(projectDir,"/bin/functions.R"), envir = .GlobalEnv)

### load only required packages
process_packages <- c(
    "dplyr",
    "GenomicRanges",
    "readr",
    NULL
)
invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

### run code

if ( params.interval_size == "null" ) {
    interval_size <- 5E6 # default value if parameter unset
} else {
    interval_size <- as.numeric(params.interval_size)
}

index_file <- list.files(pattern = "\\.fai$")

# Read in fasta index with lengths of assembly units (chromosome/scaffold/contig)
assembly_lengths <- 
    readr::read_tsv(
        index_file,
        col_names = c(
            "name", 
            "n_bases",
            "index", 
            "base_per_line",
            "byte_per_line")
    ) %>%
    dplyr::select(name, n_bases)

## split chromosomes/scaffolds/contigs into 5Mb chunks (or a threshold)

seqlengths <- assembly_lengths$n_bases
names(seqlengths) <- assembly_lengths$name
tiles <- GenomicRanges::tileGenome(seqlengths, tilewidth = interval_size)

assembly_chunks <- 
    dplyr::bind_rows(
        tiles %>%
        as.data.frame() %>%
        dplyr::mutate(
            group_name = paste0("chunk_", group),
            # convert to 0-indexed base positions for .bed format
            start = start - 1,
            end = end - 1    
        )
    ) %>%
    dplyr::mutate(interval=paste0(seqnames, ":",start,"-",end))

# Write out intervals
assembly_chunks %>%
    dplyr::group_by(group_name) %>%
    dplyr::summarise(interval = paste0(sort(interval), collapse = ",")) %>%
    dplyr::select(interval) %>%
    readr::write_tsv(file="intervals.txt", col_names = FALSE)

}, 
finally = {
    ### save R environment if script throws error code
    if (params.rdata == "true") {save.image(file = "CREATE_INTERVALS.rda")}
})