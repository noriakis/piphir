#' loadExample
#' load example input data for profileMetagenome
#' 
#' @export
#' @return list of 16S copy number, KO copy number, BLAST result, sequence abundance
loadExample <- function() {
    seqtab <- readRDS(system.file("extdata", "example_seqtab.rds", package = "piphir"))
    kocn <- readRDS(system.file("extdata", "example_KOCN.rds", package = "piphir"))
    blast <- readRDS(system.file("extdata", "example_blast.rds", package = "piphir"))
    cn <- readRDS(system.file("extdata", "example_16SCN.rds", package = "piphir"))
    list("cn16s"=cn, "cnko"=kocn, "blast"=blast, "seqtab"=seqtab)
}