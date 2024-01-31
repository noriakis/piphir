#' alignSequences
#' 
#' Use vsearch to align representative sequences to reference sequence.
#' Will take time.
#' 
#' @param repSeqPath path to representative sequences
#' @param refPath path to reference 16S sequences
#' @param out output file of vsearch
#' @param pctid percent ID, vsearch parameter
#' @param threads thread number
#' @return output the blast output to the specified path
#' @export
alignSequences <- function(repSeqPath, refPath, out, pctid=0.97, threads=1) {
	res <- system(paste("vsearch --usearch_global",
    repSeqPath,
    "--db",
    refPath,
    "--id",
    pctid,
    '--top_hits_only',
    '--maxaccepts', '0',
    '--maxrejects', '0',
    '--uc_allhits',
    '--blast6out',
    out,
    "--threads",
    threads))
    return(res)
}