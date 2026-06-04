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
#' @param maxaccepts maxaccepts parameter in vsearch
#' @param maxrejects maxrejects parameter in vsearch
#' @return output the blast output to the specified path
#' @export
alignSequences <- function(repSeqPath, refPath, out, pctid=0.97, threads=1,
    maxaccepts=8, maxrejects=64) {
	res <- system(paste("vsearch --usearch_global",
    repSeqPath,
    "--db",
    refPath,
    "--id",
    pctid,
    '--top_hits_only',
    '--maxaccepts', maxaccepts,
    '--maxrejects', maxrejects,
    '--uc_allhits',
    '--blast6out',
    out,
    "--threads",
    threads))
    return(res)
}