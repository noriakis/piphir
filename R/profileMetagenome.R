#' profileMetagenome
#' @param KOTable precalculated KO table
#' @param blastRes output of alignSequences function (vsearch)
#' The second column must contain representative taxon ID linked to KO table.
#' @export
#' @return KO profile table
profileMetagenome <- function(taxTable, copyNumTable, KOTable, blastRes) {
    totalRead <- sum(taxTable)
    # blastRes$V2 <- blastRes$V2 %>% strsplit("\\|") %>% vapply("[", 1, FUN.VALUE="a")
    ASVs <- blastRes$V1 %>% unique()
    hitRes <- lapply(ASVs, function(asv) {
        tmp <- blastRes[blastRes$V1 == asv, ]
        hit <- unique(tmp$V2)
        list(hit, length(hit))
    })
    names(hitRes) <- ASVs
    asvCount <- lapply(hitRes, function(x) x[[2]]) %>% unlist()
    
    hitRes <- hitRes[lapply(hitRes, function(x) x[[2]]!=0) %>% unlist()]
    hitLen <- sum(lapply(hitRes, function(x) x[[2]]) %>% unlist())
    hitTaxTable <- taxTable[intersect(row.names(taxTable), names(hitRes)), ]
    hitTaxTableCounts <- sum(hitTaxTable)

    # Divide ASV read count by number of top vsearch hits
    divided <- hitTaxTable / asvCount[row.names(hitTaxTable)]
    
    ## Make taxa / abundance table for RefSeq
    convertTable <- do.call(rbind, lapply(row.names(divided), function(asv) {
        do.call(rbind, lapply(hitRes[[asv]][[1]], function(taxa) {
            c(divided[asv,], taxa)
        }))
    })) %>% data.frame()
    
    colnames(convertTable) <- c(colnames(convertTable)[1:(ncol(convertTable)-1)] , "ID")

    
    convertTable <- convertTable %>%
        mutate_at(1:(ncol(convertTable)-1), as.numeric) %>%
        group_by(ID) %>%
        summarize_at(1:(ncol(convertTable)-1), sum)
    convertTable$ID <- convertTable$ID %>% unlist()
    convertTable[["copynum"]] <- copyn[convertTable$ID,]
    normalized <- convertTable %>%
        mutate_at(2:(ncol(convertTable)-1), function(x) x / convertTable$copynum)

    keggpSubset <- keggp[normalized$ID, ]
    keggpSubset[is.na(keggpSubset)] <- 0
    keggpSubset <- as.matrix(keggpSubset)
    normalized$ID <- NULL
    normalized$copynum <- NULL
    normalized <- as.matrix(normalized)
    normKO <- t(normalized) %*% keggpSubset

    return(normKO)
}