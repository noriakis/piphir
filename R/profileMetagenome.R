#' profileMetagenome
#' 
#' the main workflow
#' 
#' @param taxTable taxonomy abundance table
#' @param copyNumTable 16S copy number table (row.names corresponds to reference tax ID)
#' @param KOTable precalculated KO table
#' @param blastRes output of alignSequences function (vsearch)
#' The second column must contain representative taxon ID linked to KO table.
#' @param full return the full output, default to FALSE
#' @param fullOutput "file" or "inmemory"
#' @param fullTemp if fullOutput is "file", {{fileTemp}}.txt will be created.
#' @param fullID subset to this ID in full mode, when the data size is large.
#' @param forceCalc force calculation when the index of normalized hit table and
#' KO copy number table does not match, by removing blast hits that not in CN table.
#' @export
#' @import dplyr tidyr
#' @return KO profile table
profileMetagenome <- function(taxTable, copyNumTable, KOTable, blastRes, full=FALSE,
    fullOutput="file", fullTemp="temporary_full", fullID=NULL, forceCalc=TRUE) {
    totalRead <- sum(taxTable)
    # blastRes$V2 <- blastRes$V2 %>% strsplit("\\|") %>% vapply("[", 1, FUN.VALUE="a")
    ASVs <- blastRes[,1] %>% unique()
    hitRes <- lapply(ASVs, function(asv) {
        tmp <- blastRes[blastRes[,1] == asv, ]
        hit <- unique(tmp[,2])
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
    convertTable[["copynum"]] <- copyNumTable[convertTable$ID, 1]
    normalized <- convertTable %>%
        mutate_at(2:(ncol(convertTable)-1), function(x) x / convertTable$copynum)
    conv <- data.frame(normalized)
    if (forceCalc) {
        nonnaID <- normalized[!is.na(normalized$copynum), ]$ID
        forceID <- intersect(row.names(KOTable), nonnaID)
        KOTable <- KOTable[forceID, ]
        normalized <- normalized %>% data.frame() %>% `row.names<-`(.[["ID"]])
        normalized <- normalized[forceID, ]
        keggpSubset <- KOTable[normalized$ID, ]
        row.names(conv) <- conv$ID
        conv <- conv[forceID, ]
    } else {
        keggpSubset <- KOTable[normalized$ID, ]
    }
    keggpSubset[is.na(keggpSubset)] <- 0
    keggpSubset <- as.matrix(keggpSubset)
    normalized$ID <- NULL
    normalized$copynum <- NULL
    normalized <- as.matrix(normalized)
    normKO <- t(normalized) %*% keggpSubset
    
    if (!full) {
	    ret <- t(normKO)
	    return(ret)	
    }
    
    conv$copynum <- NULL

    longDf <- tidyr::pivot_longer(conv, 2:ncol(conv)) %>% data.frame()
    # Cannot handle this big size data in memory
    if (fullOutput=="inmemory") {
        strat <- do.call(rbind, lapply(seq_len(nrow(longDf)), function(rn) {
            kos <- names(keggpSubset[longDf[rn, "ID"], ])
            if (!is.null(fullID)) {
                if (length(intersect(fullID, kos))>=1) {} else {return(NULL)}
            }
            tmp <- data.frame(keggpSubset[longDf[rn, "ID"], ] * longDf[rn, "value"])
            colnames(tmp) <- c("value")
            tmp[["ID"]] <- longDf[rn, "ID"]
            tmp[["sample"]] <- longDf[rn, "name"]
            tmp[["KO"]] <- row.names(tmp)
            if (!is.null(fullID)) {
                if (length(intersect(fullID, tmp$KO))>=1) {} else {return(NULL)}
            }
            return(tmp)
        }))
        if (is.null(strat)) {stop("No KO profiled.")}
        return(tidyr::pivot_wider(strat, names_from=sample))     
    } else {
        fileName <- paste0(fullTemp, ".txt")
        if (file.exists(fileName)) {stop("File exists!")}
        lapply(seq_len(nrow(longDf)), function(rn) {

          tmp <- data.frame(keggpSubset[longDf[rn, "ID"], ] * longDf[rn, "value"])
          colnames(tmp) <- c("value")
          tmp[["ID"]] <- longDf[rn, "ID"]
          tmp[["sample"]] <- longDf[rn, "name"]
          tmp[["KO"]] <- row.names(tmp)
          tmp <- tmp[tmp[["value"]]!=0, ]

          if (!is.null(fullID)) {
            if (length(intersect(fullID, tmp$KO))>=1) {} else {return(NULL)}
          }         
          write.table(tmp, fileName, col.names = !file.exists(fileName), sep="\t", append=TRUE,
            quote=FALSE, row.names=FALSE)
          return(0)
        })
        return(fileName)
    }
}