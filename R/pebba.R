#' @import gplots
#' @import clusterProfiler
NULL

#' Run Enricher
#'
#' This function runs enricher
#'
#' @param top_genes A vector of genes
#' @param all_genes An object with all genes
#' @param term2gene A data.frame with enrichment term and genes
#'
#' @keywords internal

.run_enrich <- function(top_genes, all_genes, term2gene){
    enriched <- as.data.frame(clusterProfiler::enricher(gene = top_genes,
                                    pvalueCutoff = 1,
                                    minGSSize = 1,
                                    universe = all_genes,
                                    TERM2GENE = term2gene,
                                    qvalueCutoff = 1,
                                    maxGSSize = 100))[, c(1, 6)]
    return(enriched)
}


#' Get cutoff value
#'
#' This function provides the cutoff value from PEBBA
#'
#' @param deg_list A list of DEGs
#' @param logFC_col A string indicating the column with log fold
#'           change values
#' @param pvalue_col A string indicating the column with p-values
#' @param min_genes Minimum number of genes
#' @param max_genes Maximum number of genes
#'
#' @keywords internal

.get_cutoff <- function(deg_list, logFC_col, pvalue_col, min_genes, max_genes){
    dirs <- c("down", "up")

    res <- lapply(dirs, function(direction){

    decreasing <- ifelse(direction == "down", FALSE, TRUE)

    top <- deg_list[head(order(deg_list[, logFC_col],
                               decreasing=decreasing),
                         n=max_genes),
                    c(logFC_col, pvalue_col)]
    #Add pi_value
    top$pi_value <- abs(top[, logFC_col]) * -log10(top[, pvalue_col])
    #Order pi_value
    top <- top[order(top$pi_value, decreasing=TRUE), ]
    df1 <- data.frame(matrix(0, nrow=0, ncol=4))
    for (i in seq(from=min_genes, to=max_genes, by=50)) {
        top_genes  <- top[1:i, ]
        minFC <- min(abs(top_genes[, 1]))
        maxP  <- max(top_genes[, 2])
        minP  <- -log10(maxP)
        minPi <- min(top_genes[i, 3])
        rowX  <- c(minFC, minP, minPi)
        df1 <- rbind(df1,rowX)
    }
    df1
    })
    names(res) <- dirs
    top_cut <- seq(from=min_genes, to=max_genes, by=50)
    res <- do.call("cbind", res)
    res <- cbind(top_cut, res)

    res$fc <- apply(res, 1, function(x) min(x[2], x[6]) )
    res$p  <- apply(res, 1, function(x) min(x[3], x[7]) )
    res$pi <- apply(res, 1, function(x) min(x[4], x[8]) )

    names(res) <- c("TopCut", "minimum_log2fc_down", "minimum_MinuslogP_down",
                    "minimum_Pi_down", "minimum_log2fc_up", "minimum_MinuslogP_up",
                    "minimum_Pi_up", "minimum_log2fc_combined",
                    "minimum_MinuslogP_combined", "minimum_Pi_combined")

    rownames(res) <- res[, 1]
    return(res)
}

#' Get pathway
#'
#' This function gets pathways.
#'
#' @param merge_p Something to merge
#' @param term2gene A data frame containing genes and terms
#' @param all_genes An object with all genes
#' @param deg_list A list with DEGs
#' @param gene_col A string indicating the column with genes
#' @param logFC_col A string indicating the column with log fold-change values
#' @param pvalue_col A string indicating the column with p-values
#' @param direction The direction. One of "up", "down" or "any"
#' @param min_genes Minimum number of genes
#' @param max_genes Maximum number of genes
#' @param p_cut P-value cut
#'
#' @keywords internal

.get_pathway <- function(merge_p, term2gene, all_genes, deg_list,
                        gene_col, logFC_col, pvalue_col, direction,
                        min_genes=50, max_genes=3000, p_cut=0.01){
    #Rank based on fold-change
    #direction = "up" or "down"
    #Top number of genes
    #top_n = "3000"
    #Minimum number of genes
    #min_genes = "50"
    #Max number of genes
    #max_genes = "3000"
    #Adj P-value cutoff for ORA results. Remove pathways not significant in any round
    #p_cut = 0.01
    #Order Merge Table by
    #order_p = "NG", "p", "P", "first", "times", "ES3"
    if(tolower(direction) == "up"){
        top <- deg_list[head(order(deg_list[, logFC_col], decreasing=TRUE), n=max_genes), ]
    }else if(tolower(direction) == "down"){
        top <- deg_list[head(order(deg_list[, logFC_col], decreasing=FALSE), n=max_genes), ]
    }else if(tolower(direction) == "any"){
        top <- deg_list[head(order(abs(deg_list[, logFC_col]), decreasing=FALSE), n=max_genes), ]
    }else{
        stop("Invalid direction argument")
    }

    # add pi_value
    top$pi_value = abs(top[, logFC_col])*-log10(top[, pvalue_col])

    # order pi_value
    top <- top[order(top$pi_value, decreasing=TRUE), ]

    for (i in seq(from=min_genes, to=max_genes, by=50)) {
        top_genes  <- as.character(top[1:i, gene_col])
          pathG <- .run_enrich(top_genes, all_genes, term2gene)
          colnames(pathG) <- c("term",  i)


          merge_p <- merge(merge_p,
                      pathG,
                      by="term",
                      all=TRUE)
          merge_p[is.na(merge_p)] <- 1
    }
    rownames(merge_p) <- merge_p[, 1]
    merge_p           <- merge_p[, -1]
    merge_p2          <- log10(merge_p)*-1

    path_cut_p <- log10(p_cut)*-1

    df <- data.frame(matrix(0, nrow(merge_p2), ncol=0))
    rownames(df) <- rownames(merge_p2)
    #top cut with maximum MinuslogP
    df$NG <- as.numeric(colnames(merge_p2)[apply(merge_p2, 1, which.max)])
    df$p  <- as.numeric(apply(merge_p2, 1, max))
    df$P  <- as.numeric(apply(merge_p2, 1, sum))

    #How many columns above path_cut_p (freq)
    df$times <- as.numeric(apply(merge_p2, 1, function(x) length(which(x > path_cut_p))))/ncol(merge_p2)

    #If the pathway has times > 0
    #First column above path_cut_p
    df$first <- apply(merge_p2, 1, function(x) ifelse (length(which(x > path_cut_p)) >0,
                                                as.numeric(colnames(merge_p2)[min(which(x > path_cut_p))]),
                                                0))

    #ES3
    df$ES3 <- (1 - exp(-df$p))/(1 + 0.1*sqrt(df$NG))

    #order
    merge_p2 <- merge_p2[order(df[, "first"], decreasing=TRUE), ]

    colnames(df) <- c(paste("TopCut_highestMinuslogP", "_", direction, sep=""),
                      paste("maximum_MinuslogP", "_", direction, sep=""),
                      paste("sum_MinuslogP", "_", direction, sep=""),
                      paste("times_significant", "_", direction, sep=""),
                      paste("FirstTopCut_significant", "_", direction, sep=""),
                      paste("PEBBA_score", "_", direction, sep=""))

    newList <- list("data.frame" = df, "data.frame" = merge_p2)
    return(newList)
}


#' Cutoff pathway
#'
#' This function takes a table and returns a dataframe with
#' several different types of information about pathways.
#'
#' @param path_table A table with some kind of information.
#' @param p_cut P-value cut.
#' @param direction The direction.
#'
#' @keywords internal

.cutoff_path <- function(path_table, p_cut, direction){
    df <- data.frame(matrix(0, nrow=ncol(path_table), ncol=0))
    rownames(df) <- colnames(path_table)

    df$MaxR  <- as.numeric(apply(path_table, 2, max))
    df$SumR  <- as.numeric(apply(path_table, 2, sum))
    path_cut_p <- log10(p_cut)*-1
    #How many pathways above path_cut_p (freq)
    df$times <- as.numeric(apply(path_table, 2,
                               function(x) length(which(x > path_cut_p))))/nrow(path_table)
    colnames(df) <- c(paste("maximum_MinuslogP", "_", direction, sep=""),
                    paste("sum_MinuslogP", "_", direction, sep=""),
                    paste("times_significant", "_", direction, sep=""))
    return(df)
}


