read_gmt_hier <- function(fname){
    res <- list(genes=list(), desc=list())
    gmt <- file(fname)
    gmt_lines <- readLines(gmt)
    close(gmt)
    gmt_list <- lapply(gmt_lines, function(x) unlist(strsplit(x, split="\t")))
    gmt_names <- sapply(gmt_list, '[', 1)
    gmt_desc <- lapply(gmt_list, '[', 2)
    gmt_genes <- lapply(gmt_list, function(x){x[3:length(x)]})
    names(gmt_desc) <- names(gmt_genes) <- gmt_names
    res <- do.call(rbind, lapply(names(gmt_genes),
                                 function(n) cbind.data.frame(term=n, gene=gmt_genes[[n]], hier=gmt_desc[[n]], stringsAsFactors=FALSE)))
    res$term <- as.factor(res$term)
    path_desc <- as.data.frame(cbind(gmt_names, gmt_desc))
    return(list(res, path_desc))
}
