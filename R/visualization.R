#' @import plotly
#' @import iheatmapr

.create_heatmap <- function(cut_path_list, direction, term2gene, path_desc, p_cut){
    x <- .get_heatmap(cut_path_list, direction)
    cl <- .get_count_per_path(term2gene, path_desc)
    idx <- sapply(rownames(x), function(path){
        which(cl$Pathway == path)
    })
    idx <- unlist(idx)
    cl_path <- unlist(cl$gmt_desc[idx])
    ngene <- as.numeric(cl$Genes_per_pathway[idx])

    nsig <- .get_count_sig_path(x, p_cut)
    h <- .get_iheatmapr(x, cl, idx, cl_path, ngene, nsig, max_groups=10)
    return(h)
}

.get_heatmap <- function(cut_path_list, direction){
    path_df <- data.matrix(cut_path_list[[direction]]$path)
    return(path_df)
}

.get_count_per_path <- function(term2gene, path_desc){
    genes_per_path <- table(term2gene$term)
    res <- as.data.frame(cbind(names(genes_per_path), as.numeric(genes_per_path)))
    names(res) <- c("Pathway", "Genes_per_pathway")

    res <- merge(res, path_desc, by.x="Pathway", by.y="gmt_names")
    return(res)
}

.get_count_sig_path <- function(x, p_cut){
    nsig <- apply(x, 2, function(col){
        sum(col < p_cut)
    })
    return(nsig)
}

.get_iheatmapr <- function(x, cl, idx, cl_path, ngene, nsig, max_groups=10){
    h <-
        main_heatmap(x, name = "Enriching<br>value") %>%
        add_col_labels(font = list(size = 6)) %>%
        add_row_title("Pathways") %>%
        add_col_title("<b>Pathways enriched to different gene set sizes</b>", side = "top")
    if(length(unique(cl_path)) <= max_groups){
        h <- add_row_groups(cl_path, name = "Pathways<br>category")
    }

    h <- h %>%
        add_row_barplot(x = ngene,
                        tracename = "Pathway size",
                        layout = list(title = "Number of<br>genes")) %>%
        add_col_plot(y = nsig,
                     tracename = "Number of significantly<br>enriched pathways",
                     layout = list(title = "n sig. enrich. paths"),
                     side = "bottom")
    return(h)
}

