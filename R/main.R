#' @importFrom tools file_path_sans_ext
NULL

#' PEBBA analysis
#'
#' This function executes the PEBBA analysis.
#'
#' @param file_in The file or data.frame to execute the analysis on
#' @param gmt_file The name of the gmt file containing terms and genes
#' @param gene_col A string indicating the column with genes (Default: "Gene.symbol")
#' @param logFC_col A string indicating the column with log fold-change
#'        values (Default: "logFC").
#' @param pvalue_col A string indicating the column with p-values (Default: "P.Value")
#' @param min_genes Minimum number of genes (Default: 50)
#' @param max_genes Maximum number of genes (Default: 3000)
#' @param p_cut P-value cutoff (Default: 0.2)
#' @param verbose Logical. If TRUE (default), will display analysis progress messages.
#' @param analysis_name The name to give to analysis results.
#'        (If NULL, defaults to parameter \code{file_in} without the extensions and file path;
#'        if a data.frame, defaults to 'PEBBA_analysis' when left NULL)
#' @param results_dir The path into which results should be saved (Default: "Results").
#' @param force Whether or not to overrwrite an existing results directory (Default: FALSE).
#'
#' @return Tables and interactive heatmaps with PEBBA results
#'
#' @rdname pebba
#' @examples
#' # Run PEBBA analyses
#' data(example_data)
#' gmt_file <- system.file("extdata", "pathways.gmt", package = "PEBBA")
#' pebba(example_data, gmt_file)
#' @export

pebba <- function(file_in, gmt_file, gene_col="Gene.symbol",
                  logFC_col="logFC", pvalue_col="P.Value",
                  min_genes=100, max_genes=1500,
                  p_cut=0.2, verbose=TRUE,
                  analysis_name=NULL, results_dir="Results",
                  force=FALSE){

    # Validating inputs
    if(min_genes < 50 | min_genes > 2900){
        stop("Variable min_genes must be between 50 and 2900 genes")
    }
    if(max_genes < 100 | max_genes > 3000){
        stop("Variable max_genes must be between 100 and 3000 genes")
    }
    if(p_cut < 0.00001 | p_cut > 1){
        stop("Variable p_cut must be between 0.00001 and 1")
    }


    # Preparing files and workspace
    ## Disable scientifc notation
    options(scipen=999)

    ## Create a results directory
    if(dir.exists(results_dir)){
        if(!force){
            stop("Stopping analysis: ", results_dir,
                 " already exists! Use force=TRUE to overwrite.")
        }
    }else{
        dir.create(results_dir)
        dir.create(file.path(results_dir, "Tables"))
        dir.create(file.path(results_dir, "Heatmaps"))
    }

    if(is.null(analysis_name)){
        analysis_name <- "PEBBA_analysis"
    }

    ## Get information from all unique terms
    gmt_res <- read_gmt_hier(gmt_file)
    term2gene <- gmt_res[[1]]
    path_desc <- gmt_res[[2]]

    merge_p <- data.frame(unique(term2gene[1]))

    if(is.character(file_in)){
        deg_list <- read.csv(file_in, header = TRUE, sep = "\t")
        if(is.null(analysis_name)){
            analysis_name <- tools::file_path_sans_ext(basename(file_in))
        }
    }else if(is.data.frame(file_in)){
        deg_list <- file_in
    }

    ## Remove rows that do not have a valid gene symbol
    deg_list <- deg_list[which(deg_list[, gene_col]!=""), ]
    ## Make logFC and p-value columns numeric
    deg_list[, logFC_col] <- as.numeric(deg_list[, logFC_col])
    deg_list[, pvalue_col] <- as.numeric(deg_list[, pvalue_col])

    ## Get background genes as a character vector
    ## Empty values (non-annotated genes) will be removed
    all_genes <- as.character(deg_list[, gene_col])

    # Get cutoff values -------------------------------------------------------

    if(verbose) message("Getting cutoff")
    ## Get info about p-value and log2fc cutoff used on each top segments
    table_cut <- .get_cutoff(deg_list, logFC_col, pvalue_col, min_genes, max_genes)

    dirs <- c("up", "down", "any")

    cut_path_list <- lapply(dirs, function(direction){
        if(verbose) message(direction)
        if(verbose) message("Getting pathways")
        list_p <- .get_pathway(merge_p, term2gene, all_genes,
                            deg_list, gene_col, logFC_col,
                            pvalue_col, direction,
                            min_genes, max_genes, p_cut)

        df <- list_p[[1]]
        path <- list_p[[2]]

        if(verbose) message("Getting pathway cutoff")
        cut_path <- .cutoff_path(path, p_cut, direction)
        res <- list(cut_path, df, path)
        names(res) <- c("cut_path", "df", "path")
        res
    })
    names(cut_path_list) <- dirs

    ## Save heatmaps
    if(verbose) message("Saving heatmaps")
    sapply(dirs, function(direction){
        heat <- .create_heatmap(cut_path_list, direction, term2gene, path_desc, p_cut)
        f_out <- paste(analysis_name, "Pathways", paste0(direction, ".html"), sep="_")
        save_iheatmap(heat, filename=f_out)
        system(paste("mv", f_out, file.path(results_dir, "Heatmaps", f_out)))
    })

    # Results -----------------------------------------------------------------
    if(verbose) message("Combining results 1")
    df_combined <- .combine_cut(cut_path_up=cut_path_list$up$cut_path,
                               cut_path_down=cut_path_list$down$cut_path,
                               cut_path_any=cut_path_list$any$cut_path,
                               table_cut=table_cut)
    if(verbose) message("Combining results 2")
    df_c <- .combine_df(df_up   = cut_path_list$up$df,
                       df_down = cut_path_list$down$df,
                       df_any  = cut_path_list$any$df)

    if(verbose) message("Exporting data")
    .export_data(file_in=file_in, df_combined=df_combined, df_c=df_c,
                path_up=cut_path_list$up$path,
                path_down=cut_path_list$down$path,
                path_any=cut_path_list$any$path,
                analysis_name=analysis_name,
                results_dir=results_dir)
}

#' Combine cut results
#'
#' This function takes the cut_path_up, cut_path_down and
#' cut_path_any tables and returns the results tables.
#'
#' @param cut_path_up The cut_path_up table
#' @param cut_path_down The cut_path_down table
#' @param cut_path_any The cut_path_any table
#' @param table_cut The table_cut parameter
#'
#' @keywords internal
#' @noRd

.combine_cut <- function(cut_path_up, cut_path_down, cut_path_any, table_cut){
    #Combine cut_path_up and cut_path_down
    temp_df <- cbind(cut_path_up, cut_path_down, cut_path_any)
    temp_df$max <- apply(temp_df, 1, function(x) mean(x[1], x[4]) )
    temp_df$sum <- apply(temp_df, 1, function(x) mean(x[2], x[5]) )
    temp_df$times <- apply(temp_df, 1, function(x) mean(x[3], x[6]) )
    colnames(temp_df)[ncol(temp_df)-2] <- "maximum_MinuslogP_meanUPandDOWN"
    colnames(temp_df)[ncol(temp_df)-1] <- "sum_MinuslogP_meanUPandDOWN"
    colnames(temp_df)[ncol(temp_df)] <- "times_significant_meanUPandDOWN"
    df_combined <- cbind(table_cut, temp_df)
    return(df_combined)
}

#' Combine df_up and df_down
#'
#' This function combines the results of the df_up and df_down tables.
#'
#' @param df_up The df_up table
#' @param df_down The df_down table
#' @param df_any The df_any table
#'
#' @keywords internal
#' @noRd

.combine_df <- function(df_up, df_down, df_any){

    #Combine df_up and df_down
    df_c <- data.frame(matrix(0, nrow=nrow(df_up), ncol=0))
    rownames(df_c) <- rownames(df_up)
    df_c$Pathway   <- rownames(df_c)
    df_c$top   <- unlist(lapply(1:nrow(df_up), function(x) mean(df_up[x, 1], df_down[x, 1])))
    df_c$max   <- unlist(lapply(1:nrow(df_up), function(x) mean(df_up[x, 2], df_down[x, 2])))
    df_c$sum   <- unlist(lapply(1:nrow(df_up), function(x) mean(df_up[x, 3], df_down[x, 3])))
    df_c$times <- unlist(lapply(1:nrow(df_up), function(x) mean(df_up[x, 4], df_down[x, 4])))
    df_c$first <- unlist(lapply(1:nrow(df_up), function(x) mean(df_up[x, 5], df_down[x, 5])))
    df_c$fair  <- unlist(lapply(1:nrow(df_up), function(x) mean(df_up[x, 6], df_down[x, 6])))

    df_c$topMin   <- unlist(lapply(1:nrow(df_up), function(x) min(df_up[x, 1], df_down[x, 1])))
    df_c$maxMax   <- unlist(lapply(1:nrow(df_up), function(x) max(df_up[x, 2], df_down[x, 2])))
    df_c$sumMax   <- unlist(lapply(1:nrow(df_up), function(x) max(df_up[x, 3], df_down[x, 3])))
    df_c$timesMax <- unlist(lapply(1:nrow(df_up), function(x) max(df_up[x, 4], df_down[x, 4])))
    df_c$firstMin <- unlist(lapply(1:nrow(df_up), function(x) min(df_up[x, 5], df_down[x, 5])))
    df_c$fairMax  <- unlist(lapply(1:nrow(df_up), function(x) max(df_up[x, 6], df_down[x, 6])))

    colnames(df_c) <- c("Pathways",
                     "TopCut_highestMinuslogP_meanUPandDOWN",
                     "maximum_MinuslogP_meanUPandDOWN",
                     "sum_MinuslogP_meanUPandDOWN",
                     "times_significant_meanUPandDOWN",
                     "FirstTopCut_significant_meanUPandDOWN",
                     "PEBBA_score_meanUPandDOWN",
                     "TopCut_highestMinuslogP_minUPandDOWN",
                     "maximum_MinuslogP_maxUPandDOWN",
                     "sum_MinuslogP_maxUPandDOWN",
                     "times_significant_maxUPandDOWN",
                     "FirstTopCut_significant_minUPandDOWN",
                     "PEBBA_score_maxUPandDOWN")

    df_c <- cbind(df_c, df_up, df_down, df_any)
    return(df_c)
}


#' Export data
#'
#' This function takes the results of the analyses and returns txt files.
#'
#' @param file_in The input file
#' @param df_combined The result of the .combine_cut function
#' @param df_c The result of the .combine_df function
#' @param path_up Up pathways
#' @param path_down Down pathways
#' @param path_any Any pathways
#' @param analysis_name The name to be given to the analysis (Default: "Analysis1").
#' @param results_dir The folder into which results should be saved. A subdirectory
#'           called "Heatmaps" will be created inside this folder.
#'
#' @keywords internal
#' @noRd

.export_data <- function(file_in, df_combined, df_c, path_up,
                        path_down, path_any, analysis_name,
                        results_dir){
    # Exporting datasets

    f_out <- paste(analysis_name, "Cutoffs.txt", sep="_")
    write.table(df_combined, file = file.path(results_dir, "Tables", f_out),
                sep = "\t", row.names = FALSE)

    f_out <- paste(analysis_name, "PathwayMetrics.txt", sep="_")
    write.table(df_c, file = file.path(results_dir, "Tables", f_out),
                sep = "\t",row.names = FALSE)

    pathways <- as.character(rownames(path_up))
    path_up <- cbind(pathways, path_up)
    f_out <- paste(analysis_name, "PathwayUp.txt", sep="_")
    write.table(path_up, file = file.path(results_dir, "Tables", f_out),
                sep = "\t", row.names = FALSE)

    pathways <- as.character(rownames(path_down))
    path_down <- cbind(pathways, path_down)
    f_out <- paste(analysis_name, "PathwayDown.txt", sep="_")
    write.table(path_down, file = file.path(results_dir, "Tables", f_out),
                sep = "\t", row.names = FALSE)

    pathways <- as.character(rownames(path_any))
    path_any  <- cbind(pathways, path_any)
    f_out <- paste(analysis_name, "PathwayAny.txt", sep="_")
    write.table(path_any, file = file.path(results_dir, "Tables", f_out),
                sep = "\t", row.names = FALSE)
}














