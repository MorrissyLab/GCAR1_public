library(ActivePathways)
library(GSEABase)
library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(viridis)
library(circlize)
library(tidyr)
library(gprofiler2)
library(purrr)
library(tibble)
library(lsa)
library(RColorBrewer)
library(GSVA)
library(Seurat)
library(msigdbr)
library(stringr)
library(parallel)
library(gridExtra)
library(fastcluster)
library(cowplot)
library(ggrepel)
library(parallel)
library(data.table)
library(future)
library(future.apply)
library(doParallel)
library(rstatix)

setwd("~/Documents/GitHub/GCAR1_public")
SAVE_FOLDER <- paste(getwd(), "/02_Visium_HD/02_analysis/tmp/", sep = "")
PROGRAM_ANNOTATION <- NULL
RANK <- 15

# MGS Parameters and Switches ----
SPECTRA <- "02_Visium_HD/04_datasets/02_cnmf_outputs/cNMF_5_30_5_ST.gene_spectra_score.k_15.dt_2_0.txt"
GMT_FILE <- "02_Visium_HD/04_datasets/04_marker_gene_score/Supplementary_Table_7.gmt"
GMT_FILE_ALIAS <- "Supplementary_Table_7_GMT"

GENE_RANK_PLOTTING <- T
COMPUTE_MGS <- T
TOP_GENES_SPECTRA <- T
gPROFILER_RESULTS <- T

# Program Analysis Parameters and Switches ----
usage <- "02_Visium_HD/04_datasets/02_cnmf_outputs/cNMF_5_30_5_ST.usages.k_15.dt_2_0.consensus.txt"
usage <- read.table(usage, header = TRUE, sep = "")
colnames(usage) <- paste("GEP", seq(1, length(usage)), sep="")
usage_norm <- t(apply(usage, 1, function (x) x/sum(x)))
usage_norm <- usage_norm[complete.cases(usage_norm),]
program_annotation <- seq(1, RANK)
program_annotation <- ifelse(program_annotation <= 9, paste("P", program_annotation, sep = "0"), paste("P", program_annotation, sep = ""))
program_annotations <- c("Fibroblasts", "T_PCD_insulin", "T_lipidbio", "T_hyp_glycolysis",
                         "T_lipidox_Wnt", "VSMC_pericytes", "Endothelial_cells", "TAM_IFN",
                         "CD8", "CD4", "T_autophagy", "TAM_lipid", "TAM_APP",
                         "Lung_cells", "Muscle")
program_annotations <- paste(program_annotation, program_annotations, sep = "_")
program_annotations <- program_annotations

PROGRAM_ENRICHMENT <- T
PROGRAM_GENE_EXPRESSION <- T

# Neighbourhood Analysis Parameters and Switches ----
SAMPLE_COLS <- c("pNA" = "#9B5DE5", "bS1" = "#00F5D4", "bS2" = "#0068E8")
# SAMPLE_COLS <- c("pNA" = "#9B5DE5", "bS1" = "#00F5D4", "bS2" = "#0068E8", "xeno" = "#d5c58a")

NEIGHBORHOOD_ANALYSIS_PROGRAM_PAIRS <- T
NEIGHBORHOOD_ANALYSIS_SUMMARY <- T
NEIGHBORHOOD_ANALYSIS_SUMMARY_PROGRAMS <- c("9_10")

# Gene expression analysis parameters ----
GENE_EXPRESSION <- T
GENE_EXPRESSION_LIST <- c("GCAR", "GPNMB", "CD274", "PDCD1")
GENE_EXPRESSION_NORMALIZE <- list("1" = c("GPNMB"),
                                  "2" = c("CD3D", "CD3E", "CD3G", "CD4", "CD8A"))
# The following can be either "any" or "all". It corresponds to normalize to the number of spots
# with "any" expression of the genes or with expression of "all" the genes 
GENE_EXPRESSION_NORMALIZE_METHOD <- "all"

INSIDE_OUTSIDE_T_CELL <- T
INSIDE_OUTSIDE_NICHE <- T

INSIDE_OUTSIDE_GENES <- data.frame(category = c(rep("Pair_1", 2), 
                                                rep("Pair_2", 3), 
                                                rep("Pair_3", 3),
                                                rep("Pair_4", 2),
                                                rep("Pair_5", 2)),
                                   gene = c("PDCD1", "CD274", 
                                            "CTLA4", "CD80", "CD86",
                                            "TIGIT", "PVR", "NECTIN2",
                                            "KLRB1", "CLEC2D",
                                            "FAS", "FASLG"))
INSIDE_OUTSIDE_PROGRAMS_LIST <- c("9", "10")

# Save global variable names ----
program_var_names <- c(ls(), "program_var_names")

# General Functions ----
create_folder <- function(path){
  if (file.exists(path)){
  } else {
    dir.create(file.path(path), recursive = TRUE)
  }
}

gmt_to_list <- function(gmt_obj){
  geneset_list <- list()
  
  for (pathway in gmt_obj){
    geneset_name <- pathway$id
    geneset_genes <- pathway$genes
    
    if (geneset_name %in% names(geneset_list)){
      geneset_name_repeat_suffix <- length(which(names(geneset_list) == geneset_name)) + 1
      geneset_name <- paste(geneset_name, geneset_name_repeat_suffix, sep = ".")
    }
    
    geneset_list[[geneset_name]] <- geneset_genes
  }
  
  return(geneset_list)
}

# MGS Functions ----
load_spectra <- function(spectra_file){
  gene_score <- read.table(spectra_file, sep="\t", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
  gene_score <- t(as.matrix(gene_score))
}

ranked_spectra <- function(spectra_df, num_genes = "all", program_annotation = NULL){
  gene_score <- spectra_df
  
  # Set the max genes to include for each GEP
  if (num_genes == "all"){
    num_genes <- length(gene_score[,1])
  }
  
  # Add the top genes from each program from the selected rank to the GEP_df
  GEP_df <- data.frame(V1 = 1:num_genes)
  for (program in 1:ncol(gene_score)){
    sorted_spectra_df <- data.frame(genes = rownames(gene_score), value = gene_score[, program])
    sorted_spectra_df <- sorted_spectra_df[order(match(sorted_spectra_df$value, sort(sorted_spectra_df$value, decreasing = TRUE))),]
    
    top_genes <- sorted_spectra_df[1:num_genes, 1]
    GEP_df[program] <- top_genes
  }
  
  # Rename the columns for the GEP_df
  if (is.null(program_annotation)){
    program_annotation <- seq(1, ncol(gene_score))
    program_annotation <- ifelse(program_annotation <= 9, paste("GEP", program_annotation, sep = "_0"), paste("GEP", program_annotation, sep = "_"))
  }
  colnames(GEP_df) <- program_annotation
  
  return(GEP_df)
} 

marker_gene_score <- function(gep_df, gene_set, scale = FALSE, detailed = FALSE, filter = FALSE,
                              species = "Human", rank = 30){
  gene_score <- gep_df
  rownames(gene_score) <- paste("Usage_", seq(1, nrow(gene_score)), sep="")
  gene_score <- t(gene_score)
  gene_score <- data.frame(gene_score)
  
  final_df <- data.frame()
  
  if (detailed == FALSE){
    # For each pathway, get usage of the pathway genes in the GEP.
    if (!detailed){
      for (pathway in 1:length(geneset)){
        pathway_name <- geneset[[pathway]]$id
        marker_genes <- geneset[[pathway_name]]$genes
        
        # The human gene IDs are all upper case, but mouse ones are in title case.
        if (species != "Human"){
          marker_genes <- tolower(marker_genes)
          marker_genes <- tools::toTitleCase(marker_genes)
        }
        marker_genes <- gsub(" ", "", unique(marker_genes))
        
        pathway_gene_df <- data.frame(cell_type = pathway_name, n_markers = length(marker_genes))
        sub1 <- data.frame(Index = 1)
        
        for (GEP in 1:nrow(gene_score)){
          int_N <- length(intersect(gep_df[, GEP], marker_genes))
          ranks <- 1/(which(gep_df[, GEP] %in% intersect(gep_df[, GEP], marker_genes)))
          ranks_sum <- sum(ranks)
          
          marker_score <- int_N*ranks_sum
          marker_score <- marker_score/log10(length(marker_genes))
          
          sub2 <- data.frame(topN_marker_score = marker_score)
          colnames(sub2) <- paste("K", rank, "_GEP", GEP, "_", colnames(sub2), sep="")
          colnames(sub2) <- gsub("topN", paste("top", dim(gep_df)[1], sep=""), colnames(sub2))
          
          sub1 <- cbind(sub1, sub2)
        }
        sub1 <- cbind(pathway_gene_df, sub1)
        final_df <- rbind(final_df, sub1)
      }
      
      rownames_names <- names(geneset)
      rownames_names <- make.unique(rownames_names)
      rownames(final_df) <- rownames_names
    }
  }
  else {
    for (pathway in 1:length(geneset)){
      dataset_name = geneset[[pathway]]$name
      pathway_name <- geneset[[pathway]]$id
      
      marker_genes <- geneset[[pathway_name]]$genes
      
      # The human gene IDs are all upper case, but mouse ones are in title case.
      if (species != "Human"){
        marker_genes <- tolower(marker_genes)
        marker_genes <- tools::toTitleCase(marker_genes)
      }
      marker_genes <- gsub(" ", "", unique(marker_genes))
      
      marker_genes <- setdiff(marker_genes, "")
      
      pathway_gene_df <- data.frame(dataset = dataset_name, cell_type = pathway_name, n_markers = length(marker_genes))
      sub1 <- data.frame(Index = 1)
      
      for (GEP in 1:nrow(gene_score)){
        int_N <- length(intersect(gep_df[, GEP], marker_genes))
        ranks <- 1/(which(gep_df[, GEP] %in% intersect(gep_df[, GEP], marker_genes)))
        ranks_sum <- sum(ranks)
        
        marker_score <- int_N*ranks_sum
        marker_score <- marker_score/log10(length(marker_genes))
        
        intersect_genes = paste(intersect(gep_df[, GEP], marker_genes ), collapse="|")
        intersect_genes_rank = paste(which(gep_df[, GEP] %in% intersect(gep_df[,GEP], marker_genes )), collapse="|")
        
        sub2 <- data.frame(n_Int = int_N, genes_Int = intersect_genes, ranks_Int = intersect_genes_rank, topN_marker_score = marker_score)
        colnames(sub2) <- paste("K", rank, "_GEP", GEP, "_", colnames(sub2), sep="")
        colnames(sub2) <- gsub("topN", paste("top", dim(gep_df)[1], sep=""), colnames(sub2))
        
        pathway_gene_df <- cbind(pathway_gene_df, sub2)
      }
      max_score <- max(pathway_gene_df[,grep("_marker_score", colnames(pathway_gene_df), value = TRUE)])
      pathway_gene_df <- pathway_gene_df %>% add_column(maxScore = max_score, ".after" = "n_markers")
      final_df <- rbind(final_df, pathway_gene_df)
    }
  }
  
  if (filter){
    # If there's a cell_type, that is not present in the GEPs, remove the row
    if(length(which(rowSums(final_df) == 0)) > 0){
      MARKER_SCORES2 <- final_df[-which(rowSums(final_df) == 0),]
    }
    else {
      MARKER_SCORES2 <- final_df
    }
  }
  
  if (scale){
    final_df <- scale(t(final_df))
  }
  
  return(final_df)
}

gene_position <- function(gene_by_position_df, gene_list){
  # Gets the rank of select genes across every column (program/cluster)
  gene_position_df <- data.frame(gene_name = NA, row = NA, col = NA)
  
  for (gene in gene_list){
    temp_locations <- which(gep_df == gene, arr.ind = TRUE)
    temp_locations <- as.data.frame(temp_locations)
    temp_locations$gene_name <- gene
    gene_position_df <- rbind(gene_position_df, temp_locations)
  }
  
  gene_position_df <- gene_position_df[-1,]
  
  gene_position_df$row_inverse <- (gene_position_df$row ** -1)
  gene_position_df$row_inverse_logged <- log(gene_position_df$row_inverse)
  
  return(gene_position_df)
}

find_element_indices <- function(df, element){
  sapply(df, function(column) {
    index <- which(column == element)
    if (length(index) == 0) NA else index
  })
}

# Plotting specific functions
gmt_heatmap <- function(final_df, type = 'cnmf', 
                        file_name = "", file_height = 10, file_width = 25){
  # Remove the text based columns, empty rows, and empty columns
  plot_df <- final_df[-c(1,3)]
  plot_df <- plot_df[complete.cases(plot_df),]
  plot_df <- plot_df[apply(plot_df, 1, sum) != 0,]
  plot_df <- plot_df[apply(plot_df, 2, sum) != 0]
  
  col_markers <- colorRamp2(breaks = range(plot_df$n_markers), hcl_palette = "Prgn")
  haR <- HeatmapAnnotation(GMT_nr_marker = plot_df$n_markers, which = 'row', 
                           col = list(GMT_nr_marker = col_markers))
  
  plot_df <- plot_df[-1]
  plot_df <- as.matrix(plot_df)
  
  col_matrix <- colorRamp2(breaks = range(0, 0.02, 0.04, 0.06, 0.08, 0.1), 
                           hcl_palette = "Blues", reverse = TRUE)
  
  if (type == "cnmf"){
    ht_clustered <- Heatmap(plot_df, right_annotation = haR, column_title_rot = 45,
                            cluster_columns = TRUE, cluster_rows = TRUE, col = col_matrix)
    
    ht_col_clustered <- Heatmap(plot_df, right_annotation = haR, column_title_rot = 45,
                                cluster_columns = TRUE, cluster_rows = FALSE, col = col_matrix)
    
    ht_row_clustered <- Heatmap(plot_df, right_annotation = haR, column_title_rot = 45,
                                cluster_columns = FALSE, cluster_rows = TRUE, col = col_matrix)
    
    ht_unclustered <- Heatmap(plot_df, right_annotation = haR, column_title_rot = 45,
                              cluster_columns = FALSE, cluster_rows = FALSE, col = col_matrix)
  }
  
  # Add an additional annotation on the bottom w/ the number of genes per cluster.
  # This is only meant for a gene list from FindMarkers from seurat clusters.
  else{
    haB_genes <- colSums(!is.na(spectra_df))
    col_genes <- colorRamp2(breaks = range(haB_genes), hcl_palette = "Prgn")
    haB <- HeatmapAnnotation(nr_genes = haB_genes, which = 'column',
                             col = list(nr_genes = col_genes))
    
    ht_clustered <- Heatmap(plot_df, right_annotation = haR, column_title_rot = 45, bottom_annotation = haB,
                            cluster_columns = TRUE, cluster_rows = TRUE, col = col_matrix)
    
    ht_col_clustered <- Heatmap(plot_df, right_annotation = haR, column_title_rot = 45, bottom_annotation = haB,
                                cluster_columns = TRUE, cluster_rows = FALSE, col = col_matrix)
    
    ht_row_clustered <- Heatmap(plot_df, right_annotation = haR, column_title_rot = 45, bottom_annotation = haB,
                                cluster_columns = FALSE, cluster_rows = TRUE, col = col_matrix)
    
    ht_unclustered <- Heatmap(plot_df, right_annotation = haR, column_title_rot = 45, bottom_annotation = haB,
                              cluster_columns = FALSE, cluster_rows = FALSE, col = col_matrix)
  }
  
  pdf(file_name, height = file_height, width = file_width)
  draw(ht_clustered, padding = unit(c(3, 3, 3, 3), "cm"))
  draw(ht_col_clustered, padding = unit(c(3, 3, 3, 3), "cm"))
  draw(ht_row_clustered, padding = unit(c(3, 3, 3, 3), "cm"))
  draw(ht_unclustered, padding = unit(c(3, 3, 3, 3), "cm"))
  dev.off()
}

gene_position_linegraph <- function(gep_gene_position_df, 
                                    plot_title = "Temp Plot Title (K5)", 
                                    plot_x_title = "Gene Expression Program",
                                    file_name = "temp.pdf", 
                                    file_height = 10, file_width = 25){
  
  gep_gene_position <- gep_gene_position_df
  
  gg_r <- ggplot(gep_gene_position, aes(x=col, y=row)) +
    geom_line() + scale_x_continuous(breaks=seq(0, max(gep_gene_position$col), 4)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + theme_minimal() + 
    labs(title = plot_title, x = plot_x_title, y = "Gene Rank") +
    facet_wrap(. ~ gene_name, nrow = length(unique(gep_gene_position$gene_name)))
  gg_ri <- ggplot(gep_gene_position, aes(x=col, y=row_inverse)) +
    geom_line() + scale_x_continuous(breaks=seq(0, max(gep_gene_position$col), 4)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + theme_minimal() + 
    labs(title = plot_title, x = plot_x_title, y = "Inverse Gene Rank") +
    facet_wrap(. ~ gene_name, nrow = length(unique(gep_gene_position$gene_name)))
  gg_ril <- ggplot(gep_gene_position, aes(x=col, y=row_inverse_logged)) +
    geom_line() + scale_x_continuous(breaks=seq(0, max(gep_gene_position$col), 4)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + theme_minimal() + 
    labs(title = plot_title, x = plot_x_title, y = "Natural Log Inverse of Gene Rank") +
    facet_wrap(. ~ gene_name, nrow = length(unique(gep_gene_position$gene_name)))
  
  pdf(file_name, height = file_height, width = file_width)
  plot(gg_r)
  plot(gg_ri)
  plot(gg_ril)
  dev.off()
}

gene_position_heatmap <- function(gep_gene_position_df, heatmap_col = "row_inverse",
                                  file_name = "", file_height = 10, file_width = 25){
  gep_gene_position_wide <- gep_gene_position_df
  gep_gene_position_wide <- gep_gene_position_wide %>%
    pivot_wider(names_from = col, values_from = heatmap_col, id_cols = gene_name)
  gep_gene_position_wide <- as.data.frame(gep_gene_position_wide)
  row.names(gep_gene_position_wide) <- gep_gene_position_wide$gene_name
  gep_gene_position_wide <- gep_gene_position_wide[,-1]
  gep_gene_position_wide <- as.matrix(gep_gene_position_wide)
  
  scaled_mat <- t(scale(t(gep_gene_position_wide)))
  
  ht_clustered <- Heatmap(scaled_mat, cluster_columns = TRUE, cluster_rows = TRUE)
  ht_col_clustered <- Heatmap(scaled_mat, cluster_columns = TRUE, cluster_rows = FALSE)
  ht_row_clustered <- Heatmap(scaled_mat, cluster_columns = FALSE, cluster_rows = TRUE)
  ht_unclustered <- Heatmap(scaled_mat, cluster_columns = FALSE, cluster_rows = FALSE)
  
  pdf(file_name, height = file_height, width = file_width)
  draw(ht_clustered, padding = unit(c(3, 3, 3, 3), "cm"))
  draw(ht_col_clustered, padding = unit(c(3, 3, 3, 3), "cm"))
  draw(ht_row_clustered, padding = unit(c(3, 3, 3, 3), "cm"))
  draw(ht_unclustered, padding = unit(c(3, 3, 3, 3), "cm"))
  dev.off()
}

gProfiler_results <- function(GEP_df, source_filter = TRUE, term_filter = TRUE){
  gep_df_list <- as.list(GEP_df)
  
  # Run gProfiler
  gprofiler_output <- gost(gep_df_list, organism = "mmusculus", ordered_query = TRUE, multi_query = TRUE)
  
  # Convert the pVals from lists to their own dataframe with columns as programs
  gprofiler_pVal <- t(as.data.frame(gprofiler_output$result$p_values))
  gprofiler_pVal <- as.data.frame(gprofiler_pVal)
  
  # Set the column names to the program name inputs, -log10 of the pVals
  colnames(gprofiler_pVal) <- names(gep_df_list)
  gprofiler_pVal <- -log10(gprofiler_pVal)
  gprofiler_pVal <- as.data.frame(gprofiler_pVal)
  
  # Merge the pathway information with the pVals
  gprofiler_pVal_info <- data.frame(source = gprofiler_output$result$source, 
                                    term_name = gprofiler_output$result$term_name, 
                                    term_size = gprofiler_output$result$term_size)
  gprofiler_merged_pVal_info <- cbind(gprofiler_pVal_info, gprofiler_pVal)
  
  # Filter the output data frame to only GO:BP
  if (source_filter){
    gprofiler_merged_pVal_info <- gprofiler_merged_pVal_info %>% filter(source == "GO:BP")
  }
  # Filter the output data frame to only pathways with > 10 and < 2000 terms
  if (term_filter){
    gprofiler_merged_pVal_info <- gprofiler_merged_pVal_info %>% filter(term_size > 10) %>% filter(term_size < 2000)
  }
  
  # Set the rownames to be the term names and make sure they're unique
  rownames(gprofiler_merged_pVal_info) <- make.unique(gprofiler_merged_pVal_info$term_name)
  gprofiler_merged_pVal_info$term_name <- NULL
  
  # For each term, identify which program has the greatest significance for it
  gProf_res_sub2 <- gprofiler_merged_pVal_info
  gProf_res_sub2 <- data.frame(gProf_res_sub2) %>% add_column("group" = "", .before = colnames(gProf_res_sub2)[1])
  for(i in 1:nrow(gProf_res_sub2)){
    gProf_res_sub2$group[i] <- names(which.max(gProf_res_sub2[i,][4:ncol(gProf_res_sub2)]))
  }
  
  # Create a vector/df to reorder the terms based on the most significant term
  groups <- unique(gProf_res_sub2$group)
  Terms_order <- data.frame() 
  for(i in 1:length(groups)){
    gProf_res_sub3 <- gProf_res_sub2 %>% filter(group == groups[i]) 
    gProf_res_sub3 <- gProf_res_sub3[which(colnames(gProf_res_sub3) %in% c("group", groups[i]))]
    gProf_res_sub3 <- gProf_res_sub3[order(gProf_res_sub3[,2], decreasing = TRUE),]
    Terms_order <- rbind(Terms_order, data.frame(V1 = rownames(gProf_res_sub3), V2 = groups[i]))
  }
  
  # Reorder the terms in the dataframe using the vector/df
  gprofiler_merged_pVal_info_reordered <- gprofiler_merged_pVal_info[Terms_order$V1,]
  
  # Remove the term information
  gprofiler_merged_pVal_info_reordered <- as.data.frame(gprofiler_merged_pVal_info_reordered)
  gprofiler_merged_pVal_info_reordered_data <- gprofiler_merged_pVal_info_reordered[,-c(1,2)]
  gprofiler_merged_pVal_info_reordered_data <- as.matrix(gprofiler_merged_pVal_info_reordered_data)
  
  # Assign certain variables to the GlobalEnv, so they can be accessed for other functions later on
  assign("gProfiler_output", gprofiler_output, envir = .GlobalEnv)
  assign("gProfiler_merged_pVal_info", gprofiler_merged_pVal_info, envir = .GlobalEnv)
  assign("gProfiler_merged_pVal_info_reordered", gprofiler_merged_pVal_info_reordered, envir = .GlobalEnv)
  
  return(gprofiler_merged_pVal_info_reordered_data)
}

gProfiler_results_plots <- function(gProfiler_merged_pVal_info_reordered_df, pathway = "all", save_path){
  ht_opt$message = FALSE
  
  if (pathway != "all"){
    pathway_sub <- gsub(":", "", pathway)
    save_folder <- paste(save_path, pathway_sub, "/", sep = "")
    create_folder(save_folder)
  }
  else{
    pathway_sub <- "combinedPathways"
    save_folder <- save_path
  }
  
  base_file_name <- paste(save_folder, "K", RANK, "_gProfiler_", sep = "")
  output_file_name <- paste(base_file_name, pathway_sub, "_top1000.pdf", sep = "")
  capped_output_file_name <- paste(base_file_name, "capped_", pathway_sub, "_top1000.pdf", sep = "")
  capped_output_diff_file_name <- paste(base_file_name, "capped_diff_", pathway_sub, "_top1000.pdf", sep = "")

  # Retrieve the specific output ----
  if (pathway != "all"){
    output <- gProfiler_merged_pVal_info_reordered %>% filter(source == pathway) 
  }
  else{
    output <- gProfiler_merged_pVal_info_reordered
  }
  output <- output[,-c(1,2)]
  
  capped_output <- output
  capped_output[capped_output > 10] <- 10
  
  # Remove the pathways in the output that are the same across all the programs
  capped_output_diff <- as.data.frame(capped_output)
  capped_output_diff$total <- rowSums(capped_output_diff)
  capped_output_diff <- capped_output_diff[capped_output_diff$total != max(capped_output_diff),]
  capped_output_diff <- capped_output_diff[,-length(capped_output_diff)]
  
  # Convert all inputs to a matrix
  final_outputs <- list()
  final_outputs[[1]] <- as.matrix(output)
  final_outputs[[2]] <- as.matrix(capped_output)
  final_outputs[[3]] <- as.matrix(capped_output_diff)
  
  # Plotting ----
  final_outputs_save_folders <- c(output_file_name, capped_output_file_name, capped_output_diff_file_name)
  col_fun <- brewer.pal(name = "Blues", n = 9)
  
  for (output_index in seq(1:length(final_outputs))){
    temp <- final_outputs[[output_index]]
    temp_save_name <- final_outputs_save_folders[output_index]
    
    ht_clustered <- Heatmap(temp, show_row_names = TRUE, 
                            cluster_columns = TRUE, cluster_rows = TRUE,
                            col = col_fun, height = unit(0.5, "cm")*nrow(temp))
    ht_col_clustered <- Heatmap(temp, show_row_names = TRUE, 
                                cluster_columns = TRUE, cluster_rows = FALSE,
                                col = col_fun, height = unit(0.5, "cm")*nrow(temp))
    ht_row_clustered <- Heatmap(temp, show_row_names = TRUE, 
                                cluster_columns = FALSE, cluster_rows = TRUE,
                                col = col_fun, height = unit(0.5, "cm")*nrow(temp))
    ht_unclustered <- Heatmap(temp, show_row_names = TRUE, 
                              cluster_columns = FALSE, cluster_rows = FALSE,
                              col = col_fun, height = unit(0.5, "cm")*nrow(temp))
    pdf(temp_save_name, width = 25, height = (0.1969*nrow(capped_output)) + 1.3386)
    draw(ht_clustered, padding = unit(c(3, 3, 3, 3), "cm"))
    draw(ht_col_clustered, padding = unit(c(3, 3, 3, 3), "cm"))
    draw(ht_row_clustered, padding = unit(c(3, 3, 3, 3), "cm"))
    draw(ht_unclustered, padding = unit(c(3, 3, 3, 3), "cm"))
    dev.off()
  }
}

GEP_correlation_singlerank <- function(GS_dataset1, GS_dataset2){
  GS1 <- GS_dataset1
  GS2 <- GS_dataset2
  
  correlation_mat <- cor(GS1, GS2, method = "spearman")
  return(correlation_mat)
}

# Program Analysis Functions ----
program_enrichment <- function(usage_norm, threshold = 0.1, file_name, file_height = 6, file_width = 8, program_annotation = NULL){
  cutoff <- as.data.frame(usage_norm)
  cutoff$sample <- paste(sapply(strsplit(row.names(cutoff), "_"), "[", 1))
  
  # Get a df that contains the sample, # of total sample spots, # of the sample spots above the threshold per program (cpc), and a proportion.
  # The proportion is the cpc / # of total sample spots.
  p1_df <- cutoff %>%
    group_by(sample) %>%
    summarise(across(1:(ncol(cutoff) - 1), ~ sum(. > threshold), .names = "GEP{col}"),
              total_spots = n()) %>%
    pivot_longer(cols = starts_with("GEP"), names_to = "program", values_to = "cpc") %>%
    mutate(program = as.numeric(gsub("GEP", "", program)),
           proportion = cpc / total_spots) %>%
    ungroup()
  
  # Do the same as above but, this time, combine the 2 biopsy samples.
  p2_df <- cutoff
  p2_df$sample <- sub("bS1", "bS", p2_df$sample); p2_df$sample <- sub("bS2", "bS", p2_df$sample)
  p2_df <- p2_df %>%
    group_by(sample) %>%
    summarise(across(1:(ncol(p2_df) - 1), ~ sum(. > threshold), .names = "GEP{col}"),
              total_spots = n()) %>%
    pivot_longer(cols = starts_with("GEP"), names_to = "program", values_to = "cpc") %>%
    mutate(program = as.numeric(gsub("GEP", "", program)),
           proportion = cpc / total_spots) %>%
    ungroup()
  
  # If the programs have been annotated, then change the program names from GEP... to the annotated names.
  if (!is.null(program_annotation)){
    p1_df$program <- rep(program_annotations, length(unique(p1_df$sample)))
    p2_df$program <- rep(program_annotations, length(unique(p2_df$sample))) 
  } else {
    p1_df$program <- ifelse(p1_df$program <= 9, paste("GEP", p1_df$program, sep = "0"), paste("GEP", p1_df$program, sep = ""))
    p2_df$program <- ifelse(p2_df$program <= 9, paste("GEP", p2_df$program, sep = "0"), paste("GEP", p2_df$program, sep = ""))
  }
  
  # Change the sample names and make them factors, so they order themselves on the chart
  p1_df$sample[p1_df$sample == "pNA"] <- "Primary"
  p1_df$sample[p1_df$sample == "bS1"] <- "Biopsy #1"
  p1_df$sample[p1_df$sample == "bS2"] <- "Biopsy #2"
  p1_df$sample[p1_df$sample == "xeno"] <- "Xenograft"
  p2_df$sample[p2_df$sample == "pNA"] <- "Primary"
  p2_df$sample[p2_df$sample == "bS"] <- "Biopsy"
  p2_df$sample[p2_df$sample == "xeno"] <- "Xenograft"
  
  p1_df$sample <- factor(p1_df$sample, levels = c("Primary", "Biopsy #1", "Biopsy #2", "Xenograft"))
  p2_df$sample <- factor(p2_df$sample, levels = c("Primary", "Biopsy", "Xenograft"))
  
  p1_df <- as.data.frame(p1_df); p2_df <- as.data.frame(p2_df)
  colnames(p1_df)[1] <- "Sample"; colnames(p2_df)[1] <- "Sample"
  
  # Plot options
  cols <- c("#224760", "#f5c43e", "#d73b1f", "#d5c58a")
  cols_v2 <- cols[-3]
  plot_title <- paste("K", ncol(usage_norm), "Program Enrichment", sep = " ")
  
  pdf(file = file_name, height = file_height, width = file_width)
  p1 <- ggplot(p1_df, aes(x = program, y = cpc, fill = Sample)) + 
    geom_col(colour = "black", position = "fill") + theme_classic2() +
    labs(title = plot_title, x = "Gene Expression Program", y = "Proportion of Spots", color = "Sample",
         subtitle = paste("Threshold - ", threshold, sep = "")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = "black"), 
          axis.text.y = element_text(colour="black"),
          legend.position = "top", plot.margin = unit(c(1,10,1,1), "pt"), text = element_text(size = 20)) + 
    theme(plot.margin = unit(c(1,1,1,1), "cm")) + scale_y_continuous(expand = c(0,0)) + 
    scale_fill_manual(values = cols)
  p2 <- ggplot(p1_df, aes(x = program, y = proportion, fill = Sample)) + 
    geom_bar(stat="identity", color = "black", position = "dodge", width = .8) + theme_classic2() +
    labs(title = plot_title, x = "Gene Expression Program", y = "Proportion of Spots", color = "Sample",
         subtitle = paste("Threshold - ", threshold, sep = "")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = "black"), 
          axis.text.y = element_text(colour="black"),
          legend.position = "top", plot.margin = unit(c(1,10,1,1), "pt"), text = element_text(size = 20)) + 
    theme(plot.margin = unit(c(1,1,1,1), "cm")) + scale_y_continuous(expand = c(0,0)) + 
    scale_fill_manual(values = cols)
  p3 <- ggplot(p2_df, aes(x = program, y = proportion, fill = Sample)) + 
    geom_bar(stat="identity", color = "black", position = "dodge", width = .8) + theme_classic2() +
    labs(title = plot_title, x = "Gene Expression Program", y = "Proportion of Spots", color = "Sample",
         subtitle = paste("Threshold - ", threshold, sep = "")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = "black"), 
          axis.text.y = element_text(colour="black"),
          legend.position = "top", plot.margin = unit(c(1,10,1,1), "pt"), text = element_text(size = 20)) + 
    theme(plot.margin = unit(c(1,1,1,1), "cm")) + scale_y_continuous(expand = c(0,0)) + 
    scale_fill_manual(values = cols_v2)
  print(p1)
  print(p2)
  print(p3)
  dev.off()
}

# Neighbourhood Analysis Functions ----
seurat_object_coords <- function(){
  # Get spot centroids from the 3 samples
  pNA_spots <- as.data.frame(primary_no_addon@images$slice1.024um@boundaries$centroids@coords)
  bS1_spots <- as.data.frame(biopsy_s1@images$slice1.024um@boundaries$centroids@coords)
  bS2_spots <- as.data.frame(biopsy_s2@images$slice1.024um@boundaries$centroids@coords)
  xeno_spots <- as.data.frame(xenograft_no_addon@images$slice1.024um@boundaries$centroids@coords)
  
  # Set the rownames to the cell IDs
  rownames(pNA_spots) <- primary_no_addon@images$slice1.024um@boundaries$centroids@cells
  rownames(bS1_spots) <- biopsy_s1@images$slice1.024um@boundaries$centroids@cells
  rownames(bS2_spots) <- biopsy_s2@images$slice1.024um@boundaries$centroids@cells
  rownames(xeno_spots) <- xenograft_no_addon@images$slice1.024um@boundaries$centroids@cells
  
  # Add a prefix to the rownames, so they match what's in the merged seurat object
  rownames(pNA_spots) <- paste("pNA", rownames(pNA_spots), sep = "_")
  rownames(bS1_spots) <- paste("bS1", rownames(bS1_spots), sep = "_")
  rownames(bS2_spots) <- paste("bS2", rownames(bS2_spots), sep = "_")
  rownames(xeno_spots) <- paste("xeno", rownames(xeno_spots), sep = "_")
  
  # Take the union of the spots, rename the columns, and make a column to specify the sample
  union_spots <- do.call("rbind", list(pNA_spots, bS1_spots, bS2_spots, xeno_spots))
  colnames(union_spots) <- c("row", "col")
  union_spots$sample <- paste(sapply(strsplit(rownames(union_spots), "_"), "[", 1))
  
  return(union_spots)
}

program_neighbours <- function(program_spots, cell_centroids, max_radius = 3){
  program_spots <- program_spots
  radius <- 75
  distances <- seq(0, max_radius) * radius
  
  # Convert cell_centroids to a matrix for faster indexing
  cell_centroids_matrix <- as.matrix(cell_centroids[, c("col", "row")])
  cell_centroid_samples <- cell_centroids$sample  # Store sample column separately for indexing
  
  # Pre-compute bounds for each radius level to avoid repetitive calculations
  bounds <- lapply(distances, function(k) list(lower = -k, upper = k))
  
  # Parallelize across `program_spots` if possible
  cl <- makeCluster(parallelly::availableCores() - 1)  # Use available cores minus one
  clusterExport(cl, c("program_spots", "cell_centroids_matrix", "cell_centroid_samples", 
                      "bounds", "distances", "radius", "max_radius"),
                envir=environment())
  clusterEvalQ(cl, library(data.table))
  
  final_results <- parLapply(cl, seq_along(program_spots), function(i) {
    spot <- cell_centroids_matrix[program_spots[i], ]
    sample <- sub("_.*", "", rownames(cell_centroids)[i])
    
    union_spots_subset <- cell_centroids_matrix[cell_centroid_samples == sample, , drop = FALSE]
    subset_indices <- which(cell_centroid_samples == sample)
    
    # Store intermediate results for each radius
    spot_results <- vector("list", length(bounds))
    
    for (j in seq_along(bounds)) {
      k_bounds <- bounds[[j]]
      
      # Apply bounds to subset neighbors
      mask_x <- (union_spots_subset[, "col"] >= (spot["col"] + k_bounds$lower)) &
        (union_spots_subset[, "col"] <= (spot["col"] + k_bounds$upper))
      mask_y <- (union_spots_subset[, "row"] >= (spot["row"] + k_bounds$lower)) &
        (union_spots_subset[, "row"] <= (spot["row"] + k_bounds$upper))
      
      neighbors <- subset_indices[mask_x & mask_y]
      
      # Record result
      spot_results[[j]] <- data.frame(
        sample = sample,
        spot_of_interest = program_spots[i],
        r = distances[j],
        type = "circular",
        neighbours = I(list(rownames(cell_centroids)[neighbors])),
        nr_neighbours = length(neighbors)
      )
    }
    rbindlist(spot_results)
  })
  
  stopCluster(cl)  # Stop the parallel cluster
  
  # Combine all parallelized results
  final_neigh_df <- rbindlist(final_results)
  
  return(final_neigh_df)
}

remove_overlapping_neighbours <- function(spot_neighbours_df){
  # Convert to data.table for efficient updates
  spot_neighbours_dt <- as.data.table(spot_neighbours_df)
  
  # Ensure radii are unique and sorted
  radii <- sort(unique(as.numeric(spot_neighbours_dt$r)))
  
  # Set up parallel cluster backend
  num_cores <- parallelly::availableCores() - 1  # Leave one core free
  plan(cluster, workers = num_cores)
  
  # Split unique spots into chunks for parallel processing
  unique_spots <- unique(spot_neighbours_dt$spot_of_interest)
  spot_chunks <- split(unique_spots, cut(seq_along(unique_spots), num_cores, labels = FALSE))
  
  # Define function to process each chunk
  process_chunk <- function(spots_chunk) {
    # Initialize a list to store neighbors encountered for each spot of interest
    neighbors_encountered <- vector("list", length(spots_chunk))
    names(neighbors_encountered) <- spots_chunk
    
    # Process each spot of interest in the chunk
    for (soi in spots_chunk) {
      # Initialize encountered neighbors
      neighbors_encountered[[soi]] <- character(0)
      
      # Loop over each radius
      for (rad in radii) {
        # Logical indexing to select rows by spot of interest and radius
        soi_rows <- spot_neighbours_dt[spot_of_interest == soi & r == rad]
        
        # Apply unique and setdiff to the neighbors
        current_neighbors <- unique(unlist(soi_rows$neighbours))
        unique_neighbors <- setdiff(current_neighbors, neighbors_encountered[[soi]])
        
        # Update encountered neighbors list
        neighbors_encountered[[soi]] <- c(neighbors_encountered[[soi]], unique_neighbors)
        
        # Update the data.table with unique neighbors and count
        spot_neighbours_dt[spot_of_interest == soi & r == rad, `:=`(
          neighbours = list(unique_neighbors),
          nr_neighbours = length(unique_neighbors)
        )]
        
        # Debug: Check if unique_neighbors is being correctly assigned
        if (length(unique_neighbors) == 0) {
          message("Empty unique_neighbors for spot_of_interest = ", soi, ", radius = ", rad)
        }
      }
    }
    return(spot_neighbours_dt[spot_of_interest %in% spots_chunk, ])
  }
  
  # Apply processing in parallel using future_lapply with cluster strategy
  results <- future_lapply(spot_chunks, process_chunk)
  
  # Combine results into a single data.table
  final_result <- rbindlist(results)
  
  # Return final result as a data.frame if needed
  return(as.data.frame(final_result))
}

# MGS: Figure 5d-g; Extended Figure 7b. ----
# 1.1 Make the T-cell and PD1/PDL1 gene rank graph ----
if (GENE_RANK_PLOTTING){
  # 1.0 Clear existing non-function variables, set the save folders (create them if they doesn't exist) ----
  rm(list = setdiff(setdiff(ls(), lsf.str()), program_var_names))
  save_path <- paste(SAVE_FOLDER, "annotation/K", RANK, "/genes_ranked/", sep = "")
  create_folder(save_path)
  create_folder(paste(save_path, "T_Cell_Subunits/", sep = ""))
  create_folder(paste(save_path, "PD1_PDL1/", sep = ""))
  create_folder(paste(save_path, "resistance_genes/", sep = ""))
  
  # 1.1 Calculate all the gene ranks per program from the spectra ----
  spectra_df <- load_spectra(SPECTRA)
  gep_df <- ranked_spectra(spectra_df)
  
  indices <- data.frame("PD1" = find_element_indices(gep_df, "PDCD1"),
                        "PDL1" = find_element_indices(gep_df, "CD274"),
                        "CTLA4" = find_element_indices(gep_df, "CTLA4"),
                        "CD80" = find_element_indices(gep_df, "CD80"),
                        "CD86" = find_element_indices(gep_df, "CD86"),
                        "TIGIT" = find_element_indices(gep_df, "TIGIT"),
                        "PVR" = find_element_indices(gep_df, "PVR"),
                        "NECTIN2" = find_element_indices(gep_df, "NECTIN2"),
                        "KLRB1" = find_element_indices(gep_df, "KLRB1"),
                        "CLEC2D" = find_element_indices(gep_df, "CLEC2D")) 
  
  # 1.2 Retrieve specific gene positions ----
  gene_position_t_cell <- gene_position(gep_df, c("CD3D", "CD3E", "CD3G", "CD8A", "CD4"))
  gene_position_resistance <- gene_position(gep_df, c("PDCD1", "CD274"))
  gene_position_resistance_extened <- gene_position(gep_df, c("PDCD1", "CD274", "CTLA4", "CD80", "CD86", "TIGIT", "PVR", "NECTIN2", "KLRB1", "CLEC2D"))
  
  # 1.3 Save the ranked genes per gep data frame ----
  write.table(gep_df, paste(save_path, "K", RANK, "_ranked_genes.txt", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)
  
  # 1.4 Plot the gene ranks ----
  gene_position_linegraph(gene_position_t_cell, plot_title = paste("T-Cell Subunits (K", RANK, ")", sep = ""), 
                          plot_x_title = "Gene Expression Program",
                          file_name = paste(save_path, "T_Cell_Subunits/K", RANK, "_CD3_subunits_spike_graph.pdf", sep = ""), 
                          file_height = 10, file_width = 25)
  gene_position_heatmap(gene_position_t_cell, heatmap_col = "row_inverse",
                        file_name = paste(save_path, "T_Cell_Subunits/K", RANK, "_CD3_subunits_heatmap.pdf", sep = ""), 
                        file_height = 10, file_width = 25)
  
  # For CD274 and PDCD1
  gene_position_linegraph(gene_position_resistance, plot_title = paste("K", RANK, " - PDCD1 and CD274 Spike Graph", sep = ""), 
                          plot_x_title = "Gene Expression Program",
                          file_name = paste(save_path, "PD1_PDL1/K", RANK, "_PDCD1_CD274_spike_graph.pdf", sep = ""), 
                          file_height = 10, file_width = 25)
  gene_position_heatmap(gene_position_resistance, heatmap_col = "row_inverse",
                        file_name = paste(save_path, "PD1_PDL1/K", RANK, "_PDCD1_CD274_heatmap.pdf", sep = ""),
                        file_height = 10, file_width = 25)
  
  # For all the resistance genes
  write.table(indices, file = paste(save_path, "/resistance_genes/resistance_gene_ranks.tsv", sep = ""), quote = FALSE, sep = "\t")
  gene_position_linegraph(gene_position_resistance_extened, plot_title = paste("K", RANK, " - Resistance Genes Spike Graph", sep = ""), 
                          plot_x_title = "Gene Expression Program",
                          file_name = paste(save_path, "resistance_genes/K", RANK, "_resistance_genes_spike_graph.pdf", sep = ""), 
                          file_height = 10, file_width = 25)
  gene_position_heatmap(gene_position_resistance_extened, heatmap_col = "row_inverse",
                        file_name = paste(save_path, "resistance_genes/K", RANK, "_resistance_genes_heatmap.pdf", sep = ""),
                        file_height = 10, file_width = 25)
}

# 1.2 Compute MGS for the top 1000 genes ----
if (COMPUTE_MGS){
  # 2.0 Clear existing non-function variables, set the save folders (create them if they doesn't exist) ----
  rm(list = setdiff(setdiff(ls(), lsf.str()), program_var_names))
  save_path <- paste(SAVE_FOLDER, "annotation/K", RANK, "/MarkerGeneScore/", sep = "")
  create_folder(save_path)
  create_folder(paste(save_path, GMT_FILE_ALIAS, sep = ""))
  
  # 2.1 Calculate all the gene ranks per program from the spectra ----
  spectra_df <- load_spectra(SPECTRA)
  gep_df <- ranked_spectra(spectra_df, num_genes = 1000)
  
  # 2.2 Load the GMT file ---- 
  geneset <- read.GMT(GMT_FILE)
  
  # 2.3 Compute the regular and detailed marker gene scores ----
  final_df_simple <- marker_gene_score(gep_df, geneset, scale = FALSE, 
                                       detailed = FALSE, filter = FALSE,
                                       species = "Human", rank = RANK)
  final_df_detailed <- marker_gene_score(gep_df, geneset, scale = FALSE, 
                                detailed = TRUE, filter = FALSE,
                                species = "Human", rank = RANK)
  
  # 2.4 Save the output tables ----
  write.table(final_df_simple, paste(save_path, GMT_FILE_ALIAS, "/K", RANK, "_MGS_top1000_genes.txt", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(final_df_detailed, paste(save_path, GMT_FILE_ALIAS, "/K", RANK, "_MGS_top1000_genes_detailed.txt", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)

  # Plot a heatmap of the GMT gene scores. This should be ran with the top 1000 genes.
  gmt_heatmap(final_df_simple, 
              file_name = paste(save_path, GMT_FILE_ALIAS, "/K", RANK, "_MGS_top1000_genes.pdf", sep = ""),
              file_height = 25, file_width = 25)
}

# 1.3 Plot the top 10 genes per program based on their spectra score ----
if (TOP_GENES_SPECTRA){
  # 3.0 Clear existing non-function variables, set the save folders (create them if they doesn't exist) ----
  rm(list = setdiff(setdiff(ls(), lsf.str()), program_var_names))
  save_path <- paste(SAVE_FOLDER, "annotation/K", RANK, "/genes_ranked/", sep = "")
  create_folder(save_path)
  
  # 3.1 Calculate all the gene ranks per program from the spectra ----
  spectra_df <- load_spectra(SPECTRA)
  gep_df <- ranked_spectra(spectra_df, num_genes = 1000, program_annotation = PROGRAM_ANNOTATION)
  
  # 3.2 Get the top 10 genes per program ranked based on the spectra ---- 
  top_genes <- gep_df[c(1:10),]
  top_genes <- top_genes %>%
    pivot_longer(everything(), names_to = "gep", values_to = "gene") %>%
    mutate(program = as.numeric(gsub("GEP_", "", gep))) %>%
    arrange(program)
  top_gene_scores <- t(spectra_df[c(top_genes$gene),]); top_gene_scores <- scale(top_gene_scores)
  
  top_gene_scores_annot <- unique(top_genes$program)
  haT <- HeatmapAnnotation(program_of_top_scoring_gene = top_genes$gep,
                           which = "column")
  haR <- HeatmapAnnotation(program = top_gene_scores_annot,
                           col = list(program = colorRamp2(range(top_gene_scores_annot), 
                                                           hcl_palette = "Blue-Red", reverse = TRUE)),
                           which = "row")
  
  col_fun <- colorRamp2(breaks = range(top_gene_scores), hcl_palette = "RdYlBu", reverse = TRUE)
  ht <- Heatmap(top_gene_scores, column_title_rot = 45, 
                right_annotation = haR, top_annotation = haT,
                cluster_columns = FALSE, cluster_rows = FALSE, 
                col = col_fun, name = "Scaled Spectra Scores for \nTop 10 Scoring Genes Per Program")
  ht_row_clust <- Heatmap(top_gene_scores, column_title_rot = 45, 
                          right_annotation = haR, top_annotation = haT,
                          cluster_columns = FALSE, cluster_rows = TRUE, 
                          col = col_fun, name = "Scaled Spectra Scores for \nTop 10 Scoring Genes Per Program")
  ht_col_clust <- Heatmap(top_gene_scores, column_title_rot = 45, 
                          right_annotation = haR, top_annotation = haT,
                          cluster_columns = TRUE, cluster_rows = FALSE, 
                          col = col_fun, name = "Scaled Spectra Scores for \nTop 10 Scoring Genes Per Program")
  ht_clust <- Heatmap(top_gene_scores, column_title_rot = 45, 
                      right_annotation = haR, top_annotation = haT,
                      cluster_columns = TRUE, cluster_rows = TRUE, 
                      col = col_fun, name = "Scaled Spectra Scores for \nTop 10 Scoring Genes Per Program")
  pdf(paste(save_path, "K", RANK, "_Top_10_Scoring_Genes.pdf", sep = ""), width = 50, height = 10)
  draw(ht, padding = unit(c(3, 3, 3, 3), "cm"))
  draw(ht_row_clust, padding = unit(c(3, 3, 3, 3), "cm"))
  draw(ht_col_clust, padding = unit(c(3, 3, 3, 3), "cm"))
  draw(ht_clust, padding = unit(c(3, 3, 3, 3), "cm"))
  dev.off() 
}

# 1.4 Generate gProfiler results ----
if (gPROFILER_RESULTS){
  # 4.0 Clear existing non-function variables, set the save folders (create them if they doesn't exist) ----
  rm(list = setdiff(setdiff(ls(), lsf.str()), program_var_names))
  save_path <- paste(SAVE_FOLDER, "annotation/K", RANK, "/gProfiler/", sep = "")
  create_folder(save_path)
  
  # 4.1 Calculate all the gene ranks per program from the spectra ----
  spectra_df <- load_spectra(SPECTRA)
  gep_df <- ranked_spectra(spectra_df, num_genes = 10, program_annotation = PROGRAM_ANNOTATION)
  
  # 4.2 Run gProfiler on the gep_df ----
  output <- gProfiler_results(gep_df, source_filter = FALSE)
  
  # 4.3 Create a gost plot of the results ----
  pdf(paste(save_path, "K", RANK, "_gostplot_top1000.pdf", sep = ""), height = 50, width = 7)
  gostplot(gProfiler_output, capped = TRUE, interactive = FALSE)
  dev.off()
 
  # 4.4 Create the pathways plot of the results ----
  pathways <- unique(gProfiler_merged_pVal_info_reordered$source)
  pathways <- c("all", pathways)
  for (pathway in pathways){
    print(pathway)
    gProfiler_results_plots(gProfiler_merged_pVal_info_reordered, 
                            pathway = pathway,
                            save_path) 
  }
  
  # 4.5 Make the selected ranks annotation file ----
  sra_file <- data.frame(Rank = RANK, GEP = seq(1:RANK), 
                         Label = "", Annotation_notes_SM = "", 
                         Annotation = "", 
                         Comment = "", GO_BP = "")
  
  sra_file$GO_BP <- apply(output, 2, function(col) {
    top_500 <- order(col, decreasing = TRUE)[1:500]
    top_500_pathways <- rownames(output)[top_500]
    paste(top_500_pathways, collapse = " | ")
  })
  write.table(sra_file, file = paste(SAVE_FOLDER, "annotation/K", RANK, "/SelectedRankAnnotation.txt", sep = ""), 
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# Program Analysis: Figure 5d-g ----
# Load the Seurat Objects ----
rm(list = setdiff(setdiff(ls(), lsf.str()), program_var_names))

# Classical loading
# primary_no_addon <- Load10X_Spatial(data.dir = primary_no_addon_path, bin.size = c(24))
# biopsy_s1 <- Load10X_Spatial(data.dir = biopsy_s1_path, bin.size = c(24))
# biopsy_s2 <- Load10X_Spatial(data.dir = biopsy_s2_path, bin.size = c(24))
# xenograft_no_addon <- Load10X_Spatial(data.dir = xenograft_no_addon_path, bin.size = c(24))

# Reviewer loading
primary_no_addon <- readRDS("02_Visium_HD/04_datasets/01_seurat_objects/primary_sample_seurat_object.rds")
biopsy_s1 <- readRDS("02_Visium_HD/04_datasets/01_seurat_objects/biopsy_sample_1_seurat_object.rds")
biopsy_s2 <- readRDS("02_Visium_HD/04_datasets/01_seurat_objects/biopsy_sample_2_seurat_object.rds")
xenograft_no_addon <- readRDS("02_Visium_HD/04_datasets/01_seurat_objects/xenograft_sample_seurat_object.rds")

# combined <- merge(primary_no_addon, y=c(biopsy_s1, biopsy_s2),
#                   add.cell.ids = c("pNA", "bS1", "bS2"),
#                   project = "combined_samples")
combined <- merge(primary_no_addon, y=c(biopsy_s1, biopsy_s2, xenograft_no_addon),
                  add.cell.ids = c("pNA", "bS1", "bS2", "xeno"),
                  project = "combined_samples")
combined_layers <- JoinLayers(combined)
combined_layers <- subset(combined_layers, nCount_Spatial.024um > 0)

program_var_names <- c(ls(), "program_var_names")

# 2.1 Program Enrichment  ----
if (PROGRAM_ENRICHMENT){
  # 1.0 Clear existing non-function variables, set the save folders (create them if they doesn't exist) ----
  rm(list = setdiff(setdiff(ls(), lsf.str()), program_var_names))
  save_path <- paste(SAVE_FOLDER, "K", RANK, "/program_enrichment/", sep = "")
  create_folder(save_path)
  
  # 1.1 Basic program enrichment across samples ----
  program_enrichment(usage_norm, threshold = 0.05, 
                     file_name = paste(save_path, "/K", RANK, "_Program_Enrichment.pdf", sep = ""), 
                     file_height = 10, file_width = 12,
                     program_annotation = program_annotation)
  
}

# 2.2 Panel Gene Expression  ----
if (PROGRAM_GENE_EXPRESSION){
  # 2.0 Clear existing non-function variables, set the save folders (create them if they doesn't exist) ----
  rm(list = setdiff(setdiff(ls(), lsf.str()), program_var_names))
  save_path <- paste(SAVE_FOLDER, "K", RANK, "/program_gene_expression/", sep = "")
  create_folder(save_path)
  
  # 2.1 Make a plot showcasing gene expression across programs for resistance genes ----
  binarized_usage <- ifelse(usage_norm < 0.05, 0, 1)
  
  resistance_gene_panel <- data.frame(category = c("ASPS",
                                                   rep("Pair_1", 2), 
                                                   rep("Pair_2", 3), 
                                                   rep("Pair_3", 3),
                                                   rep("Pair_4", 2),
                                                   rep("Pair_5", 2)),
                                 gene = c("GPNMB", 
                                          "PDCD1", "CD274", 
                                          "CTLA4", "CD80", "CD86",
                                          "TIGIT", "PVR", "NECTIN2",
                                          "KLRB1", "CLEC2D",
                                          "FAS", "FASLG"))
  gene_panel <- resistance_gene_panel
  combined_count_table <- as.data.frame(t(as.matrix(combined_layers@assays$Spatial.024um$counts[unlist(gene_panel$gene) ,])))
  
  binarized_usage_spots <- rownames(binarized_usage[rowSums(binarized_usage) > 0, ])
  combined_count_table_subset <- combined_count_table[, colnames(combined_count_table) %in% gene_panel$gene]
  combined_count_table_subset <- combined_count_table_subset[rownames(combined_count_table_subset) %in% binarized_usage_spots,]
  
  # 2.1.1 Make a plot that is tissue agnostic ----
  corr_matrix <- cor(combined_count_table_subset, usage_norm)
  corr_matrix_cosine <- cosine(as.matrix(cbind(combined_count_table_subset, usage_norm)))
  corr_matrix_cosine <- corr_matrix_cosine[c(1:(nrow(corr_matrix))), c((nrow(corr_matrix)+1):ncol(corr_matrix_cosine))]
  
  gene_order <- gene_panel$gene[gene_panel$gene %in% rownames(corr_matrix)]
  corr_matrix <- corr_matrix[match(gene_order, rownames(corr_matrix)),]
  corr_matrix_cosine <- corr_matrix_cosine[match(gene_order, rownames(corr_matrix_cosine)),]
  
  split <- factor(gene_panel$category[gene_panel$gene %in% rownames(corr_matrix)],
                  levels = unique(gene_panel$category[gene_panel$gene %in% rownames(corr_matrix)]))
  
  ht <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                cluster_columns = FALSE, cluster_rows = FALSE,
                row_split = split, row_title_rot = 0)
  ht_row <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                    cluster_columns = FALSE, cluster_rows = TRUE,
                    row_split = split, row_title_rot = 0)
  ht_col <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                    cluster_columns = TRUE, cluster_rows = FALSE,
                    row_split = split, row_title_rot = 0)
  ht_both <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                     cluster_columns = TRUE, cluster_rows = TRUE,
                     row_split = split, row_title_rot = 0)
  ht_cosine <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                       cluster_columns = FALSE, cluster_rows = FALSE,
                       row_split = split, row_title_rot = 0)
  ht_cosine_row <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                           cluster_columns = FALSE, cluster_rows = TRUE,
                           row_split = split, row_title_rot = 0)
  ht_cosine_col <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                           cluster_columns = TRUE, cluster_rows = FALSE,
                           row_split = split, row_title_rot = 0)
  ht_cosine_both <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                            cluster_columns = TRUE, cluster_rows = TRUE,
                            row_split = split, row_title_rot = 0)
  
  p1 <- grid.grabExpr(draw(ht, padding = unit(c(3, 3, 3, 3), "cm"))); p2 <- grid.grabExpr(draw(ht_cosine, padding = unit(c(3, 3, 3, 3), "cm")))
  p3 <- grid.grabExpr(draw(ht_row, padding = unit(c(3, 3, 3, 3), "cm"))); p4 <- grid.grabExpr(draw(ht_cosine_row, padding = unit(c(3, 3, 3, 3), "cm")))
  p5 <- grid.grabExpr(draw(ht_col, padding = unit(c(3, 3, 3, 3), "cm"))); p6 <- grid.grabExpr(draw(ht_cosine_col, padding = unit(c(3, 3, 3, 3), "cm")))
  p7 <- grid.grabExpr(draw(ht_both, padding = unit(c(3, 3, 3, 3), "cm"))); p8 <- grid.grabExpr(draw(ht_cosine_both, padding = unit(c(3, 3, 3, 3), "cm")))
  
  pdf(paste(save_path, "K", RANK, "_Resistance_Panel_gene_expression_program_usage_correlation.pdf", sep = ""), width = 40, height = 10)
  print(plot_grid(p1, p2))
  print(plot_grid(p3, p4))
  print(plot_grid(p5, p6))
  print(plot_grid(p7, p8))
  dev.off()
  
  # 2.1.2 Make a plot that is tissue specific ----
  tissues_list <- unique(paste(sapply(strsplit(row.names(usage_norm), "_"), "[", 1)))
  tissues_list <- unique(gsub("bS1|bS2", "bS", tissues_list))
  
  for (tissue in tissues_list){
    combined_count_table_subset_tissue <- combined_count_table_subset[grep(tissue, rownames(combined_count_table_subset)), ]
    usage_norm_subset <- usage_norm[grep(tissue, rownames(usage_norm)), ]
    
    corr_matrix <- cor(combined_count_table_subset_tissue, usage_norm_subset)
    corr_matrix_cosine <- cosine(as.matrix(cbind(combined_count_table_subset_tissue, usage_norm_subset)))
    
    corr_matrix_cosine <- corr_matrix_cosine[c(1:(nrow(corr_matrix))), c((nrow(corr_matrix)+1):ncol(corr_matrix_cosine))]
    
    gene_order <- gene_panel$gene[gene_panel$gene %in% rownames(corr_matrix)]
    corr_matrix <- corr_matrix[match(gene_order, rownames(corr_matrix)),]
    
    corr_matrix_cosine <- corr_matrix_cosine[match(gene_order, rownames(corr_matrix_cosine)),]
    
    split <- factor(gene_panel$category[gene_panel$gene %in% rownames(corr_matrix)],
                    levels = unique(gene_panel$category[gene_panel$gene %in% rownames(corr_matrix)]))
    
    ht <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                  cluster_columns = FALSE, cluster_rows = FALSE,
                  row_split = split, row_title_rot = 0)
    ht_row <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                      cluster_columns = FALSE, cluster_rows = TRUE,
                      row_split = split, row_title_rot = 0)
    ht_col <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                      cluster_columns = TRUE, cluster_rows = FALSE,
                      row_split = split, row_title_rot = 0)
    ht_both <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                       cluster_columns = TRUE, cluster_rows = TRUE,
                       row_split = split, row_title_rot = 0)
    ht_cosine <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                         cluster_columns = FALSE, cluster_rows = FALSE,
                         row_split = split, row_title_rot = 0)
    ht_cosine_row <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                             cluster_columns = FALSE, cluster_rows = TRUE,
                             row_split = split, row_title_rot = 0)
    ht_cosine_col <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                             cluster_columns = TRUE, cluster_rows = FALSE,
                             row_split = split, row_title_rot = 0)
    ht_cosine_both <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                              cluster_columns = TRUE, cluster_rows = TRUE,
                              row_split = split, row_title_rot = 0)
    
    p1 <- grid.grabExpr(draw(ht, padding = unit(c(3, 3, 3, 3), "cm"))); p2 <- grid.grabExpr(draw(ht_cosine, padding = unit(c(3, 3, 3, 3), "cm")))
    p3 <- grid.grabExpr(draw(ht_row, padding = unit(c(3, 3, 3, 3), "cm"))); p4 <- grid.grabExpr(draw(ht_cosine_row, padding = unit(c(3, 3, 3, 3), "cm")))
    p5 <- grid.grabExpr(draw(ht_col, padding = unit(c(3, 3, 3, 3), "cm"))); p6 <- grid.grabExpr(draw(ht_cosine_col, padding = unit(c(3, 3, 3, 3), "cm")))
    p7 <- grid.grabExpr(draw(ht_both, padding = unit(c(3, 3, 3, 3), "cm"))); p8 <- grid.grabExpr(draw(ht_cosine_both, padding = unit(c(3, 3, 3, 3), "cm")))
    
    pdf(paste(save_path, "K", RANK, "_", tissue, "_Resistance_Panel_gene_expression_program_usage_correlation.pdf", sep = ""), width = 40, height = 10)
    print(plot_grid(p1, p2))
    print(plot_grid(p3, p4))
    print(plot_grid(p5, p6))
    print(plot_grid(p7, p8))
    dev.off()
  }
  
  # 2.2 For resistance w/o GPNMB ----
  resistance_gene_panel <- data.frame(category = c(rep("Pair_1", 2), 
                                                   rep("Pair_2", 3), 
                                                   rep("Pair_3", 3),
                                                   rep("Pair_4", 2),
                                                   rep("Pair_5", 2)),
                                      gene = c("PDCD1", "CD274", 
                                               "CTLA4", "CD80", "CD86",
                                               "TIGIT", "PVR", "NECTIN2",
                                               "KLRB1", "CLEC2D",
                                               "FAS", "FASLG"))
  
  gene_panel <- resistance_gene_panel
  combined_count_table <- as.data.frame(t(as.matrix(combined_layers@assays$Spatial.024um$counts[unlist(gene_panel$gene) ,])))
  
  binarized_usage_spots <- rownames(binarized_usage[rowSums(binarized_usage) > 0, ])
  combined_count_table_subset <- combined_count_table[, colnames(combined_count_table) %in% gene_panel$gene]
  combined_count_table_subset <- combined_count_table_subset[rownames(combined_count_table_subset) %in% binarized_usage_spots,]
  
  # 2.2.1 Make a plot that is tissue agnostic ----
  corr_matrix <- cor(combined_count_table_subset, usage_norm)
  corr_matrix_cosine <- cosine(as.matrix(cbind(combined_count_table_subset, usage_norm)))
  corr_matrix_cosine <- corr_matrix_cosine[c(1:(nrow(corr_matrix))), c((nrow(corr_matrix)+1):ncol(corr_matrix_cosine))]
  
  gene_order <- gene_panel$gene[gene_panel$gene %in% rownames(corr_matrix)]
  corr_matrix <- corr_matrix[match(gene_order, rownames(corr_matrix)),]
  corr_matrix_cosine <- corr_matrix_cosine[match(gene_order, rownames(corr_matrix_cosine)),]
  
  split <- factor(gene_panel$category[gene_panel$gene %in% rownames(corr_matrix)],
                  levels = unique(gene_panel$category[gene_panel$gene %in% rownames(corr_matrix)]))
  
  ht <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                cluster_columns = FALSE, cluster_rows = FALSE,
                row_split = split, row_title_rot = 0)
  ht_row <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                    cluster_columns = FALSE, cluster_rows = TRUE,
                    row_split = split, row_title_rot = 0)
  ht_col <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                    cluster_columns = TRUE, cluster_rows = FALSE,
                    row_split = split, row_title_rot = 0)
  ht_both <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                     cluster_columns = TRUE, cluster_rows = TRUE,
                     row_split = split, row_title_rot = 0)
  ht_cosine <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                       cluster_columns = FALSE, cluster_rows = FALSE,
                       row_split = split, row_title_rot = 0)
  ht_cosine_row <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                           cluster_columns = FALSE, cluster_rows = TRUE,
                           row_split = split, row_title_rot = 0)
  ht_cosine_col <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                           cluster_columns = TRUE, cluster_rows = FALSE,
                           row_split = split, row_title_rot = 0)
  ht_cosine_both <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                            cluster_columns = TRUE, cluster_rows = TRUE,
                            row_split = split, row_title_rot = 0)
  
  p1 <- grid.grabExpr(draw(ht, padding = unit(c(3, 3, 3, 3), "cm"))); p2 <- grid.grabExpr(draw(ht_cosine, padding = unit(c(3, 3, 3, 3), "cm")))
  p3 <- grid.grabExpr(draw(ht_row, padding = unit(c(3, 3, 3, 3), "cm"))); p4 <- grid.grabExpr(draw(ht_cosine_row, padding = unit(c(3, 3, 3, 3), "cm")))
  p5 <- grid.grabExpr(draw(ht_col, padding = unit(c(3, 3, 3, 3), "cm"))); p6 <- grid.grabExpr(draw(ht_cosine_col, padding = unit(c(3, 3, 3, 3), "cm")))
  p7 <- grid.grabExpr(draw(ht_both, padding = unit(c(3, 3, 3, 3), "cm"))); p8 <- grid.grabExpr(draw(ht_cosine_both, padding = unit(c(3, 3, 3, 3), "cm")))
  
  pdf(paste(save_path, "K", RANK, "_resistance_gene_expression_program_usage_correlation.pdf", sep = ""), width = 40, height = 10)
  print(plot_grid(p1, p2))
  print(plot_grid(p3, p4))
  print(plot_grid(p5, p6))
  print(plot_grid(p7, p8))
  dev.off()
  
  # 2.2.2 Make a plot that is tissue specific ----
  tissues_list <- unique(paste(sapply(strsplit(row.names(usage_norm), "_"), "[", 1)))
  tissues_list <- unique(gsub("bS1|bS2", "bS", tissues_list))
  
  for (tissue in tissues_list){
    combined_count_table_subset_tissue <- combined_count_table_subset[grep(tissue, rownames(combined_count_table_subset)), ]
    usage_norm_subset <- usage_norm[grep(tissue, rownames(usage_norm)), ]
    
    corr_matrix <- cor(combined_count_table_subset_tissue, usage_norm_subset)
    corr_matrix_cosine <- cosine(as.matrix(cbind(combined_count_table_subset_tissue, usage_norm_subset)))
    
    corr_matrix_cosine <- corr_matrix_cosine[c(1:(nrow(corr_matrix))), c((nrow(corr_matrix)+1):ncol(corr_matrix_cosine))]
    
    gene_order <- gene_panel$gene[gene_panel$gene %in% rownames(corr_matrix)]
    corr_matrix <- corr_matrix[match(gene_order, rownames(corr_matrix)),]
    
    corr_matrix_cosine <- corr_matrix_cosine[match(gene_order, rownames(corr_matrix_cosine)),]
    
    split <- factor(gene_panel$category[gene_panel$gene %in% rownames(corr_matrix)],
                    levels = unique(gene_panel$category[gene_panel$gene %in% rownames(corr_matrix)]))
    
    ht <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                  cluster_columns = FALSE, cluster_rows = FALSE,
                  row_split = split, row_title_rot = 0)
    ht_row <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                      cluster_columns = FALSE, cluster_rows = TRUE,
                      row_split = split, row_title_rot = 0)
    ht_col <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                      cluster_columns = TRUE, cluster_rows = FALSE,
                      row_split = split, row_title_rot = 0)
    ht_both <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                       cluster_columns = TRUE, cluster_rows = TRUE,
                       row_split = split, row_title_rot = 0)
    ht_cosine <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                         cluster_columns = FALSE, cluster_rows = FALSE,
                         row_split = split, row_title_rot = 0)
    ht_cosine_row <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                             cluster_columns = FALSE, cluster_rows = TRUE,
                             row_split = split, row_title_rot = 0)
    ht_cosine_col <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                             cluster_columns = TRUE, cluster_rows = FALSE,
                             row_split = split, row_title_rot = 0)
    ht_cosine_both <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                              cluster_columns = TRUE, cluster_rows = TRUE,
                              row_split = split, row_title_rot = 0)
    
    p1 <- grid.grabExpr(draw(ht, padding = unit(c(3, 3, 3, 3), "cm"))); p2 <- grid.grabExpr(draw(ht_cosine, padding = unit(c(3, 3, 3, 3), "cm")))
    p3 <- grid.grabExpr(draw(ht_row, padding = unit(c(3, 3, 3, 3), "cm"))); p4 <- grid.grabExpr(draw(ht_cosine_row, padding = unit(c(3, 3, 3, 3), "cm")))
    p5 <- grid.grabExpr(draw(ht_col, padding = unit(c(3, 3, 3, 3), "cm"))); p6 <- grid.grabExpr(draw(ht_cosine_col, padding = unit(c(3, 3, 3, 3), "cm")))
    p7 <- grid.grabExpr(draw(ht_both, padding = unit(c(3, 3, 3, 3), "cm"))); p8 <- grid.grabExpr(draw(ht_cosine_both, padding = unit(c(3, 3, 3, 3), "cm")))
    
    pdf(paste(save_path, "K", RANK, "_", tissue, "_resistance_gene_expression_program_usage_correlation.pdf", sep = ""), width = 40, height = 10)
    print(plot_grid(p1, p2))
    print(plot_grid(p3, p4))
    print(plot_grid(p5, p6))
    print(plot_grid(p7, p8))
    dev.off()
  }
}

# Neighbourhood Analysis: Figure h; Extended Figure 7e. ----
# Create or just load the spot_neighbours df for all neighborhood analysis
rm(list = setdiff(setdiff(ls(), lsf.str()), program_var_names))

union_cell_coordinate <- seurat_object_coords()
# union_cell_coordinate <- union_cell_coordinate[rownames(union_cell_coordinate) %in% rownames(usage_norm),]
combined_metadata_table <- as.data.frame(combined_layers@meta.data)
options(future.globals.maxSize= 891289600)

# Check if the spot neighbours file exists, if not, then make it and save it 
if(!file.exists("02_Visium_HD/04_datasets/05_bin_neighbours/spot_neighbours_df.rds")){
  spot_neighbours_df <- program_neighbours(rownames(union_cell_coordinate), union_cell_coordinate)
  saveRDS(spot_neighbours_df, file = paste(SAVE_FOLDER, "spot_neighbours_df.rds", sep = ""))
} else{spot_neighbours_df <- readRDS("02_Visium_HD/04_datasets/05_bin_neighbours/spot_neighbours_df.rds")}
if(!file.exists("02_Visium_HD/04_datasets/05_bin_neighbours/spot_neighbours_df_wo_overlaps.rds")){
  spot_neighbours_df_wo_overlaps <- remove_overlapping_neighbours(spot_neighbours_df)
  saveRDS(spot_neighbours_df_wo_overlaps, file = paste(SAVE_FOLDER, "spot_neighbours_df_wo_overlaps.rds", sep = ""))
} else{spot_neighbours_df_wo_overlaps <- readRDS("02_Visium_HD/04_datasets/05_bin_neighbours/spot_neighbours_df_wo_overlaps.rds")}

program_var_names <- c(ls(), "program_var_names")

# 3.1 Neighbourhood of each pair of programs ----
if (NEIGHBORHOOD_ANALYSIS_PROGRAM_PAIRS){
  rm(list = setdiff(setdiff(ls(), lsf.str()), program_var_names))
  usage_threshold <- 0.05
  Binarized_threshold <- 0.05
  
  binarized_gep <- usage_norm
  binarized_gep[binarized_gep >= Binarized_threshold] <- 1; binarized_gep[binarized_gep < Binarized_threshold] <- 0
  
  # Parallelized ----
  # Subset to just program 9 and 10 for ease of running
  loaded_packages <- setdiff(loadedNamespaces(), c("base", "stats", "utils", "graphics", "grDevices", "methods", "datasets"))
  combinations <- combn(seq_len(ncol(usage_norm)), 2)
  program_combinations <- data.frame(program_one = combinations[1, ], program_two = combinations[2, ])
  program_combinations <- program_combinations[85,]
  
  # Set up the parallel backend with the number of available cores
  num_cores <- 1  # leave one core free
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  foreach(combination = seq_len(nrow(program_combinations)), .packages = loaded_packages) %dopar% {
    program_1 <- program_combinations[combination, 1]
    program_2 <- program_combinations[combination, 2]
    
    gep_pos_p1_name <- paste("GEP", program_1, sep = "_")
    gep_pos_p2_name <- paste("GEP", program_2, sep = "_")
    gep_post_p1_p2_name <- paste("GEP", program_1, program_2, sep = "_")
    print(gep_post_p1_p2_name)
    
    # 1.0 Clear existing non-function variables, set the save folders (create them if they doesn't exist) ----
    save_path <- paste(SAVE_FOLDER, "K", RANK, "/neighbourhood_analysis/Programs/", gep_post_p1_p2_name, "/", sep = "")
    create_folder(save_path)
    
    # 1.1 Get the GEP+ cells, the spot coordinates from the 3 samples, and their intersection ----
    gep_pos_spots_p1 <- rownames(usage_norm[usage_norm[, program_1] > usage_threshold,])
    gep_pos_spots_p2 <- rownames(usage_norm[usage_norm[, program_2] > usage_threshold,])
    
    # Get the spot names and merge them with the union_cell_coordinate dataframe
    # This gets the unique spot IDs. So, it only returns a single spot ID even if a spot belongs to both GEPs. 
    # This is intended and we will be differentiate what GEPs each spot belongs to in the heatmap using the annotation tracks.
    program_pair_spots <- rownames(usage_norm[usage_norm[, program_1] > usage_threshold | usage_norm[, program_2] > usage_threshold,])
    gep_pos_spots <- intersect(program_pair_spots, rownames(union_cell_coordinate))
    
    # 1.2 Subset the spot neighbours ----
    spot_neighbours_df_subset <- spot_neighbours_df_wo_overlaps[spot_neighbours_df_wo_overlaps$spot_of_interest %in% gep_pos_spots,]
    
    # 1.4 Make the spot_neighbour_programs_df ----
    # The final_df will have additional columns that are named on the unique program-radii combinations
    # Cbind the program_neighbors_df with an empty matrix w/ ncol = unique programs * unique radii
    program_radius_comb_df <- expand.grid(colnames(usage_norm), unique(spot_neighbours_df_subset$r))
    program_radius_comb_df$GEP_r_comb <- paste(program_radius_comb_df$Var1, program_radius_comb_df$Var2, sep = "_")
    temp_matrix <- matrix(0, 
                          nrow = 1, 
                          ncol = RANK*length(unique(spot_neighbours_df_subset$r)), 
                          dimnames = list(NULL, program_radius_comb_df$GEP_r_comb))
    spot_neighbour_programs_df <- cbind(spot_neighbours_df_subset, temp_matrix)
    
    # 1.5 Get the binarized usage sum per program for all a spots neighbours ----
    # Iterate i over each row (spot-radius combination) of the program_neighbours_df
    for (i in seq_len(nrow(spot_neighbour_programs_df))){
      # Get the neighbors for the current row (spot-radius combination)
      neighbors <- unlist(spot_neighbour_programs_df$neighbours[i])
      
      # Get the binary gene expression programs for the neighbors
      neighbor_programs <- binarized_gep[rownames(binarized_gep) %in% neighbors, , drop = FALSE]
      
      # Compute the column sums for the neighbors
      neighbor_programs_sum <- colSums(neighbor_programs) / nrow(neighbor_programs)
      
      # Update spot_neighbour_programs_df with the calculated sums using column name matching
      colnames_to_paste <- paste(colnames(neighbor_programs), spot_neighbour_programs_df$r[i], sep = "_")
      spot_neighbour_programs_df[i, c(colnames_to_paste)] <- neighbor_programs_sum
    }
    
    # 1.6 Remove the spot-radii combs, thereby flattening the df ----
    spot_neighbour_programs_df_flat <- data.frame()
    for (soi in unique(spot_neighbour_programs_df$spot_of_interest)){
      soi_program_subset <- spot_neighbour_programs_df[spot_neighbour_programs_df$spot_of_interest == soi,]
      
      gep_sums <- soi_program_subset[, c(7:ncol(soi_program_subset))]
      gep_sums <- colSums(gep_sums)
      
      row_data <- c(soi_program_subset[1, c(1, 2, 4)], c(soi_program_subset$nr_neighbours), gep_sums)
      row_data <- as.data.frame(row_data)
      colnames(row_data)[4:7] <- c("nr_neigh_r0", "nr_neigh_r1", "nr_neigh_r2", "nr_neigh_r3")
      
      if (length(spot_neighbour_programs_df_flat) == 0){
        spot_neighbour_programs_df_flat <- row_data
      }
      else{
        spot_neighbour_programs_df_flat <- rbind(spot_neighbour_programs_df_flat, row_data) 
      }
    }
    
    # 1.7 Modify the spot_neighbour_programs_df_flat for plotting ----
    # Add a tissue column and organize the columns
    spot_neighbour_programs_df_flat <- as.data.frame(spot_neighbour_programs_df_flat)
    spot_neighbour_programs_df_flat[, c(4:7)] <- as.numeric(unlist(spot_neighbour_programs_df_flat[, c(4:7)]))
    spot_neighbour_programs_df_flat$tissue <- spot_neighbour_programs_df_flat$sample
    spot_neighbour_programs_df_flat$tissue <- sub("bS1", "bS", spot_neighbour_programs_df_flat$tissue)
    spot_neighbour_programs_df_flat$tissue <- sub("bS2", "bS", spot_neighbour_programs_df_flat$tissue)
    spot_neighbour_programs_df_flat <- spot_neighbour_programs_df_flat[, c(1:7, ncol(spot_neighbour_programs_df_flat), 8:(ncol(spot_neighbour_programs_df_flat) - 1))]
    spot_neighbour_programs_df_flat[is.na(spot_neighbour_programs_df_flat)] <- 0
    
    # Make the annotation and plotting df
    spot_neighbour_programs_df_flat_annot <- spot_neighbour_programs_df_flat[,1:8]
    spot_neighbour_programs_df_flat_plot <- t(spot_neighbour_programs_df_flat[, c(9:ncol(spot_neighbour_programs_df_flat))])
    colnames(spot_neighbour_programs_df_flat_plot) <- c(spot_neighbour_programs_df_flat$spot_of_interest)
    
    # Set an annotation column for the program that the cell belongs to
    gep_pos_p1_annot <- ifelse(spot_neighbour_programs_df_flat_annot$spot_of_interest %in% gep_pos_spots_p1, "TRUE", "FALSE")
    gep_pos_p2_annot <- ifelse(spot_neighbour_programs_df_flat_annot$spot_of_interest %in% gep_pos_spots_p2, "TRUE", "FALSE")
    
    spot_neighbour_programs_df_flat_annot <- cbind(spot_neighbour_programs_df_flat_annot, gep_pos_p1_annot)
    spot_neighbour_programs_df_flat_annot <- cbind(spot_neighbour_programs_df_flat_annot, gep_pos_p2_annot)
    
    spot_neighbour_programs_df_flat_annot$GEP <- ifelse(spot_neighbour_programs_df_flat_annot$spot_of_interest %in% gep_pos_spots_p1, gep_pos_p1_name, gep_pos_p2_name)
    spot_neighbour_programs_df_flat_annot$GEP[spot_neighbour_programs_df_flat_annot$spot_of_interest %in% intersect(gep_pos_spots_p1, gep_pos_spots_p2)] <- gep_post_p1_p2_name
    
    # 1.8 Plotting ----
    # 1.8.1 Create the annotations ----
    GEP_col_names <- setNames(c("#891A53", "#1E7D42", "#5456C9"), c(gep_pos_p1_name, gep_pos_p2_name, gep_post_p1_p2_name))
    haT <- HeatmapAnnotation(sample = spot_neighbour_programs_df_flat_annot$sample,
                             GEP = spot_neighbour_programs_df_flat_annot$GEP,
                             nr_neigh_r0 = spot_neighbour_programs_df_flat_annot$nr_neigh_r0,
                             nr_neigh_r1 = spot_neighbour_programs_df_flat_annot$nr_neigh_r1,
                             nr_neigh_r2 = spot_neighbour_programs_df_flat_annot$nr_neigh_r2,
                             nr_neigh_r3 = spot_neighbour_programs_df_flat_annot$nr_neigh_r3,
                             col = list(sample = SAMPLE_COLS,
                                        GEP = GEP_col_names,
                                        nr_neigh_r0 = colorRamp2(c(range(spot_neighbour_programs_df_flat_annot$nr_neigh_r0), range(spot_neighbour_programs_df_flat_annot$nr_neigh_r0)+1), hcl_palette = "YlOrRd", reverse = TRUE),
                                        nr_neigh_r1 = colorRamp2(range(spot_neighbour_programs_df_flat_annot$nr_neigh_r1), hcl_palette = "YlOrRd", reverse = TRUE),
                                        nr_neigh_r2 = colorRamp2(range(spot_neighbour_programs_df_flat_annot$nr_neigh_r2), hcl_palette = "YlOrRd", reverse = TRUE),
                                        nr_neigh_r3 = colorRamp2(range(spot_neighbour_programs_df_flat_annot$nr_neigh_r3), hcl_palette = "YlOrRd", reverse = TRUE)),
                             which = "col")
    
    program_radius_comb_df$Var2 <- as.numeric(program_radius_comb_df$Var2)
    program_radius_comb_df$program <- gsub("GEP", "", program_radius_comb_df$Var1)
    program_radius_comb_df$program <- as.numeric(program_radius_comb_df$program)
    haL <- HeatmapAnnotation(radius = program_radius_comb_df$Var2, program = program_radius_comb_df$program,
                             col = list(radius = colorRamp2(range(program_radius_comb_df$Var2), 
                                                            hcl_palette = "Blue-Yellow", reverse = TRUE),
                                        program = colorRamp2(range(program_radius_comb_df$program), 
                                                             hcl_palette = "Inferno", reverse = TRUE)),
                             which = "row")
    
    # Set the variables that the plots will be split by
    sample_split <- names(SAMPLE_COLS)
    sample_split <- unique(gsub("bS1|bS2", "bS", sample_split))
    col_spliter <- factor(spot_neighbour_programs_df_flat_annot$tissue, levels = sample_split)
    col_spliter_GEP <- as.factor(spot_neighbour_programs_df_flat_annot$GEP)
    row_spliter <- factor(program_radius_comb_df$Var1)
    col_fun = colorRamp2(range(spot_neighbour_programs_df_flat_plot), hcl_palette = "Blues", reverse = TRUE)
    
    # Set the program annotation names
    rownames(spot_neighbour_programs_df_flat_plot) <- paste(program_annotations, c(rep(0, RANK), rep(75, RANK), rep(150, RANK), rep(225, RANK)), sep = "_")
    spot_neighbour_programs_df_flat_plot <- as.matrix(spot_neighbour_programs_df_flat_plot)
    
    # 1.8.3 Plot but w/o the gene expression tracks ----
    ht_combined <- Heatmap(spot_neighbour_programs_df_flat_plot, left_annotation = haL, top_annotation = haT,
                           name = "Mean of Neighbour Program",
                           cluster_columns = TRUE, cluster_rows = FALSE, col = col_fun, 
                           show_row_names = TRUE, show_column_names = FALSE,
                           column_split = col_spliter, row_split = row_spliter)
    ht_combined_GEP_split <- Heatmap(spot_neighbour_programs_df_flat_plot, left_annotation = haL, top_annotation = haT,
                                     name = "Mean of Neighbour Program",
                                     cluster_columns = TRUE, cluster_rows = FALSE, col = col_fun, 
                                     show_row_names = TRUE, show_column_names = FALSE,
                                     column_split = col_spliter_GEP, row_split = row_spliter)
    ht_combined_row_clust <- Heatmap(spot_neighbour_programs_df_flat_plot, left_annotation = haL, top_annotation = haT,
                                     name = "Mean of Neighbour Program",
                                     cluster_columns = TRUE, cluster_rows = TRUE, col = col_fun, 
                                     show_row_names = TRUE, show_column_names = FALSE,
                                     column_split = col_spliter)
    ht_combined_row_clust_GEP_split <- Heatmap(spot_neighbour_programs_df_flat_plot, left_annotation = haL, top_annotation = haT,
                                               name = "Mean of Neighbour Program",
                                               cluster_columns = TRUE, cluster_rows = TRUE, col = col_fun, 
                                               show_row_names = TRUE, show_column_names = FALSE,
                                               column_split = col_spliter_GEP)
    ht_combined_col_clust <- Heatmap(spot_neighbour_programs_df_flat_plot, left_annotation = haL, top_annotation = haT,
                                     name = "Mean of Neighbour Program",
                                     cluster_columns = TRUE, cluster_rows = FALSE, col = col_fun, 
                                     show_row_names = TRUE, show_column_names = FALSE,
                                     row_split = row_spliter)
    ht_combined_clust <- Heatmap(spot_neighbour_programs_df_flat_plot, left_annotation = haL, top_annotation = haT,
                                 name = "Mean of Neighbour Program",
                                 cluster_columns = TRUE, cluster_rows = TRUE, col = col_fun, 
                                 show_row_names = TRUE, show_column_names = FALSE)
    pdf(paste(save_path, gep_post_p1_p2_name, "_T", usage_threshold, "_Binarized_T", Binarized_threshold, "_Mean_Neighbour_Programs_No_Overlaps_No_GE_Tracks.pdf", sep = ""), width = 60, height = 30)
    draw(ht_combined, padding = unit(c(3, 3, 3, 3), "cm"))
    draw(ht_combined_GEP_split, padding = unit(c(3, 3, 3, 3), "cm"))
    draw(ht_combined_row_clust, padding = unit(c(3, 3, 3, 3), "cm"))
    draw(ht_combined_row_clust_GEP_split, padding = unit(c(3, 3, 3, 3), "cm"))
    draw(ht_combined_col_clust, padding = unit(c(3, 3, 3, 3), "cm"))
    draw(ht_combined_clust, padding = unit(c(3, 3, 3, 3), "cm"))
    dev.off()
    
    # 1.8.4 Save the Rdata image, so plotting changes can be easily performed ----
    save(list = ls(environment()), file = paste0(save_path, gep_post_p1_p2_name, "_neighbourhood_save.RData"))
  }
  stopCluster(cl)
}

# 3.2 Custom neighbourhood splitting and dot plot ----
if (NEIGHBORHOOD_ANALYSIS_SUMMARY){
  rm(list = setdiff(setdiff(ls(), lsf.str()), program_var_names))
  
  for (program_selection in NEIGHBORHOOD_ANALYSIS_SUMMARY_PROGRAMS){
    # 2.1 Clear existing non-function variables, set the save folders (create them if they doesn't exist) ----
    save_path <- paste(SAVE_FOLDER, "K", RANK, "/neighbourhood_analysis/Programs/GEP_", 
                       program_selection, "/GEP_", program_selection, sep = "")
    usage_threshold <- 0.05
    Binarized_threshold <- 0.05
    save_path_root <- paste(save_path, "_T", 
                            usage_threshold, "_Binarized_T", Binarized_threshold, sep = "")
    
    load(paste(save_path, "_neighbourhood_save.RData", sep = ""))
    
    # 2.2 Make the heatmap, extract Kmeans column order, and then further differentiate large niches ----
    ht_combined <- Heatmap(spot_neighbour_programs_df_flat_plot,
                           left_annotation = haL, top_annotation = haT,
                           name = "Mean of Neighbour Program",
                           cluster_columns = TRUE, cluster_rows = FALSE, col = col_fun,
                           show_row_names = TRUE, show_column_names = FALSE,
                           column_split = col_spliter, row_split = row_spliter, column_km = 2)
    ht_combined_drawn <- draw(ht_combined, padding = unit(c(3, 3, 3, 3), "cm"))
    heatmap_col_order <- column_order(ht_combined_drawn)
    heatmap_col_order <- heatmap_col_order[c(1,3,4)]
    
    temp <- spot_neighbour_programs_df_flat_plot[c(6,7), grepl("bS1|bS2", colnames(spot_neighbour_programs_df_flat_plot))]
    heatmap_col_order_names <- colnames(spot_neighbour_programs_df_flat_plot)[unlist(heatmap_col_order)]
    temp <- temp[,!colnames(temp) %in% unlist(heatmap_col_order_names)]
    heatmap_col_order[["2,bS_6_7_positive"]] <- unique(c(colnames(temp[,temp[1,] > 0]),
                                                         colnames(temp[,temp[2,] > 0])))
    heatmap_col_order[["3,bS_remainder"]] <- colnames(temp[, !colnames(temp) %in% unlist(heatmap_col_order)])
    heatmap_col_order[["2,bS_6_7_positive"]] <- which(colnames(spot_neighbour_programs_df_flat_plot) %in% heatmap_col_order[["2,bS_6_7_positive"]])
    heatmap_col_order[["3,bS_remainder"]] <- which(colnames(spot_neighbour_programs_df_flat_plot) %in% heatmap_col_order[["3,bS_remainder"]])
    
    heatmap_col_order_final <- unlist(lapply(seq_along(heatmap_col_order), function(i) {
      # Get the name and values of the current list element
      name <- names(heatmap_col_order)[i]
      values <- heatmap_col_order[[i]]
      
      # Return the name repeated for each value
      rep(name, length(values))
    }))[order(unlist(heatmap_col_order))]
    
    # 2.3 Plot the final heatmap with the more granular niches ----
    ht_combined <- Heatmap(spot_neighbour_programs_df_flat_plot, 
                           left_annotation = haL, top_annotation = haT,
                           name = "Mean of Neighbour Program",
                           cluster_columns = TRUE, cluster_rows = FALSE, col = col_fun,
                           show_row_names = TRUE, show_column_names = FALSE,
                           column_split = heatmap_col_order_final, row_split = row_spliter)
    pdf(paste(save_path_root, "_Mean_Neighbour_Programs_No_Overlaps_No_GE_Tracks_Custom_Split_Custom.pdf", sep = ""), width = 60, height = 30)
    draw(ht_combined, padding = unit(c(3, 3, 3, 3), "cm"))
    dev.off()
    
    # 2.4 Make a Dot plot with these niches ----
    # Subset to the data of interest
    combined_layers_subset <- combined_layers
    combined_layers_subset@meta.data$spot <- rownames(combined_layers_subset@meta.data)
    combined_layers_subset@meta.data$sample <- paste(sapply(strsplit(combined_layers_subset@meta.data$spot, "_"), "[", 1))
    combined_layers_subset@meta.data$sample <- gsub("bS1|bS2", "bS", combined_layers_subset@meta.data$sample)
    combined_layers_subset <- subset(combined_layers_subset, 
                                     subset = spot %in% spot_neighbour_programs_df_flat$spot_of_interest)
    
    # Because the order of the rows to col transpose is the same, we can just set the col indices as rows
    setdiff(colnames(spot_neighbour_programs_df_flat_plot), spot_neighbour_programs_df_flat$spot_of_interest)
    spot_neighbour_programs_df_flat_mod <- spot_neighbour_programs_df_flat
    spot_neighbour_programs_df_flat_mod$clust <- "0"
    
    for (split_name in names(heatmap_col_order)){
      split_name_reformatted <- sub("^([0-9]+),(.*)$", "\\2,\\1", split_name)
      split_name_reformatted <- sub(",", "_", split_name_reformatted)
      spot_neighbour_programs_df_flat_mod[heatmap_col_order[[split_name]], "clust"] <- split_name_reformatted
    }
    combined_layers_subset@meta.data <- cbind(combined_layers_subset@meta.data, 
                                              spot_neighbour_programs_df_flat_mod[9:ncol(spot_neighbour_programs_df_flat_mod)])
    cols_to_plot <- colnames(combined_layers_subset@meta.data)[6:(ncol(combined_layers_subset@meta.data)-1)]
    
    # 2.4.2 Create the dot plot ----
    pdf(paste(save_path_root, "_neighbourhood_dot_plot_cluster_split_Custom.pdf", sep = ""), height = 14, width = 6)
    DotPlot(combined_layers_subset, 
            features = cols_to_plot,
            group.by = "clust",
              scale = F) + coord_flip() + RotatedAxis()
    for (radii in unique(gsub(".*_", "", cols_to_plot))){
      cols_to_plot_subset <- cols_to_plot[grepl(paste("_", radii, sep = ""), cols_to_plot)]
      
      print(DotPlot(combined_layers_subset, 
              features = cols_to_plot_subset,
              group.by = "clust",
                scale = F) + coord_flip() + RotatedAxis())
    }
    dev.off()
    
    # 2.5 Save the spot id's per niche ----
    niche_spot_ids <- list()
    for (niche in names(heatmap_col_order)){
     spot_names <-  colnames(spot_neighbour_programs_df_flat_plot)
     spot_names_indicies <- heatmap_col_order[[niche]]
     spot_names <- spot_names[spot_names_indicies]
     
     niche_spot_ids[[niche]] <- spot_names
    }
    
    saveRDS(niche_spot_ids, file = paste(save_path_root, "_niche_spot_ids_split_custom_V2.rds", sep = ""))
    saveRDS(spot_neighbour_programs_df_flat_plot, file = paste(save_path_root, "_spot_neighbour_programs_df_flat_plot.rds", sep = ""))
  }
}

# 4. Gene Expression Analysis ----
# 4.1 Inside-Outside T-Cell Gene Expression Analysis ----
if (INSIDE_OUTSIDE_T_CELL){
  # 1.0 Clear existing non-function variables, set the save folders (create them if they doesn't exist) ----
  rm(list = setdiff(setdiff(ls(), lsf.str()), program_var_names))
  save_path <- paste(SAVE_FOLDER, "K", RANK, "/inside_outside/T_Cell/", sep = "")
  create_folder(save_path)
  
  # 1.1 Get the T-cells, the spot coordinates from the 3 samples, and their intersection ----
  T_cell_so <- subset(x = combined_layers, subset = CD3D | CD3E | CD3G | CD4 | CD8A | GCAR, slot = "counts")
  T_cell_spots <- rownames(T_cell_so@meta.data)
  
  # Get the spot names and merge them with the union_cell_coordinate dataframe
  T_cell_spots <- intersect(T_cell_spots, rownames(union_cell_coordinate))
  
  # 1.2 Subset the spot neighbours ----
  spot_neighbours_df_subset <- spot_neighbours_df[spot_neighbours_df$spot_of_interest %in% T_cell_spots,]
  
  # 1.3 Calculate the within and outside cells across radii ----
  combined_metadata_table_duplicate <- combined_metadata_table
  new_columns <- as.vector(unique(as.character(spot_neighbours_df_subset$r)))
  new_columns <- paste("dist", new_columns, sep = "_")
  combined_metadata_table_duplicate[, new_columns] <- "outside"
  
  # Set the neighbours of a spot to be "within", else make them "outside"
  for (radii in unique(spot_neighbours_df_subset$r)){
    # Get all the neighbours at a given radius
    radii_neighs <- unlist(spot_neighbours_df_subset[spot_neighbours_df_subset$r == radii, "neighbours"])

    # Set the column name to set cells to "within"
    radii_name <- paste("dist", radii, sep = "_")
    
    # Set the cells to "within"
    combined_metadata_table_duplicate[rownames(combined_metadata_table_duplicate) %in% radii_neighs, radii_name] <- "inside"
  }
  
  # 1.4 Merge the spots with specific gene counts ----
  count_table_subset <- as.data.frame(t(as.matrix(combined_layers@assays$Spatial.024um$counts[INSIDE_OUTSIDE_GENES$gene ,])))
  
  within_outside_df <- merge(combined_metadata_table_duplicate, count_table_subset, by = "row.names")
  colnames(within_outside_df)[1] <- "spot"
  
  within_outside_df$tissue <- paste(sapply(strsplit(within_outside_df$spot, "_"), "[", 1))
  within_outside_df$tissue <- sub("bS1", "bS", within_outside_df$tissue)
  within_outside_df$tissue <- sub("bS2", "bS", within_outside_df$tissue)
  
  within_outside_df_spread <- within_outside_df %>%
    group_by(spot) %>%
    pivot_longer(c(dist_0, dist_75, dist_150, dist_225), names_to = "radius", values_to = "location") %>%
    pivot_longer(INSIDE_OUTSIDE_GENES$gene, names_to = "gene", values_to = "expression")
  
  # 1.5 Histogram of the fraction of PDL1+ bins in the inside vs outside bins ----
  histo_plot_data <- within_outside_df_spread %>%
    group_by(tissue, radius, location, gene) %>%
    summarise(total_spots = n(), nr_gene_pos_bins = sum(expression > 0), gene_pos_bin_mean_exp = mean(expression[expression > 0])) %>%
    mutate(proportion_gene_pos_bins = nr_gene_pos_bins / total_spots, 
           proportion_gene_pos_bin_mean_exp = gene_pos_bin_mean_exp / nrow(combined_metadata_table_duplicate)) %>% 
    mutate(percent_gene_pos_bins = proportion_gene_pos_bins * 100, 
           percent_gene_pos_bin_mean_exp = proportion_gene_pos_bin_mean_exp * 100) %>%
    ungroup()
  histo_plot_data[is.na(histo_plot_data)] <- 0
  histo_plot_data$split <- paste(histo_plot_data$tissue, histo_plot_data$location)
  histo_plot_data$radius <- factor(histo_plot_data$radius, levels = c("dist_0", "dist_75", "dist_150", "dist_225"))
  
  # Modifications for plotting ----
  histo_plot_data$split <- sub("bS", "Biopsy", histo_plot_data$split)
  histo_plot_data$split <- sub("pNA", "Primary", histo_plot_data$split)
  histo_plot_data$split <- sub("xeno", "Xenograft", histo_plot_data$split)
  histo_plot_data$split <- sub("inside", "Inside", histo_plot_data$split)
  histo_plot_data$split <- sub("outside", "Outside", histo_plot_data$split)
  colnames(histo_plot_data)[12] <- "Split"
  
  histo_plot_data$radius <- sub("dist_0", "0", histo_plot_data$radius)
  histo_plot_data$radius <- sub("dist_75", "1", histo_plot_data$radius)
  histo_plot_data$radius <- sub("dist_150", "2", histo_plot_data$radius)
  histo_plot_data$radius <- sub("dist_225", "3", histo_plot_data$radius)
  histo_plot_data$radius <- factor(histo_plot_data$radius, levels = c("0", "1", "2", "3"))
  
  histo_plot_data$tissue <- sub("bS", "Biopsy", histo_plot_data$tissue)
  histo_plot_data$tissue <- sub("pNA", "Primary", histo_plot_data$tissue)
  histo_plot_data$tissue <- sub("xeno", "Xenograft", histo_plot_data$tissue)
  histo_plot_data$tissue <- factor(histo_plot_data$tissue, levels = c("Primary", "Biopsy", "Xenograft"))
  
  histo_plot_data$Split_V2 <- paste(histo_plot_data$location, histo_plot_data$radius, sep = "_") 
  
  # 1.6 Run the stats ----
  stats <- compare_means(proportion_gene_pos_bins ~ Split, histo_plot_data, 
                         group.by = "gene", method = "t.test")
  
  # 1.7 Plotting ----
  # 1.7.1 Plot the histogram ----
  # Get the range of the heatmap w/ p-values
  ht_col_range <- formatC(stats$p, format = "e", digits = 2)
  ht_col_range <- as.numeric(gsub(".*e[+|-]", "", ht_col_range))
  if (min(ht_col_range) != 0){
    if (min(ht_col_range) == max(ht_col_range)){
      col_fun <- colorRamp2(c(0, max(ht_col_range, na.rm = TRUE)+1), hcl_palette = "Blues", reverse = TRUE) 
    }else{
      col_fun <- colorRamp2(c(0, max(ht_col_range, na.rm = TRUE)), hcl_palette = "Blues", reverse = TRUE) 
    }
  }else{
    col_fun <- colorRamp2(range(ht_col_range, na.rm = TRUE), hcl_palette = "Blues", reverse = TRUE)
  }
  
  # Make the plots
  pdf(paste(save_path, "Inside_Outside_Gene_Expression_Histogram_T_Cell_Nominal_P.pdf", sep = ""), width = 30, height = 10)
  for (i in seq_len(nrow(INSIDE_OUTSIDE_GENES))){
    pair <- INSIDE_OUTSIDE_GENES$category[i]
    gene <- INSIDE_OUTSIDE_GENES$gene[i]
    
    # Subset the histogram data to the gene of interest
    temp_plot_data <- histo_plot_data[histo_plot_data$gene == gene,]
    
    # Create the histogram
    gp <- ggplot(temp_plot_data, aes(x = radius, y = proportion_gene_pos_bins, fill = Split_V2)) +
      geom_bar(stat="identity", color = "black", position = "dodge", width = .8) + theme_classic2() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, colour="black"), 
            axis.text.y=element_text(colour="black"), text = element_text(size = 20), 
            legend.position = "bottom") +
      labs(x = "Distance From CD3D/E/G, CD4, CD8A, and GCAR Bin", 
           y = "Proportion of Gene+ Spots\nof Total Spots Per Condition",
           title = paste(gene, pair, sep = ": "),
           subtitle = "Fraction of Adjacent and Distant Bins From CD3D/E/G, CD4, CD8A, and GCAR+ Cells") + 
      scale_y_continuous(expand = c(0,0)) +  facet_wrap( ~ tissue + location, ncol = 6) +
      scale_fill_manual(values = c("#2b465e", "#607486", "#95a2af", "#cad1d7", 
                                   "#d2c691", "#ddd4ad", "#e9e2c8", "#f4f1e3")) +
      guides(fill = guide_legend(ncol = 4, byrow = TRUE))
    
    # Create the P-Value Heatmap in the same figure
    temp_stats_data <- as.data.frame(stats[stats$gene == gene, c(3, 4, 5)])
    temp_stats_data_rev <- temp_stats_data
    colnames(temp_stats_data_rev)[1:2] <- c("group2", "group1")
    temp_stats_data <- rbind(temp_stats_data, temp_stats_data_rev)
    
    temp_stats_data_spread <- dcast(temp_stats_data, group1 ~ group2, value.var = "p")
    rownames(temp_stats_data_spread) <- temp_stats_data_spread$group1
    temp_stats_data_spread <- temp_stats_data_spread[-1]
    temp_stats_data_spread[is.na(temp_stats_data_spread)] <- 0
    temp_stats_data_spread <- as.matrix(temp_stats_data_spread)
    temp_stats_data_spread <- formatC(temp_stats_data_spread, format = "e", digits = 2)
    
    temp_stats_data_spread <- gsub(".*e[+|-]", "", temp_stats_data_spread)
    
    temp_stats_data_spread <- as.data.frame(temp_stats_data_spread)
    stats_row_name_store <- rownames(temp_stats_data_spread)
    
    temp_stats_data_spread <- sapply(temp_stats_data_spread, as.numeric)
    rownames(temp_stats_data_spread) <- stats_row_name_store
    
    temp_stats_data_spread <- as.matrix(temp_stats_data_spread)
    diag(temp_stats_data_spread) <- NA
    
    data_cols <- colnames(temp_stats_data_spread)
    data_cols_order <- c(data_cols[3], data_cols[1], data_cols[5], 
                         data_cols[4], data_cols[2], data_cols[6])
    data_cols_order <- data_cols_order[!is.na(data_cols_order)]
    
    temp_stats_data_spread <- temp_stats_data_spread[data_cols_order, data_cols_order]
    
    # Make the heatmap
    htp_name <- paste(gene, ": ", pair, ":\nAbsolute Value of The P-Value Scientific Notation Exponent\n[T.test]", sep = "")
    htp <- Heatmap(temp_stats_data_spread,
                   cluster_columns = FALSE, cluster_rows = FALSE,
                   show_row_names = TRUE, show_column_names = TRUE,
                   col = col_fun, na_col = "black")
    htp_grb <- grid.grabExpr(draw(htp, column_title=htp_name, padding = unit(c(3, 3, 3, 3), "cm")))

    grid.arrange(gp, htp_grb, ncol = 2)
  }
  dev.off()

}

# 4.2 Inside-Outside Niche Gene Expression Analysis ----
if (INSIDE_OUTSIDE_NICHE){
  # 2.0 Clear existing non-function variables, set the save folders (create them if they doesn't exist) ----
  rm(list = setdiff(setdiff(ls(), lsf.str()), program_var_names))
  save_path <- paste(SAVE_FOLDER, "K", RANK, "/inside_outside/niche/custom_split/", sep = "")
  create_folder(save_path)
  
  # 2.1 Load the data ----
  niche_spot_ids <- readRDS(file = "02_Visium_HD/02_analysis/tmp/K15/neighbourhood_analysis/Programs/GEP_9_10/GEP_9_10_T0.05_Binarized_T0.05_niche_spot_ids_split_custom_V2.rds")
  spot_neighbour_programs_df_flat_plot <- readRDS(file = "02_Visium_HD/02_analysis/tmp/K15/neighbourhood_analysis/Programs/GEP_9_10/GEP_9_10_T0.05_Binarized_T0.05_spot_neighbour_programs_df_flat_plot.rds")
  
  # Subset to the bS1 and bS2 spots and remove any spots that belong to niche 1
  niche_spot_ids_custom <- niche_spot_ids[c(1,4,5)]
  
  # 2.2 Contingency table ----
  # Retrieve the R0 data
  temporary_df <- data.frame()
  for (niche in names(niche_spot_ids_custom)){
    # 4.1 Get the T-cells, the spot coordinates from the 3 samples, and their intersection ----
    niche_spots <- niche_spot_ids_custom[[niche]]
    dp_spots <- niche_spots
    
    # Get the spot names and merge them with the union_cell_coordinate dataframe
    dp_spots <- intersect(dp_spots, rownames(union_cell_coordinate))
    
    # 4.2 Subset the spot neighbours ----
    spot_neighbours_df_subset <- spot_neighbours_df[spot_neighbours_df$spot_of_interest %in% dp_spots,]
    
    # 4.3 Calculate the inside and outside cells across radii ----
    combined_metadata_table_duplicate <- combined_metadata_table
    new_columns <- as.vector(unique(as.character(spot_neighbours_df_subset$r)))
    new_columns <- paste("dist", new_columns, sep = "_")
    combined_metadata_table_duplicate[, new_columns] <- "outside"
    
    # Set the neighbours of a spot to be "inside", else make them "outside"
    for (radii in unique(spot_neighbours_df_subset$r)){
      # Get all the neighbours at a given radius
      radii_neighs <- unlist(spot_neighbours_df_subset[spot_neighbours_df_subset$r == radii, "neighbours"])
      
      # Set the column name to set cells to "inside"
      radii_name <- paste("dist", radii, sep = "_")
      
      # Set the cells to "inside"
      combined_metadata_table_duplicate[rownames(combined_metadata_table_duplicate) %in% radii_neighs, radii_name] <- "inside"
    }
    
    # Remove the other niche spots from this analysis
    niche_spots_to_remove <- unlist(niche_spot_ids_custom[names(niche_spot_ids_custom) != niche])
    combined_metadata_table_duplicate <- combined_metadata_table_duplicate[!rownames(combined_metadata_table_duplicate) 
                                                                           %in% niche_spots_to_remove,]
    
    # 4.4 Merge the spots with specific gene counts ----
    count_table_subset <- as.data.frame(t(as.matrix(combined_layers@assays$Spatial.024um$counts[INSIDE_OUTSIDE_GENES$gene ,])))
    
    within_outside_df <- merge(combined_metadata_table_duplicate, count_table_subset, by = "row.names")
    colnames(within_outside_df)[1] <- "spot"
    
    within_outside_df$tissue <- paste(sapply(strsplit(within_outside_df$spot, "_"), "[", 1))
    within_outside_df$tissue <- sub("bS1", "bS", within_outside_df$tissue)
    within_outside_df$tissue <- sub("bS2", "bS", within_outside_df$tissue)
    
    # Only keep the tissue of interest
    within_outside_df <- within_outside_df[within_outside_df$tissue == strsplit(sub(",", "_", niche), "_")[[1]][2],]
    
    # 4.4.1 Spread the data ----
    within_outside_df_spread <- within_outside_df %>%
      group_by(spot) %>%
      pivot_longer(c(dist_0, dist_75, dist_150, dist_225), names_to = "radius", values_to = "location") %>%
      pivot_longer(INSIDE_OUTSIDE_GENES$gene, names_to = "gene", values_to = "expression")
    
    # 4.5 Histogram of the fraction of PDL1+ bins in the inside vs outside bins ----
    histo_plot_data <- within_outside_df_spread %>%
      group_by(tissue, radius, location, gene) %>%
      summarise(total_spots = n(), nr_gene_pos_bins = sum(expression > 0), gene_pos_bin_mean_exp = mean(expression[expression > 0])) %>%
      mutate(proportion_gene_pos_bins = nr_gene_pos_bins / total_spots, 
             proportion_gene_pos_bin_mean_exp = gene_pos_bin_mean_exp / nrow(combined_metadata_table_duplicate)) %>% 
      mutate(percent_gene_pos_bins = proportion_gene_pos_bins * 100, 
             percent_gene_pos_bin_mean_exp = proportion_gene_pos_bin_mean_exp * 100) %>%
      ungroup()
    histo_plot_data[is.na(histo_plot_data)] <- 0
    colnames(histo_plot_data)[1] <- "tissue"
    histo_plot_data$split <- paste(histo_plot_data$tissue, histo_plot_data$location)
    histo_plot_data$radius <- factor(histo_plot_data$radius, levels = c("dist_0", "dist_75", "dist_150", "dist_225"))
    
    # Modifications for plotting ----
    histo_plot_data$split <- sub("bS", "Biopsy", histo_plot_data$split)
    histo_plot_data$split <- sub("pNA", "Primary", histo_plot_data$split)
    histo_plot_data$split <- sub("xeno", "Xenograft", histo_plot_data$split)
    histo_plot_data$split <- sub("inside", "Inside", histo_plot_data$split)
    histo_plot_data$split <- sub("outside", "Outside", histo_plot_data$split)
    colnames(histo_plot_data)[12] <- "Split"
    
    histo_plot_data$radius <- sub("dist_0", "0", histo_plot_data$radius)
    histo_plot_data$radius <- sub("dist_75", "1", histo_plot_data$radius)
    histo_plot_data$radius <- sub("dist_150", "2", histo_plot_data$radius)
    histo_plot_data$radius <- sub("dist_225", "3", histo_plot_data$radius)
    histo_plot_data$radius <- factor(histo_plot_data$radius, levels = c("0", "1", "2", "3"))
    
    histo_plot_data$tissue <- sub("bS", "Biopsy", histo_plot_data$tissue)
    histo_plot_data$tissue <- sub("pNA", "Primary", histo_plot_data$tissue)
    histo_plot_data$tissue <- sub("_", " Niche ", histo_plot_data$tissue)
    histo_plot_data$tissue <- sub("xeno", "Xenograft", histo_plot_data$tissue)
    
    histo_plot_data$Split_V2 <- paste(histo_plot_data$location, histo_plot_data$radius, sep = "_")
    
    histo_plot_data <- histo_plot_data[histo_plot_data$radius == 0,]
    histo_plot_data$niche <- niche
    
    temporary_df <- rbind(temporary_df, histo_plot_data)
  }
  
  # 2.3 Data Manipulation ----
  temporary_df_final <- temporary_df[,c(14,3,4,6,5)]
  colnames(temporary_df_final)[4] <- "positive"
  temporary_df_final$negative <- temporary_df_final$total_spots - temporary_df_final$positive
  temporary_df_final <- temporary_df_final[,c(1:4,6,5)]
  temporary_df_final <- temporary_df_final[temporary_df_final$location == "inside",]
  temporary_df_final$niche <- sub(",", "_", temporary_df_final$niche)
  temporary_df_final$niche <- paste(temporary_df_final$niche, temporary_df_final$location, sep = "_")
  write.table(temporary_df_final, paste(save_path, "R0_Final_Table_Custom_Split.tsv", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
  
  # 2.4 Plotting ----
  f_pw_df <- data.frame()
  pdf(paste(save_path, "R0_Gene_Expression_Fisher_Exact_Nominal_P_Custom_Split.pdf", sep = ""), width = 15, height = 15)
  for (i in seq_len(nrow(INSIDE_OUTSIDE_GENES))){
    pair <- INSIDE_OUTSIDE_GENES$category[i]
    gene <- INSIDE_OUTSIDE_GENES$gene[i]
    
    foo <- temporary_df_final
    foo$niche <- sub(",", "_", foo$niche)
    foo$niche <- paste(foo$niche, foo$location, sep = "_")
    foo <- foo[foo$gene == gene,]
    foo <- as.data.frame(foo); rownames(foo) <- foo$niche
    foo <- foo[,-c(1,2,3,6)]
    
    temp <- fisher.test(foo, workspace = 2e9)
    temp_ph <- pairwise_fisher_test(foo, p.adjust.method = "fdr")
    
    # Create the P-Value Heatmap figure ----
    temp_stats_data <- temp_ph
    temp_stats_data_rev <- temp_stats_data
    colnames(temp_stats_data_rev)[1:2] <- c("group2", "group1")
    temp_stats_data <- rbind(temp_stats_data, temp_stats_data_rev)
    temp_temp_ph <- temp_stats_data
    temp_temp_ph$gene <- gene
    f_pw_df <- rbind(f_pw_df, temp_temp_ph)
    
    temp_stats_data_spread <- dcast(temp_stats_data, group1 ~ group2, value.var = "p")
    rownames(temp_stats_data_spread) <- temp_stats_data_spread$group1
    temp_stats_data_spread <- temp_stats_data_spread[-1]
    temp_stats_data_spread[is.na(temp_stats_data_spread)] <- 0
    temp_stats_data_spread <- as.matrix(temp_stats_data_spread)
    temp_stats_data_spread <- formatC(temp_stats_data_spread, format = "e", digits = 2)
    
    temp_stats_data_spread <- gsub(".*e[+|-]", "", temp_stats_data_spread)

    temp_stats_data_spread <- as.data.frame(temp_stats_data_spread)
    stats_row_name_store <- rownames(temp_stats_data_spread)
    
    temp_stats_data_spread <- sapply(temp_stats_data_spread, as.numeric)
    rownames(temp_stats_data_spread) <- stats_row_name_store

    temp_stats_data_spread <- as.matrix(temp_stats_data_spread)
    diag(temp_stats_data_spread) <- NA
    
    # Make the heatmap ----
    if (min(range(temp_stats_data_spread, na.rm = TRUE)) == max(range(temp_stats_data_spread, na.rm = TRUE))){
      col_fun <- colorRamp2(c(min(temp_stats_data_spread, na.rm = TRUE), min(temp_stats_data_spread, na.rm = TRUE) + 1), hcl_palette = "Blues", reverse = TRUE)
    } else{col_fun <- colorRamp2(range(temp_stats_data_spread, na.rm = TRUE), hcl_palette = "Blues", reverse = TRUE)}
    
    htp_name <- paste(gene, "- ", pair, ":",
                      "\nAbsolute Value of The P-Value Scientific Notation Exponent [Capped at 10] \n[Post-Hoc Pairwise Fisher Exact P-Value]", 
                      "\nFisher Exact P-Value: ", temp$p.value, sep = "")
    htp <- Heatmap(temp_stats_data_spread,
                   cluster_columns = FALSE, cluster_rows = FALSE,
                   show_row_names = TRUE, show_column_names = TRUE,
                   col = col_fun, na_col = "black")

    draw(htp, column_title = htp_name, padding = unit(c(3, 3, 3, 3), "cm"))
  }
  dev.off()
  write.table(f_pw_df, paste(save_path, "R0_Fischer_Pairwise_Stat_Table_Custom_Split.tsv", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
  
  temporary_df_final$prop <- (temporary_df_final$positive / temporary_df_final$total_spots) * 100
  
  temporary_df_final$niche[temporary_df_final$niche == "1_bS_inside"] <- "1 (fibroblast-rich)"
  temporary_df_final$niche[temporary_df_final$niche == "2_bS_6_7_positive_inside"] <- "2 (endothelial / pericytes)"
  temporary_df_final$niche[temporary_df_final$niche == "3_bS_remainder_inside"] <- "3 (tumor)"
  
  temporary_df_final <- merge(temporary_df_final, INSIDE_OUTSIDE_GENES, by = "gene")
  
  cols <- c("#43914e", "#3f509e", "#cf853a")
  names(cols) <- unique(temporary_df_final$niche)
  pdf(paste(save_path, "Proportion_Gene_Positive_Bins_Custom_Split.pdf", sep = ""), width = 30, height = 8)
  gp1 <- ggplot(temporary_df_final, aes(x = niche, y = prop, color = niche, fill = niche)) + 
    geom_bar(stat="identity") + 
    theme_classic2() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    labs(title = "Percent of Bins Positive For Gene Expression Across Neighborhoods.") + 
    scale_color_manual(values = cols) + scale_fill_manual(values = cols) +
    scale_y_continuous(expand = c(0, 0)) + 
    facet_wrap(~ category + gene, scales = "free_y", ncol = length(unique(temporary_df_final$gene)))
  print(gp1)
  dev.off()
}
