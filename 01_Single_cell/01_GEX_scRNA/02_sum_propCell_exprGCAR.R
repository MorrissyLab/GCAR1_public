################################################
#### [GCAR1] scRNA-seq data analysis:       ####
#### summarize propCell_exprGCAR            ####
#### - figure: [ supplementary figure #10 ] ####
################################################
# version: 0.1.1 (2024-12-13)
# by hyojinsong
today <- Sys.Date()
# ---

#### packages ####
mypackages <- c(
  "BiocManager", "ggplot2", "plyr", "dplyr", "reshape2", "RColorBrewer", "pals", # "readr", "parallel",
  "Seurat", "hdf5r", "sctransform", "glmGamPoi", #for Seurat
  "ggalluvial"
)

lapply(mypackages, library, character.only = TRUE)
# ---


#### running environment ####
# mode <- "ARC"
mode <- "local"

#### running switch ####
tcell_only <- T #either all(GEX) or subset(VDJ)
if (!isTRUE(tcell_only)) {
  dtype <- "GEX-all"
} else if (isTRUE(tcell_only)) {
  dtype <- "GEX-subsetVDJ"
}

mk_plt1 <- T
mk_plt1a <- T
mk_plt1b <- T
# ---

#### project ####
prj <- "mrg7"


#### output directory ####
mydir_res <- paste0(dir_res_scm,today,"_SCim_",prj,"/")
if (!dir.exists(mydir_res)) {
  dir.create(mydir_res, recursive = T)
  print(paste0("[DONE] create a directory: ",mydir_res))
} else {
  print(paste0("[SKIP] directory exists: ",mydir_res))
}


#### input data ####
## input 1 | <obs.csv> file - from SCimilarity "constrained" approach
# "/Users/hyojinsong/work/uofc_iSARP/result/SCimilarity/SCim_GCAR1_PT01/mrg7/adata-to-csv_constrain_res2.4_VDJonly/obs.csv"
fn_obs_c <- paste0(dir_res_scm,"SCim_GCAR1_PT01/mrg7/adata-to-csv_constrain_res2.4_VDJonly/","obs.csv")
obs_c <- read.csv(fn_obs_c, header = T, row.names = 1)
head(obs_c,3)
dim(obs_c) #77651 20(GEX-subsetVDJ;nBatch3)

## from the [VDJ] data ##
## contig group: [all cells] ##
## load the [scRepertoire] "clone" data
fn_cln_all <- paste0("/Users/hyojinsong/work/uofc_iSARP/result/scRepertoire/clones_GCAR1_PT01_VDJ_all__nSmp7.csv")
cln_all <- read.csv(file = fn_cln_all, header = T, row.names = 1)
head(cln_all,3)
dim(cln_all) #81645 7(VDJ;nBatch3)
# as.data.frame(table(cln_all$group)) #checkpoint

## add a column to the "GEX" data frame
obs_c$seurat_clusters_VDJsubset <- cln_all$seurat_clusters[match(rownames(obs_c),rownames(cln_all))]
head(obs_c,3) #checkpoint



#### data summary ####
# ---
## step 1 | quantify the number of cells w/ GCAR expression ####
## make a data frame for GCAR-expressing cells
df_expr_gcar <- as.data.frame(table(obs$orig.ident, obs$expr_GCAR))
dim(df_expr_gcar) #14 3
head(df_expr_gcar)
colnames(df_expr_gcar) <- c("sample","expr_GCAR","nCell")
tail(df_expr_gcar)

## add a column: [nCell_eachSmp]
df_expr_gcar$nCell_eachSmp <- df_cells$nCell_eachSmp[match(df_expr_gcar$sample, df_cells$sample)]
head(df_expr_gcar)
subset(df_expr_gcar, sample=="Product") #11520
subset(df_expr_gcar, sample=="D15") #23420

## add a column: [propCell]
df_expr_gcar <- df_expr_gcar %>%
  group_by(sample) %>%
  mutate(
    propCell=(nCell/nCell_eachSmp)*100
  ) %>%
  as.data.frame()
head(df_expr_gcar)
subset(df_expr_gcar, sample=="Product") #11520
subset(df_expr_gcar, sample=="D15") #23420

## set the order of [sample] variables in chronological order
if (nSmp==7) {
  df_expr_gcar <- df_expr_gcar %>%
    mutate(
      order=case_when(
        sample=="Product" ~ 1,
        sample=="Tcells" ~ 2,
        sample=="D08" ~ 3,
        sample=="D15" ~ 4,
        sample=="D22" ~ 5,
        sample=="D28" ~ 6,
        sample=="D42" ~ 7
      )) %>%
    arrange(order,expr_GCAR) %>%
    as.data.frame()
} else if (nSmp==6) {
  df_expr_gcar <- df_expr_gcar %>%
    mutate(
      order=case_when(
        sample=="Product" ~ 1,
        sample=="Tcells" ~ 2,
        sample=="D08" ~ 3,
        sample=="D15" ~ 4,
        sample=="D22" ~ 5,
        sample=="D28" ~ 6
      )) %>%
    arrange(order,expr_GCAR) %>%
    as.data.frame()
} #closer: if (nSmp==7) {
head(df_expr_gcar,3)
tail(df_expr_gcar,3)

## convert the [expr_GCAR] variables as "factors"
# df_expr_gcar$expr_GCAR <- factor(df_expr_gcar$expr_GCAR, levels=c("pos","neg"))
df_expr_gcar$expr_GCAR <- factor(df_expr_gcar$expr_GCAR, levels=c("neg","pos"))
# df_expr_gcar$expr_GCAR <- as.factor(df_expr_gcar$expr_GCAR)
table(df_expr_gcar$expr_GCAR)

## modify the [sample] variable names
df_expr_gcar <- df_expr_gcar %>%
  mutate(
    sample_new=case_when(
      sample=="Product" ~ "Enriched", #"Apheresis",
      sample=="Tcells" ~ "Harvest", #"GCAR1",
      sample!="Product" & sample!="Tcells" ~ sample
    )
  ) %>%
  as.data.frame()
table(df_expr_gcar$sample_new)
## modify the order of the "sample_new" variables
if (nSmp==7) {
  df_expr_gcar$sample_new <- factor(df_expr_gcar$sample_new, 
                                    levels = c(
                                      "Enriched","Harvest",
                                      "D08","D15","D22","D28","D42"
                                    ))
} else if (nSmp==6) {
  df_expr_gcar$sample_new <- factor(df_expr_gcar$sample_new, 
                                    levels = c(
                                      "Enriched","Harvest",
                                      "D08","D15","D22","D28"
                                    ))
}
table(df_expr_gcar$sample_new)

## count the number of samples
# nSmp <- length(unique(df_expr_gcar$sample))
# ---


#### visualization ####
# ---
## plot 1: summarized proportions of the cluster-based predicted cell types in each sample (stacked bar plot) ####
if (mk_plt1) {
  
  mytitle <- paste0("[",prjid,"/",smp,"] Total number of cells")
  mysubt <- paste0("mean: ",mean_nCell_eachSmp," | median: ",median_nCell_eachSmp)
  # print(mytitle)
  # print(mysubt)
  mylab_x <- "Time point"
  mylab_y <- "Total number of cells"
  
  ggsbar1 <-
    ggplot(data = df_cells, aes(x = sample_new, y=nCell_eachSmp)) +
    theme_bw() +
    geom_bar(
      stat = "identity",
      width = .38,
      fill = blueColor_cool
    ) +
    # geom_line(
    #   aes(x = sample_new, y=nCell_eachSmp),
    #   stat="identity",
    #   group=1,
    #   col=greyColor_light
    #   ) +
    theme(
      legend.position = "none", #"right" #"bottom"
      # axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    ) +
    scale_y_continuous(expand = c(0.025,0.025)) + #to start from zero (w/o additional space on the beginning of the y-axis)
    ylab(label = mylab_y) +
    xlab(label = mylab_x) +
    labs(
      title = mytitle,
      subtitle = mysubt
    )
  ## save as PDF: [ggsbar1]
  mywidth <- 4.3
  myheight <- 3
  fn_pdf_ggsbar1 <- paste0(mydir_res,"1_sBar_totalNumCell__",prj,"_nSmp",nSmp, ".pdf")
  ggsave(filename = fn_pdf_ggsbar1, plot = ggsbar1, device = "pdf", width = mywidth, height = myheight, units = "in", dpi = 300)
  print(paste0("[DONE] save as PDF: ",fn_pdf_ggsbar1))
  
} #closer: if (mk_plt1) {
# ---

# ---
## plot 1a: proportion of GCAR-pos/neg cells in each sample - in the "merged" data (stacked bar plot) ####
if (isTRUE(mk_plt1a)) {
  
  mytitle <- paste0("[",prjid,"/",smp,"] Proportion of GCAR-pos/neg cells")
  mysubt <- paste0("nSmp: ",nSmp," | {GCAR-pos = if nCell >=1}")
  # print(mytitle)
  # print(mysubt)
  mylab_x <- "Time point" #Sample"
  mylab_y <- "Proportion of cells (%)"
  lgd_fill <- "GCAR1 expression"
  # mylab_EleA <- "GCAR1-positive"
  # mylab_EleB <- "GCAR1-negative"
  
  ggsbar1a <-
    ggplot(data = df_expr_gcar, aes(x = reorder(sample_new,order), y=propCell, fill=expr_GCAR)) +
    theme_bw() +
    geom_bar(
      stat = "identity",
      width = .55
    ) +
    # geom_flow(aes(alluvium=expr_GCAR), curve_type = "linear", alpha=0.25) +
    ## white border around the bar plots
    geom_flow(aes(alluvium=expr_GCAR), curve_type = "linear", alpha=0.25, color="white") +
    geom_col(width = .55, color = "white") +
    theme(
      # axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "bottom" #"right" #"bottom"
    ) +
    scale_fill_manual(values = exprColor,
                      ## change the order of variables in the legend
                      breaks = c("pos","neg")
                      # labels = c(mylab_EleA,mylab_EleB)
    ) +
    scale_y_continuous(expand = c(0.01,0.01)) + #to start from zero (w/o additional space on the beginning of the y-axis)
    geom_text(
      # aes(label=ifelse(expr_GCAR=="pos", round(propCell,2), "")),
      aes(label=ifelse(expr_GCAR=="pos", paste0(round(propCell,2),"\n(",nCell,"\n/",nCell_eachSmp,")"), "")),
      vjust=-.5, #neg(upward) #pos(downward)
      position = "stack",
      size = 2.5,
      color=greyColor_dark
    ) +
    ylab(label = mylab_y) +
    xlab(label = mylab_x) +
    labs(
      title = mytitle,
      subtitle = mysubt,
      fill = lgd_fill
    )
  ## save as PDF: [ggsbar1a]
  mywidth <- 5
  myheight <- 4.5
  fn_pdf_ggsbar1a <- paste0(mydir_res,"1a_sBar_propCell_exprGCAR__",prj,"_nSmp",nSmp, ".pdf")
  ggsave(filename = fn_pdf_ggsbar1a, plot = ggsbar1a, device = "pdf", width = mywidth, height = myheight, units = "in", dpi = 300)
  print(paste0("[DONE] save as PDF: ",fn_pdf_ggsbar1a))
  
} #closer: if (isTRUE(mk_plt1a)) {
# ---



# ---
## plot 1b: proportion of GCAR-pos/neg cells in each sample - in the "merged" GEX+VDJ data (grouped stacked bar plot) ####
if (isTRUE(mk_plt1b)) {
  
  ## ongoing (need-to-adjust):
  mytitle <- paste0("[",prjid,"/",smp,"] Proportion of GCAR-pos/neg cells in GEX+VDJ")
  mysubt <- paste0("nSmp: ",nSmp," | {GCAR-pos = if nCell >=1}")
  # print(mytitle)
  # print(mysubt)
  mylab_x <- "Time point" #Sample"
  mylab_y <- "Proportion of cells (%)"
  lgd_fill <- "GCAR1 expression"
  # mylab_EleA <- "GCAR1-positive"
  # mylab_EleB <- "GCAR1-negative"

  ggsbar1b <-
    ggplot(data = df_expr_gcar, aes(x = reorder(sample_new,order), y=propCell, fill=expr_GCAR)) +
    theme_bw() +
    geom_bar(
      stat = "identity",
      width = .55,
      position = position_dodge()
    ) +
    geom_col(width = .55, color = "white") +
    theme(
      strip.background = element_rect(fill = "transparent"),
      # axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "bottom" #"right" #"bottom"
    ) +
    scale_fill_manual(values = exprColor,
                      ## change the order of variables in the legend
                      breaks = c("pos","neg")
                      # labels = c(mylab_EleA,mylab_EleB)
    ) +
    scale_y_continuous(expand = c(0.01,0.01)) + #to start from zero (w/o additional space on the beginning of the y-axis)
    geom_text(
      # aes(label=ifelse(expr_GCAR=="pos", round(propCell,2), "")),
      aes(label=ifelse(expr_GCAR=="pos", paste0(round(propCell,2),"\n(",nCell,"\n/",nCell_eachSmp,")"), "")),
      vjust=-.5, #neg(upward) #pos(downward)
      position = "stack",
      size = 2.5,
      color=greyColor_dark
    ) +
    facet_wrap(~sample_new, scales = "free_x", nrow = 1) +
    ylab(label = mylab_y) +
    xlab(label = mylab_x) +
    labs(
      title = mytitle,
      subtitle = mysubt,
      fill = lgd_fill
    )
  ## save as PDF: [ggsbar1b]
  mywidth <- 5
  myheight <- 4.5
  fn_pdf_ggsbar1b <- paste0(mydir_res,"1b_gsBar_propCell_exprGCAR__",prj,"_nSmp",nSmp, ".pdf")
  ggsave(filename = fn_pdf_ggsbar1b, plot = ggsbar1b, device = "pdf", width = mywidth, height = myheight, units = "in", dpi = 300)
  print(paste0("[DONE] save as PDF: ",fn_pdf_ggsbar1b))
  
} #closer: if (isTRUE(mk_plt1b)) {
# ---

