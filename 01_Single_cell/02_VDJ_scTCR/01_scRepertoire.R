######################################
#### single cell VDJ-seq analysis ####
#### - project: [GCAR1]           ####
######################################
# version: 2.27.4 (2024-12-13)
# by hyojinsong
today <- Sys.Date()
# ---

#### PACKAGES ####
mypackages <- c(
  "BiocManager", "ggplot2", "plyr", "dplyr", "reshape2", "RColorBrewer", "pals",
  "scRepertoire",
  "ggalluvial",
  "stringr",
  "ComplexHeatmap"
  # "ggrepel"
)

# ## [scRepertoire] ##
# devtools::install_github("ncborcherding/scRepertoire", dependencies = T, force = T)
# detach(package:scRepertoire, unload = T)
# library(scRepertoire)

lapply(mypackages, library, character.only = TRUE)
# ---


#### project ####
prj <- "GCAR1_PT01_VDJ"
if (prj=="GCAR1_PT01_GEX_VDJ" | prj=="GCAR1_PT01_GEX" | prj=="GCAR1_PT01_VDJ") {
  prjid <- "GCAR1_PT01"
  sp <- "hs38_GCAR"
}
print(paste0("=== prj: [",prj,"] | prjid: [",prjid,"] | species: [",sp,"] ==="))


#### running switch ####
mk_mdata <- F

## [scRepertoire] ##
anlz_vdj <- T #[master-running-switch]
exp_cln <- T

# ---
## final ##
mk_plt_final <- T
tab_cln__sel <- T
tab_cln__sel__topClones <- F
tab_cln__sel__allClones <- T

plt_4_n6 <- T #6-sample (except for [Enriched])

mk_sum_tab <- T
# ---

## semi-final ##
## - for "cluster-of-interest"
final_plt_dns <- F
final_plt_vbox <- T

## final 
plt_vbox_ccoi <- F
plt_vbox_vdj <- T
# ---

tab_cln__sel <- T
# ===


#### directory ####
homedir <- "/Users/your_username/"
dir_data <- paste0(homedir,"data/",prjid,"/")
dir_res <- paste0(homedir,"result/")
dir_res_cr <- paste0(dir_res,"cellranger/",prjid,"/")
dir_res_scm <- paste0(dir_res,"SCimilarity/")
dir_res_scb <- paste0(dir_res,"Scrublet/")

#### [scRepertoire] ####
dir_res_scR <- paste0(dir_res,"scRepertoire/")
setwd(dir = dir_res_scR)

## output directory ##
mydir_res <- paste0(dir_res_scR,today,"_scRepertoire","_",prj,"/")
if (!dir.exists(mydir_res)) {
  dir.create(mydir_res, recursive = T)
  print(paste0("[DONE] create a directory: ",mydir_res))
} else {
  print(paste0("[SKIP] directory exists: ",mydir_res))
}


#### metadata ####
## make 
if (isTRUE(mk_mdata)) {
  ## initialize metadata
  mdata <- data.frame(
    smp_id=c("Enriched_Apheresis","GCAR1_Product","D08","D15","D22","D28","D42")
  )
  ## add columns:
  mdata <- mdata %>%
    mutate(
      smpType=case_when(
        smp_id=="Enriched_Apheresis" ~ "Enriched_Apheresis", #"untransducedCAR"
        smp_id=="GCAR1_Product" ~ "GCAR1_Product", #"CART"
        smp_id=="D08"|smp_id=="D15"|smp_id=="D22"|smp_id=="D28"|smp_id=="D42" ~ "afterInfusion"
      ),
      smp_name=case_when(
        smp_id=="Enriched_Apheresis" ~ "CD4_CD8_enriched_product_untransducedCARs",
        smp_id=="GCAR1_Product" ~ "CART_product_beforeInfusion",
        smp_id=="D08" ~ "day08_afterInfusion",
        smp_id=="D15" ~ "day15_afterInfusion",
        smp_id=="D22" ~ "day22_afterInfusion",
        smp_id=="D28" ~ "day28_afterInfusion",
        smp_id=="D42" ~ "day42_afterInfusion",
      )
    ) %>%
    as.data.frame()
  head(mdata,3)
  tail(mdata,3)
  dim(mdata) #7 3
  
  fn_mdata <- paste0(dir_data,"metadata_",prjid,"_n",nrow(mdata),".txt")
  if (!file.exists(fn_mdata)) {
    write.table(x = mdata, file = fn_mdata, sep = "\t", col.names = T, row.names = F, quote = F)
    print(paste0("[DONE] create metadata: ",fn_mdata))
  } else {
    print(paste0("[DONE] metadata TXT file exists: ",fn_mdata))
  }
  
} else if (!isTRUE(mk_mdata)) {
  ## load the metadata
  fn_mdata <- paste0(dir_data, "metadata_GCAR1_PT01_n7.txt")
  if (file.exists(fn_mdata)) {
    mdata <- read.table(file = fn_mdata, sep = "\t", header = T, row.names = NULL, quote = NULL)
    print(paste0("[DONE] load the metadata: ", fn_mdata))
  } else {
    print(paste0("[CHECK] metadata does not exist: ", fn_mdata))
  }
  
} #closer: if (isTRUE(mk_mdata)) {


## count the number of samples
# nSmp <- length(mdata$smp_id)

# ---
## load the metadata - extracted from the Seurat Object (i.e., [sobj_flt])
fn_txt_md <- paste0(homedir,"collaboration/RCCI__GCAR1/","03a_result_GEX/01_Seurat_basic/","metadata_mrg7_sobj_GCAR1_PT01_GEX-subsetVDJ_res2.4",".txt")
txt_md <- read.table(file = fn_txt_md, sep = "\t", header = T, row.names = 1)
dim(txt_md) #77651 12
print(paste0("--- nRow/nCellBarcode: [",nrow(txt_md),"] ---"))
head(txt_md,3)

## check GCAR-pos/neg cells
table(txt_md$expr_GCAR)
subset(txt_md, is.na(expr_GCAR))
# ---

# ===

#### colour palette ####
## [https://blog.r-project.org/2019/04/01/hcl-based-color-palettes-in-grdevices/]
# hcl.pals()
pal_basic <- "Temps"
pal_cln <- "Zissou 1"
pal_clnH <- "Oslo"
pal_hmp <- "Blue-Red 3" #"YlOrRd" #"Oslo" #"Blue-Red 3"
pal_aa <- "Temps" #"YlOrRd" #"Oslo" #"Blue-Red 3"
pal_pct <- "YlGnBu" #"Oslo"(previously)
pal_smp <- "Zissou 1"
## longitudinal samples ##
col_gcar_p <- "#666666" #dark_grey
col_gcar_t <- "#bf5b17" #brown
col_gcar_d08 <- "#a6d854" #light_green
col_gcar_d15 <- "#e7298a" #pink
col_gcar_d22 <- "#e6ab02" #yellow
col_gcar_d28 <- "#1b9e77" #dark_green
col_gcar_d42 <- "#386cb0" #blue
smpColor <- setNames(c(col_gcar_p,col_gcar_t,col_gcar_d08,col_gcar_d15,col_gcar_d22,col_gcar_d28,col_gcar_d42), c("Enriched_Apheresis","GCAR1_Product","D08","D15","D22","D28","D42")) #sample_new

## UMAP plot - dots ##
col_na <- "#f0f0f0" #"#d9d9d9"(Greys: 3-of-9) "#f0f0f0" (Greys: 2-of-9)


########################
#### [scRepertoire] ####
########################
## PART 1 - TCR analysis using scTCR-seq data ####
if (isTRUE(anlz_vdj)) {
  ## step 1 | load the data ####
  ## - "scRepertoire functions using the filtered_contig_annotations.csv output from the 10x Genomics Cell Ranger. This file is located in the ./outs/ directory of the VDJ alignment folder. To generate a list of contigs to use for scRepertoire: 
  ## load the filtered_contig_annotations.csv for each of the samples.
  ## make a list in the R environment."
  
  ## directory ##
  dir_res <- paste0(homedir,"result/")
  dir_res_cr <- paste0(dir_res,"cellranger/",prjid,"/")
  
  ## variable ##
  dtype_csv <- "filtered" #"raw" | "filtered"
  list_smpid <- c(unique(mdata$smp_id))
  print(list_smpid)
  
  ## count the number of samples included in the [list_smpid]
  nSmp <- length(list_smpid)
  print(paste0("--- nSmp: ",nSmp," ---"))
  
  ## load the metadata from the [Seurat Object]
  if (nSmp==7) {
    fn_md_mrg_sobj <- paste0(homedir,"collaboration/RCCI__GCAR1/03a_result_GEX/01_Seurat_basic/","metadata_mrg7_sobj_GCAR1_PT01_GEX-subsetVDJ_res2.4",".txt")
  }
  md_mrg_sobj <- read.table(file = fn_md_mrg_sobj, sep = "\t", header = T, row.names = 1)
  dim(md_mrg_sobj) #224380 6 #188417 8(mrg6) #77651 12(GEX-subsetVDJ;nBatch3)
  head(md_mrg_sobj,3)
  length(unique(md_mrg_sobj$barcode)) #198427(mrg; prev) #170235(mrg6) #74305(mrg7;GEX-subsetVDJ;nBatch3)
  
  ## check whether any NA values are included in the [expr_GCAR] column
  subset(md_mrg_sobj, is.na(expr_GCAR)) #0 #checkpoint
  
  ## add a column: [cell_barcode]
  md_mrg_sobj$cell_barcode <- rownames(md_mrg_sobj)
  head(md_mrg_sobj,3)
  dim(md_mrg_sobj) #224380 7(mrg7;prev) #188417 9(mrg6) #77651 12(mrg7;GEX-subsetVDJ;nBatch3)
  

  for (smpid in list_smpid) {
    print(paste0("=== dtype_csv: [",dtype_csv,"] | smpid: [",smpid,"]"," - (",which(list_smpid==smpid)," /",length(list_smpid),") ==="))
    
    ## input ##
    if (dtype_csv=="filtered") {
      dir_csv <- paste0(dir_res_cr,"runs_",sp,"_multi/",smpid,"/outs/per_sample_outs/",smpid,"/vdj_t/")
      fn_contig <- paste0(dir_csv, "filtered_contig_annotations.csv")
      
    } else if (dtype_csv=="raw") {
      dir_csv <- paste0(dir_res_cr,"runs_",sp,"_multi/",smpid,"/outs/multi/vdj_t/")
      fn_contig <- paste0(dir_csv, "all_contig_annotations.csv")
    }
    
    ## read the config file
    contig <- read.csv(file = fn_contig)
    dim(contig) #nVDJ_barcode #nColumn
    head(contig,1)
    print(paste0("--- nTCR_barcode: [",nrow(contig),"] ---"))
    
    ## subset the "source" metadata by each sample
    md_mrg_sobj__smp <- subset(md_mrg_sobj, orig.ident==smpid)
    print(dim(md_mrg_sobj__smp)) #25281 7
    print(head(md_mrg_sobj__smp))
    
    ## add a column: [expr_GCAR]
    contig$expr_GCAR <- md_mrg_sobj__smp$expr_GCAR[match(contig$barcode, md_mrg_sobj__smp$barcode)] #need-to-check-again
    ## check the [expr_GCAR] status
    print(as.data.frame(table(contig$expr_GCAR)))
    # 1  neg 18927
    # 2  pos  5509
    
    # print(head(as.data.frame(table(contig$expr_GCAR, contig$v_gene)),6))
    # 1  neg TRAV1-1   68
    # 2  pos TRAV1-1   23
    # 3  neg TRAV1-2  350
    # 4  pos TRAV1-2  108
    # 5  neg  TRAV10  115
    # 6  pos  TRAV10   31
    # print(head(as.data.frame(table(contig$expr_GCAR, contig$d_gene)),6))
    # 1  neg       17173
    # 2  pos        5005
    # 3  neg TRBD1  1149
    # 4  pos TRBD1   356
    # 5  neg TRBD2   290
    # 6  pos TRBD2   101
    # print(head(as.data.frame(table(contig$expr_GCAR, contig$j_gene)),6))
    # 1  neg TRAJ10  142
    # 2  pos TRAJ10   35
    # 3  neg TRAJ11  106
    # 4  pos TRAJ11   40
    # 5  neg TRAJ12  115
    # 6  pos TRAJ12   48
    # print(head(as.data.frame(table(contig$expr_GCAR, contig$c_gene)),6))
    # 1  neg         79
    # 2  pos         15
    # 3  neg  TRAC 9079
    # 4  pos  TRAC 2687
    # 5  neg TRBC1 2827
    # 6  pos TRBC1  805
    
    ## subset the contig file - based on the [expr_GCAR]
    contig_p <- subset(contig, expr_GCAR=="pos")
    contig_n <- subset(contig, expr_GCAR=="neg")
    
    ## set contig variables for each [smpid]
    if (nSmp==7) {
      if (smpid=="Enriched_Apheresis") {
        contig1_Enriched_Apheresis <- contig
        contig1_Enriched_Apheresis_p <- contig_p #zero
        contig1_Enriched_Apheresis_n <- contig_n
      } else if (smpid=="GCAR1_Product") {
        contig2_GCAR1_Product <- contig
        contig2_GCAR1_Product_p <- contig_p
        contig2_GCAR1_Product_n <- contig_n
      } else if (smpid=="D08") {
        contig3_D08 <- contig
        contig3_D08_p <- contig_p
        contig3_D08_n <- contig_n
      } else if (smpid=="D15") {
        contig4_D15 <- contig
        contig4_D15_p <- contig_p
        contig4_D15_n <- contig_n
      } else if (smpid=="D22") {
        contig5_D22 <- contig
        contig5_D22_p <- contig_p
        contig5_D22_n <- contig_n
      } else if (smpid=="D28") {
        contig6_D28 <- contig
        contig6_D28_p <- contig_p
        contig6_D28_n <- contig_n
      } else if (smpid=="D42") {
        contig7_D42 <- contig
        contig7_D42_p <- contig_p #zero
        contig7_D42_n <- contig_n
      }
      print(":")
    } #closer:  if (nSmp==7) {
  } #closer: for (smpid in list_smpid) {
  
  
  ## checkpoints:
  nrow(subset(contig4_D15, is.na(expr_GCAR))) #364(GEX-subsetVDJ;nBatch3)
  subset(md_mrg_sobj, barcode=="AAACCTGAGCACCGTC-1") #D15 #neg
  
  
  ## make a list of contigs from multiple samples: [contig_list] ####
  ## - note: check the length of each list before making a "merged" contig list
  if (prj=="GCAR1_PT01_VDJ") {
    if (nSmp==7) {
      contig_list <- list(
        contig1_Enriched_Apheresis, contig2_GCAR1_Product,
        contig3_D08, contig4_D15, contig5_D22, contig6_D28, contig7_D42
      )
      contig_list_gcar_pos <- list(
        contig1_Enriched_Apheresis_p,
        contig2_GCAR1_Product_p,
        contig3_D08_p, contig4_D15_p, contig5_D22_p, contig6_D28_p,
        contig7_D42_p
      )
      contig_list_gcar_neg <- list(
        contig1_Enriched_Apheresis_n, contig2_GCAR1_Product_n, 
        contig3_D08_n, contig4_D15_n, contig5_D22_n, contig6_D28_n, contig7_D42_n
      )
    }
    
  } else {
    paste0("[CHECK] contig files")
  } #closer: if (prj=="GCAR1_PT01_VDJ") {
  # ---
  
  
  ## step 2.1 | combine the TCR contigs: [combineTCR] ####
  ## - notes: combining contigs into "clones"
  ## - source: [https://www.borch.dev/uploads/screpertoire/articles/combining_contigs]
  # "There are varying definitions of clones or clones in the literature. 
  # For the purposes of scRepertoire, we will use "clone" and define this as the cells with shared/trackable complementarity-determining region 3 (CDR3) sequences. 
  
  # Within this definition, one might use amino acid (aa) sequences of one or both chains to define a clone. 
  # Alternatively, we could use nucleotide (nt) or the V(D)JC genes (genes) to define a clone. 
  
  # The latter genes would be a more permissive definition of “clones”, as multiple amino acid or nucleotide sequences can result from the same gene combination. 
  
  # Another option to define clone is the use of the V(D)JC and nucleotide sequence (strict). 
  # scRepertoire allows for the use of all these definitions of clones and allows for users to select both or individual chains to examine."
  
  
  fn_rds_combined_tcr <- paste0(dir_res_scR,"combinedTCR_",prj,"_all","_nSmp",nSmp,".rds")
  if (!file.exists(fn_rds_combined_tcr)) {
    gc()
    combined_tcr <- combineTCR(contig_list, 
                               samples = list_smpid, #[recommended] The labels of samples
                               # ID = list_id, #[optional] The additional sample labeling
                               removeNA = F, #FALSE(default) #TRUE(a stringent filter)
                               removeMulti = F, #FALSE(default) #TRUE(a stringent filter)
                               filterMulti = F #FALSE(default) #TRUE(isolated the top 2 expressed chains in cell barcodes with multiple chains)
                               # cells ="T-AB"
    )
    head(combined_tcr[[1]])
    as.data.frame(names(combined_tcr))
    # 1   Enriched_Apheresis
    # 2   GCAR1_Product
    # 3   D08
    # 4   D15
    # 5   D22
    # 6   D28
    # 7   D42
    
    ## add variables to group samples ####
    ## - note: "In addition, we can use the group.by parameter to look at the relative proportion of clones between groups - such as tissue."
    # as.data.frame(names(combined_tcr))
    # 1    Enriched_Apheresis
    # 2    GCAR1_Product
    # 3    D08_afterInfusion
    # 4    D15_afterInfusion
    # 5    D22_afterInfusion
    # 6    D28_afterInfusion
    # 7    D42_afterInfusion
    
    ## [smpOrder] ##
    if (nSmp==7) {
      list_smpOrder <- c(
        paste("1",list_smpid[1],sep = "_"), #Enriched_Apheresis #Enriched
        paste("2",list_smpid[2],sep = "_"), #GCAR1_Product #Harvest
        paste("3",list_smpid[3],sep = "_"), #D08
        paste("4",list_smpid[4],sep = "_"), #D15
        paste("5",list_smpid[5],sep = "_"), #D22
        paste("6",list_smpid[6],sep = "_"), #D28
        paste("7",list_smpid[7],sep = "_") #D42
      )
    } else if (nSmp==6) {
      list_smpOrder <- c(
        paste("1",list_smpid[1],sep = "_"), #Enriched_Apheresis
        paste("2",list_smpid[2],sep = "_"), #GCAR1_Product
        paste("3",list_smpid[3],sep = "_"), #D08
        paste("4",list_smpid[4],sep = "_"), #D15
        paste("5",list_smpid[5],sep = "_"), #D22
        paste("6",list_smpid[6],sep = "_") #D28
      )
    }
    print(list_smpOrder) #checkpoint
    
    combined_tcr <- addVariable(combined_tcr,
                                variable.name = "smpOrder",
                                variables = list_smpOrder
    )
    print(unique(combined_tcr[[1]]$smpOrder)) #checkpoint
    print(unique(combined_tcr[[2]]$smpOrder)) #checkpoint
    print(unique(combined_tcr[[nSmp]]$smpOrder)) #checkpoint
    
    ## [smpType] ##
    combined_tcr <- addVariable(combined_tcr,
                                variable.name = "smpType",
                                variables = c("1_pre","2_product",rep("3_post",(nSmp-2)))
    )
    print(unique(combined_tcr[[1]]$smpType)) #checkpoint
    print(unique(combined_tcr[[2]]$smpType)) #checkpoint
    print(unique(combined_tcr[[nSmp]]$smpType)) #checkpoint
    
    ## [sample_new] ##
    if (nSmp==7) {
      list_smp_new <- c(
        "Enriched_Apheresis","GCAR1_Product",
        "D08","D15","D22","D28","D42"
      )
    } #closer: if (nSmp==7) {
    combined_tcr <- addVariable(combined_tcr,
                                variable.name = "sample_new",
                                variables = list_smp_new
    )
    print(paste0("--- first sample: ",unique(combined_tcr[[1]]$sample_new))) #checkpoint
    print(paste0("--- 2nd sample: ",unique(combined_tcr[[2]]$sample_new))) #checkpoint
    print(paste0("--- last sample: ",unique(combined_tcr[[nSmp]]$sample_new))) #checkpoint
    
    ## save as RDS: [combined_tcr]
    saveRDS(object = combined_tcr, file = fn_rds_combined_tcr)
    print(paste0("[DONE] save as RDS: ",fn_rds_combined_tcr))
    
  } else {
    print(paste0("[SKIP] RDS file exists: ",fn_rds_combined_tcr))
    
    combined_tcr <- readRDS(fn_rds_combined_tcr)
    print(paste0("[DONE] load the RDS files: ",fn_rds_combined_tcr))
  } #closer: if (!file.exists(fn_rds_combined_tcr)) {
  
  
  ## step 2.1.1 | combine the TCR contigs: [combineTCR] - for GCAR-positive & -negative cells ####
  list_gcar_expr <- c("GCARpos","GCARneg")
  
  for (expr in list_gcar_expr) {
    print(paste0("--- expr: [",expr,"] ---"))
    
    if (expr=="GCARpos") {
      contig_list_gcar <- contig_list_gcar_pos
    } else if (expr=="GCARneg") {
      contig_list_gcar <- contig_list_gcar_neg
    }
    list_smpid_gcar <- list_smpid
    nSmp <- length(contig_list_gcar)
    
    fn_rds_combined_tcr_gcar <- paste0(dir_res_scR,"combinedTCR_",prj,"_",expr,"_nSmp",nSmp,".rds")
    if (!file.exists(fn_rds_combined_tcr_gcar)) {
      gc()
      combined_tcr_gcar <- combineTCR(contig_list_gcar, 
                                      samples = list_smpid_gcar , #[recommended] The labels of samples
                                      # ID = list_id, #[optional] The additional sample labeling
                                      removeNA = F, #FALSE(default) #TRUE(a stringent filter)
                                      removeMulti = F, #FALSE(default) #TRUE(a stringent filter)
                                      filterMulti = F #FALSE(default) #TRUE(isolated the top 2 expressed chains in cell barcodes with multiple chains)
                                      # cells ="T-AB"
      )
      print(head(combined_tcr_gcar[[1]]))
      as.data.frame(names(combined_tcr_gcar))
      # 1     1_Enriched_Apheresis
      # 2     2_GCAR1_Product
      # 3     3_D08
      # 4     4_D15
      # 5     5_D22
      # 6     6_D28
      # 7     7_D42
      
      ## add variables to group samples: [smpOrder],[smpType] ####
      if (nSmp==7) {
        list_smpOrder <- c(
          paste("1",list_smpid[1],sep = "_"), #Enriched_Apheresis
          paste("2",list_smpid[2],sep = "_"), #GCAR1_Product
          paste("3",list_smpid[3],sep = "_"), #D08
          paste("4",list_smpid[4],sep = "_"), #D15
          paste("5",list_smpid[5],sep = "_"), #D22
          paste("6",list_smpid[6],sep = "_"), #D28
          paste("7",list_smpid[7],sep = "_") #D42
        )
      } #closer: if (nSmp==7) {
      print(list_smpOrder) #checkpoint
      
      combined_tcr_gcar <- addVariable(combined_tcr_gcar,
                                       variable.name = "smpOrder",
                                       variables = list_smpOrder
      )
      print(unique(combined_tcr_gcar[[1]]$smpOrder)) #checkpoint
      print(unique(combined_tcr_gcar[[2]]$smpOrder)) #checkpoint
      print(unique(combined_tcr_gcar[[nSmp]]$smpOrder)) #checkpoint
      
      ## [smpType] ##
      list_smpType <- c(
        c("1_pre","2_product",rep("3_post",(nSmp-2)))
      )
            
      combined_tcr_gcar <- addVariable(combined_tcr_gcar,
                                       variable.name = "smpType",
                                       variables = list_smpType
      )
      print(unique(combined_tcr_gcar[[1]]$smpType)) #checkpoint
      print(unique(combined_tcr_gcar[[2]]$smpType)) #checkpoint
      print(unique(combined_tcr_gcar[[nSmp]]$smpType)) #checkpoint
      
      ## [sample_new] ##
      if (nSmp==7) {
        list_smp_new <- c(
          "Enriched_Apheresis","GCAR1_Product",
          "D08","D15","D22","D28","D42"
        )
      } #closer: if (nSmp==7) {
      
      
      combined_tcr_gcar <- addVariable(combined_tcr_gcar,
                                       variable.name = "sample_new",
                                       variables = list_smp_new
      )
      print(paste0("--- first sample: ",unique(combined_tcr_gcar[[1]]$sample_new))) #checkpoint
      print(paste0("--- 2nd sample: ",unique(combined_tcr_gcar[[2]]$sample_new))) #checkpoint
      print(paste0("--- last sample: ",unique(combined_tcr_gcar[[nSmp]]$sample_new))) #checkpoint
      
      ## check the VDJ data dimension (for each sample) ####
      for (i in seq(1,nSmp)) {
        print(paste0("[",names(combined_tcr_gcar)[i], "] ",expr,":"))
        print(paste0(dim(combined_tcr_gcar[[i]])))
        print(":")
      }
      
      if (expr=="GCARpos") {
        combined_tcr_gcarP <- combined_tcr_gcar
      } else if (expr=="GCARneg") {
        combined_tcr_gcarN <- combined_tcr_gcar
      }
      
      ## save as RDS: [combined_tcr]
      saveRDS(object = combined_tcr_gcar, file = fn_rds_combined_tcr_gcar)
      print(paste0("[DONE] save as RDS: ",fn_rds_combined_tcr_gcar))
      
    } else {
      print(paste0("[SKIP] RDS file exists: ",fn_rds_combined_tcr_gcar))
      
      if (expr=="GCARpos") {
        fn_rds_combined_tcr <- paste0(dir_res_scR,"combinedTCR_",prj,"_",expr,"_nSmp",nSmp,".rds")
        print(paste0("--- fn_rds_combined_tcr: <",fn_rds_combined_tcr,"> ---"))
        combined_tcr_gcarP <- readRDS(fn_rds_combined_tcr)
      } else if (expr=="GCARneg") {
        fn_rds_combined_tcr <- paste0(dir_res_scR,"combinedTCR_",prj,"_",expr,"_nSmp",nSmp,".rds")
        print(paste0("--- fn_rds_combined_tcr: <",fn_rds_combined_tcr,"> ---"))
        combined_tcr_gcarN <- readRDS(fn_rds_combined_tcr)
      }
      print(paste0("[DONE] load the RDS files: ",fn_rds_combined_tcr_gcar))
    } #closer: if (!file.exists(fn_rds_combined_tcr_gcar)) {
    
  } #closer: for (expr in list_gcar_expr) {
  # ---
  
  # ---
  ## check the VDJ data dimension (for each sample) ####
  ## GCAR-pos VDJ ##
  for (i in seq(1,nSmp)) {
    print(paste0("[",names(combined_tcr_gcarP)[i], "] GCAR-pos:"))
    print(paste0(dim(combined_tcr_gcarP[[i]])))
    print(":")
  }
  # [1] "[Enriched_Apheresis] GCAR-pos:"
  # [1] "6931" "15"  
  # [1] ":"
  # [1] "[GCAR1_Product] GCAR-pos:"
  # [1] "9111" "15"  
  # [1] ":"
  # [1] "[D08] GCAR-pos:"
  # [1] "11082" "15"   
  # [1] ":"
  # [1] "[D15] GCAR-pos:"
  # [1] "11025" "15"   
  # [1] ":"
  # [1] "[D22] GCAR-pos:"
  # [1] "17508" "15"   
  # [1] ":"
  # [1] "[D28] GCAR-pos:"
  # [1] "8648" "15"  
  # [1] ":"
  # [1] "[D42] GCAR-pos:"
  # [1] "17340" "15"   
  # [1] ":"
  
  ## GCAR-neg VDJ ##
  for (i in seq(1,nSmp)) {
    print(paste0("[",names(combined_tcr_gcarN)[i], "] GCAR-neg:"))
    print(paste0(dim(combined_tcr_gcarN[[i]])))
    print(":")
  }
  # [1] "[Enriched_Apheresis] GCAR-neg:"
  # [1] "6931" "15"  
  # [1] ":"
  # [1] "[GCAR1_Product] GCAR-neg:"
  # [1] "9111" "15"  
  # [1] ":"
  # [1] "[D08] GCAR-neg:"
  # [1] "11082" "15"   
  # [1] ":"
  # [1] "[D15] GCAR-neg:"
  # [1] "11025" "15"   
  # [1] ":"
  # [1] "[D22] GCAR-neg:"
  # [1] "17508" "15"   
  # [1] ":"
  # [1] "[D28] GCAR-neg:"
  # [1] "8648" "15"  
  # [1] ":"
  # [1] "[D42] GCAR-neg:"
  # [1] "17340" "15"   
  # [1] ":"
  # ---
  
  
  ## make a variable list
  list_var <- unique(names(combined_tcr))
  print(list_var)
  # ---
  
  
  # ---
  ## make a list of "contig groups"
  list_contigGrp <- c("all","GCARpos","GCARneg")
  
  ## make [plot 1a] for each [cgrp]
  for (cgrp in list_contigGrp) {
    print(paste0("--- contig group: [",cgrp,"] ---"))
    
    if (cgrp=="all") {
      combined_tcr_cgrp <- combined_tcr
    } else if (cgrp=="GCARpos") {
      combined_tcr_cgrp <- combined_tcr_gcarP
    } else if (cgrp=="GCARneg") {
      combined_tcr_cgrp <- combined_tcr_gcarN
    }
    nSmp <- length(combined_tcr_cgrp)
    
    ## step 3.3 | export clones: [exportClones] ####
    ## - most important step!
    if (isTRUE(exp_cln)) {
      fn_csv_clones <- paste0(dir_res_scR,"clones_",prj,"_",cgrp,"__nSmp",nSmp,".csv")
      fname_csv_clones <- paste0("clones_",prj,"_",cgrp,"__nSmp",nSmp,".csv")
      print(fname_csv_clones)
      
      if (!file.exists(fn_csv_clones)) {
        ## save as CSV
        exportClones(combined_tcr_cgrp, 
                     write.file = TRUE,
                     dir = dir_res_scR, 
                     file.name = fname_csv_clones
        )
        print(paste0("[DONE] export clones: ", fn_csv_clones))
      } else {
        print(paste0("[SKIP] CSV file exists: ", fn_csv_clones))
      }
    } #closer: if (isTRUE(exp_cln)) {
  } #closer: for (cgrp in list_contigGrp) {
  # ---
} #closer: if (isTRUE(anlz_vdj)) {
# ---


#### metadata ####
## load the [Scrublet] results
fn_agg_obs_c <- paste0(dir_res_scb,"agg_obs_scb_GCAR1_PT01_GEX_nSmp7.txt")
agg_obs_c <- read.table(fn_agg_obs_c, sep = "\t", header = T, row.names = NULL)
head(agg_obs_c,3)
dim(agg_obs_c) #244112 16(GEX; nBatch3; nSmp7)



# ---
## step 1 | calculate the length of "TRA/TRB chain gene combinations" per cell barcode ####
## load the [cln_all_up]
fn_csv_cln_all_up <- paste0(homedir,"collaboration/","RCCI__GCAR1","/03b_result_VDJ/","clones_GCAR1_PT01_VDJ_all__nSmp",nSmp,"_updated",".csv")
cln_all_up <- read.table(file = fn_csv_cln_all_up, sep = ",", header = T, row.names = NULL, quote = NULL)
head(cln_all_up,3)
dim(cln_all_up) #81645 13->22

## add columns: [len_chain1_aa], [len_chain2_aa]
cln_all_up <- cln_all_up %>%
  group_by(cell_barcode) %>%
  mutate(
    # sapply(strsplit(md_ePath$source__term_name, split = "\\_"), FUN = "[", 1)
    ## [aa] ##
    seqLen_chain1_aa=str_length(chain1_aa),
    seqLen_chain2_aa=str_length(chain2_aa),
    ## [nt] ##
    seqLen_chain1_nt=str_length(chain1_nt),
    seqLen_chain2_nt=str_length(chain2_nt),
  ) %>%
  as.data.frame()
head(cln_all_up,3)
tail(cln_all_up,3)
dim(cln_all_up) #81645 17->23


## add columns: [len_chain1_gene], [len_chain2_gene] - "length of each chain's genes"
cln_all_up$len_chain1_gene <-(str_count(cln_all_up$chain1_genes,";")+1)
cln_all_up$len_chain2_gene <- (str_count(cln_all_up$chain2_genes,";")+1)
head(cln_all_up,3)
## change the "NA" values to "zero" - for visualization
cln_all_up$len_chain1_gene[which(is.na(cln_all_up$len_chain1_gene))] <- 0
cln_all_up$len_chain2_gene[which(is.na(cln_all_up$len_chain2_gene))] <- 0
head(cln_all_up,3)

max(cln_all_up$len_chain1_gene) #2
max(cln_all_up$len_chain2_gene) #2

min(cln_all_up$len_chain1_gene) #0
min(cln_all_up$len_chain2_gene) #0

## add a column: [len_chain1_chain2_gene] - the length of chain1 & chain2 genes
dim(cln_all_up) #81645 20
head(cln_all_up,3)

as.data.frame(table(cln_all_up$combLen_chain1_chain2_gene))
#   Var1  Freq
# 1 0--1 15191
# 2 1--0  2075
# 3 1--1 51712
# 4 1--2  4941
# 5 2--1  5628
# 6 2--2  2098
nrow(subset(cln_all_up, combLen_chain1_chain2_gene=="2--2")) #2098
unique(subset(cln_all_up, combLen_chain1_chain2_gene=="2--2")$chain1_gene) #TRA
unique(subset(cln_all_up, combLen_chain1_chain2_gene=="2--2")$chain2_gene) #TRB
nrow(subset(cln_all_up, combLen_chain1_chain2_gene=="1--2")) #4941 
# ---


# ---
## summarize the [cln_all_up] data frame
sum_cln_all_up <- ddply(cln_all_up,
                        c("sample_new","combLen_chain1_chain2_gene","expr_GCAR"),
                        summarise,
                        nCombLen_chain1_chain2_gene=length(combLen_chain1_chain2_gene)
)
head(sum_cln_all_up,3)
tail(sum_cln_all_up,3)
dim(sum_cln_all_up) #115 4

## set the variables as "factors" - for the visualization
## [sample_new] ##
sum_cln_all_up$sample_new <- factor(sum_cln_all_up$sample_new, 
                                    levels = c(
                                      "Enriched_Apheresis","GCAR1_Product",
                                      "D08","D15","D22","D28","D42"
                                    ))
table(sum_cln_all_up$sample_new)
## [expr_GCAR] ##
sum_cln_all_up$expr_GCAR <- factor(sum_cln_all_up$expr_GCAR, 
                                   levels = c("pos","neg"))
as.data.frame(table(sum_cln_all_up$expr_GCAR))
# 1  pos   37
# 2  neg   42
# ---



#### visualization ####
## colour palette ##
## longitudinal samples ##
col_gcar_p <- "#666666" #dark_grey
col_gcar_t <- "#bf5b17" #brown
col_gcar_d08 <- "#a6d854" #light_green
col_gcar_d15 <- "#e7298a" #pink
col_gcar_d22 <- "#e6ab02" #yellow
col_gcar_d28 <- "#1b9e77" #dark_green
col_gcar_d42 <- "#386cb0" #blue
smpColor <- setNames(c(col_gcar_p,col_gcar_t,col_gcar_d08,col_gcar_d15,col_gcar_d22,col_gcar_d28,col_gcar_d42), c("Enriched_Apheresis","GCAR1_Product","D08","D15","D22","D28","D42")) #sample_new
greyColor_dark <- "#252525"
greyColor_light <- "#969696"


# ---
## step 2 | add a column: [select_combLen_chain1_chain2_gene] - based on the subset ####
# head(cln_all_up,1)
cln_all_up <- cln_all_up %>%
  mutate(
    select_combLen_chain1_chain2_gene=case_when(
      combLen_chain1_chain2_gene=="1--0" | combLen_chain1_chain2_gene=="1--2" | combLen_chain1_chain2_gene=="2--2" ~ "N_exclude",
      combLen_chain1_chain2_gene!="1--0" & combLen_chain1_chain2_gene!="1--2" & combLen_chain1_chain2_gene!="2--2" ~ "Y_select"
    )
  ) %>%
  as.data.frame()
head(cln_all_up,3)
dim(cln_all_up) #81645 21

## save as TXT: [cln_all_up]
fn_csv_cln_all_up <- paste0(dir_res_scR,"clones_GCAR1_PT01_VDJ_all__nSmp",nSmp,"_updated",".csv")
write.csv(x = cln_all_up, file = fn_csv_cln_all_up, col.names = T, row.names = F, quote = F)
print(paste0("[DONE] save as CSV: ",fn_csv_cln_all_up))
# ---


# ---
## step 4 | add columns to each data frame embedded in the "data frame list" ####
## make a list of "contig groups"
list_contigGrp <- c("all","GCARpos","GCARneg")

for (cgrp in list_contigGrp) {
  print(paste0("--- contig group: [",cgrp,"] ---"))
  
  if (cgrp=="all") {
    comb_tcr <- combined_tcr
  } else if (cgrp=="GCARpos") {
    comb_tcr <- combined_tcr_gcarP
  } else if (cgrp=="GCARneg") {
    comb_tcr <- combined_tcr_gcarN
  }
  nSmp <- length(combined_tcr_cgrp)
  
  list_i <- c(seq(1,nSmp,by=1))
  for (i in list_i) {
    ## add an identifier columns: [sample_clones]
    comb_tcr[[i]]$sample_clones <- paste(comb_tcr[[i]]$sample, comb_tcr[[i]]$cdr3_aa1, comb_tcr[[i]]$cdr3_aa2,sep="_")
    
    ## map the columns to the combined TCR data list: [combLen_chain1_chain2_gene], [select_combLen_chain1_chain2_gene]
    comb_tcr[[i]]$combLen_chain1_chain2_gene <- cln_all_up$combLen_chain1_chain2_gene[match(comb_tcr[[i]]$sample_clones, cln_all_up$sample_clones)]
    comb_tcr[[i]]$select_combLen_chain1_chain2_gene <- cln_all_up$select_combLen_chain1_chain2_gene[match(comb_tcr[[i]]$sample_clones, cln_all_up$sample_clones)]
    
    ## map the column: [RNA_snn_res.2.4] = Seurat cell clusters (resolution=2.4)
    comb_tcr[[i]]$RNA_snn_res.2.4 <- cln_all_up$RNA_snn_res.2.4[match(comb_tcr[[i]]$sample_clones, cln_all_up$sample_clones)]
    
    ## subset the updated "data frame list"
    comb_tcr[[i]] <- subset(comb_tcr[[i]], select_combLen_chain1_chain2_gene=="Y_select")
  }
  print(head(comb_tcr[[1]],3)) #checkpoint
  print(tail(comb_tcr[[1]],3)) #checkpoint
  
  print(head(comb_tcr[[7]],3)) #checkpoint
  # subset(comb_tcr[[1]], !is.na(select_combLen_chain1_chain2_gene)) #checkpoint
  
  ## re-assign each "updated data frame list subset" to each [cgrp] variable
  if (cgrp=="all") {
    combined_tcr__sel <- comb_tcr
  } else if (cgrp=="GCARpos") {
    combined_tcr_gcarP__sel <- comb_tcr
  } else if (cgrp=="GCARneg") {
    combined_tcr_gcarN__sel <- comb_tcr
  }
  
} #closer: for (cgrp in list_contigGrp) {
# ---


##############################
#### FINAL plots & tables ####
##############################
# ===
if (isTRUE(mk_plt_final)) {
  ## set input variables
  list_ccall <- "strict"
  chain <- "both"
  
  ## set the number of the top clones of interest
  nTopClone <- 30
  
  ## make a list of "contig groups"
  list_contigGrp <- c("GCARpos","GCARneg")
  
  for (cgrp in list_contigGrp) {
    for (ccall in list_ccall) {
      print(paste0("--- contig group: [",cgrp,"] | cloneCall mode: [",ccall,"] ---"))
      
      if (cgrp=="all") {
        comb_tcr__sel <- combined_tcr__sel
      } else if (cgrp=="GCARpos") {
        comb_tcr__sel <- combined_tcr_gcarP__sel
      } else if (cgrp=="GCARneg") {
        comb_tcr__sel <- combined_tcr_gcarN__sel
      }
      
      ## set [relab_cln] boolean variables - only for [aa] (cf. other [ccall] variables' legend details are too long to included in the plot)
      if (ccall=="aa") {
        relab_cln <- F
      } else {
        relab_cln <- T
      } #closer: if (ccall=="aa") {
      
      ## step 4 | calculate the length of "TRA/TRB chain gene combinations" per cell barcode ####
      ## export table - for each "cloneCall" mode ####
      if (isTRUE(tab_cln__sel)) {
        if (isTRUE(tab_cln__sel__topClones)) {
          ## original clone labels ##
          tab_cln_comp <- clonalCompare(comb_tcr__sel, 
                                        order.by = faclist_smpid,
                                        chain = chain,
                                        top.clones = nTopClone,
                                        cloneCall=ccall, #aa(amino acid) #nt(nucleotide) #strict(VDJC gene + CDR3 nucleotide) #gene(VDJC gene)
                                        relabel.clones = relab_cln, #note: if set to TRUE, clones will be labelled as [Clone: 1] to [Clone: n]
                                        exportTable = T, 
                                        graph = "alluvial", 
                                        palette = pal_basic
          )
          ## save as TXT: [tab_cln_comp]
          fn_txt_cln_comp <- paste0(mydir_res,"clonalCompare_",prjid,"_",cgrp,"_n",nSmp,"__",chain,"_top",nTopClone,"_selClones_",ccall,".txt")
          if (!file.exists(fn_txt_cln_comp)) {
            write.table(x = tab_cln_comp, file = fn_txt_cln_comp, sep = "\t", col.names = T, row.names = F, quote = F)
            print(paste0("[DONE] save as TXT: ",fn_txt_cln_comp))
          } else {
            print(paste0("[SKIP] TXT file exists: ",fn_txt_cln_comp))
          }
          
          ## re-labelled clones ##
          tab_cln_comp_relab <- clonalCompare(comb_tcr__sel, 
                                              chain = chain,
                                              top.clones = nTopClone,  
                                              cloneCall=ccall, #aa(amino acid) #nt(nucleotide) #strict(VDJC gene + CDR3 nucleotide) #gene(VDJC gene)
                                              relabel.clones = relab_cln, #note: if set to TRUE, clones will be labelled as [Clone: 1] to [Clone: n]
                                              exportTable = T, 
                                              graph = "alluvial", 
                                              palette = pal_basic
          )
          ## save as TXT: [tab_cln_comp_relab]
          fn_txt_cln_comp <- paste0(mydir_res,"clonalCompare_",prjid,"_",cgrp,"_n",nSmp,"__",chain,"_top",nTopClone,"_selClones_",ccall,"_reLab",".txt")
          if (!file.exists(fn_txt_cln_comp)) {
            write.table(x = tab_cln_comp_relab, file = fn_txt_cln_comp, sep = "\t", col.names = T, row.names = F, quote = F)
            print(paste0("[DONE] save as TXT: ",fn_txt_cln_comp))
          } else {
            print(paste0("[SKIP] TXT file exists: ",fn_txt_cln_comp))
          }
        } #closer: if (isTRUE(tab_cln__sel__topClones)) {
        
        
        if (isTRUE(tab_cln__sel__allClones)) {
          ## original clone labels - all clones ##
          tab_cln_comp__allClones <- clonalCompare(comb_tcr__sel, 
                                                   order.by = faclist_smpid,
                                                   chain = chain,
                                                   cloneCall=ccall, #aa(amino acid) #nt(nucleotide) #strict(VDJC gene + CDR3 nucleotide) #gene(VDJC gene)
                                                   relabel.clones = F, #note: if set to TRUE, clones will be labelled as [Clone: 1] to [Clone: n]
                                                   exportTable = T, 
                                                   graph = "alluvial", 
                                                   palette = pal_basic
          )
          ## save as TXT: [tab_cln_comp__allClones]
          fn_txt_cln_comp__allClones <- paste0(mydir_res,"clonalCompare_",prjid,"_",cgrp,"_n",nSmp,"__",chain,"_ALL","_selClones_",ccall,".txt")
          if (!file.exists(fn_txt_cln_comp__allClones)) {
            write.table(x = tab_cln_comp__allClones, file = fn_txt_cln_comp__allClones, sep = "\t", col.names = T, row.names = F, quote = F)
            print(paste0("[DONE] save as TXT: ",fn_txt_cln_comp__allClones))
          } else {
            print(paste0("[SKIP] TXT file exists: ",fn_txt_cln_comp__allClones))
          }
          
          if (cgrp=="all") {
            tab_cln_comp__allClones_all <- tab_cln_comp__allClones
          } else if (cgrp=="GCARpos") {
            tab_cln_comp__allClones_GCARpos <- tab_cln_comp__allClones
          } else if (cgrp=="GCARneg") {
            tab_cln_comp__allClones_GCARneg <- tab_cln_comp__allClones
          }
          
        } #closer: if (isTRUE(tab_cln__sel__allClones)) {
      } #closer: if (isTRUE(tab_cln__sel)) {
      
      
      # ---
      ## plot 4 | stacked bar plot - top30 TRB clones from each time point (nSmp6 = Harvest to D42) ####
      ## - figure: [ main figure #4i ]
      if (isTRUE(plt_4_n6)) {
        ## set [chain] & [ccall] variables
        chain <- "both"
        ccall <- "strict"
        
        ## set the detailed note on each [ccall] variable
        if (ccall=="strict") {
          ccall_notes <- "VDJC gene + CDR nucleotide"
        } else if (ccall=="aa") {
          ccall_notes <- "CDR3 amino acid"
        } else if (ccall=="nt") {
          ccall_notes <- "CDR3 nucleotide"
        } else if (ccall=="gene") {
          ccall_notes <- "VDJC gene"
        }
        
        ## plot dimension ##
        mywidth <- 25 #width=18
        myheight <- 8
        
        ## set [relab_cln] boolean variables - only for [aa] 
        ## (cf. other [ccall] variables' legend details are too long to included in the plot)
        if (ccall=="aa") {
          relab_cln <- F
        } else {
          relab_cln <- T
        } #closer: if (ccall=="aa") {
        
        title_cln_comp <- "Comparison of TCR clones btwn samples & changes in dynamics"
        subt_cln_comp <- paste0("patient: [",prjid,"/",cgrp,"]", " | chain: [",chain,"]", " | nClones-of-Interest: [all /top",nTopClone,"]
                        \n- cloneCall mode: ",ccall," (",ccall_notes,")")
        
        fn_pdf_cln_comp <- paste0(mydir_res, "4_clonalCompare_",prjid,"_",cgrp,"_nSmp6__",chain,"_top",nTopClone,"_selClones_",ccall,".pdf")
        pdf(file = fn_pdf_cln_comp, width = mywidth, height = myheight)
        print(paste0("[BEGIN] run 'plt_4_n6': ", Sys.time()))
        print(
          clonalCompare(comb_tcr__sel, 
                        chain = chain,
                        top.clones = nTopClone,  
                        samples = list_var[2:7],
                        cloneCall=ccall, #aa(amino acid) #nt(nucleotide) #strict(VDJC gene + CDR3 nucleotide) #gene(VDJC gene)
                        relabel.clones = relab_cln, #note: if set to TRUE, clones will be labelled as [Clone: 1] to [Clone: n]
                        graph = "alluvial",
                        palette = pal_basic
          ) + 
            theme_bw() +
            theme(
              panel.grid.major.x = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              legend.position = "right"
            ) +
            xlab(label = "") +
            labs(title = title_cln_comp, subtitle = subt_cln_comp) +
            scale_x_discrete(labels = c(list_smp_new)[2:7])
        )
        dev.off()
        print(fn_pdf_cln_comp)
        print(paste0("[DONE] run 'plt_4_n6': ", Sys.time()))
      } #closer:  if (isTRUE(plt_4_n6)) {
      # ---
      
    } #closer: for (ccall in list_ccall) {
  } #closer: for (cgrp in list_contigGrp) {
  
  
  # ---
  ## step 5 | make a summary table of the top30 clones (from each time point) with proportions ####
  if (isTRUE(mk_sum_tab)) {
    ## make a list of "contig groups"
    list_contigGrp <- c("all","GCARpos","GCARneg")
    
    for (cgrp in list_contigGrp) {
      # chain <- "TRB" #"TRB" #"TRA"
      
      list_chain <- c("both") #"TRA" "TRB" "both"
      
      for (chain in list_chain) {
        ## directory ##
        dir_data_byChain <- paste0(homedir,"collaboration/RCCI__GCAR1/03b_result_VDJ/0_top",nTopClone,"_",chain,"_clones__fromEachSmp/")
        
        if (!dir.exists(dir_data_byChain)) {
          dir.create(dir_data_byChain, recursive = T)
          print(paste0("[DONE] create a directory: ",dir_data_byChain))
        } else {
          print(paste0("[SKIP] directory exists: ",dir_data_byChain))
        }
        
        fn_topXX_selCln_comp <- paste0(dir_data_byChain,"clonalCompare_GCAR1_PT01_",cgrp,"_n",nSmp,"__",chain,"_top",nTopClone,"_selClones_",ccall,".txt")
        topXX_selCln_comp <- read.table(file = fn_topXX_selCln_comp, sep = "\t", header = T, row.names = NULL)
        head(topXX_selCln_comp,3)
        dim(topXX_selCln_comp) #582 3
        
        ## add a column: [sample_new]
        topXX_selCln_comp <- topXX_selCln_comp %>%
          mutate(
            sample_new=case_when(
              Sample=="Enriched_Apheresis" ~ "Enriched_Apheresis", #"Apheresis",
              Sample=="GCAR1_Product" ~ "GCAR1_Product", #"GCAR1",
              Sample!="Enriched_Apheresis" & Sample!="GCAR1_Product" ~ Sample
            )
          ) %>%
          as.data.frame()
        head(topXX_selCln_comp,3)
        as.data.frame(table(topXX_selCln_comp$sample_new))
        
        ## set the order of the [sample_new] variables
        topXX_selCln_comp$sample_new <- factor(x = topXX_selCln_comp$sample_new, 
                                               levels = c(
                                                 "Enriched_Apheresis","GCAR1_Product",
                                                 "D08","D15","D22","D28","D42"
                                               ))
        as.data.frame(table(topXX_selCln_comp$sample_new)) #checkpoint
        head(topXX_selCln_comp,3)
        
        
        ## transform the "top30 clones from the selected clone list" 
        # dc_lncnt_goi7 <- dcast(m_lncnt_goi7, formula = cell_id ~ gene, value.var = "logNorm_readCount")
        dc_topXX_selCln_comp <- dcast(topXX_selCln_comp, formula = sample_new ~ clones, value.var = "Proportion")
        head(dc_topXX_selCln_comp[,1:10],3)
        dim(dc_topXX_selCln_comp) #7 121
        ## set row names
        rownames(dc_topXX_selCln_comp) <- dc_topXX_selCln_comp$sample_new
        dc_topXX_selCln_comp <- subset(dc_topXX_selCln_comp, select=-c(sample_new))
        head(dc_topXX_selCln_comp[,1:10],3)
        dim(dc_topXX_selCln_comp) #7 120
        
        ## transpose the "transformed" data frame
        dc_topXX_selCln_comp_tx <- as.data.frame(t(dc_topXX_selCln_comp))
        dim(dc_topXX_selCln_comp_tx) #120 7
        head(dc_topXX_selCln_comp_tx,3)
        
        # ---
        ## replicate the transposed & transformed data frame to check the number of samples with the top30 clones across all the sample set
        chk_topXX_selCln_acrossSmp <- dc_topXX_selCln_comp_tx
        ## replace NA values to "zero" & non-NA values to "one"
        chk_topXX_selCln_acrossSmp[!is.na(chk_topXX_selCln_acrossSmp)] <- 1
        chk_topXX_selCln_acrossSmp[is.na(chk_topXX_selCln_acrossSmp)] <- 0
        head(chk_topXX_selCln_acrossSmp,3)
        ## count the number of samples with the top-ranked clones
        chk_topXX_selCln_acrossSmp$nSmp_of_7 <- rowSums(chk_topXX_selCln_acrossSmp)
        head(chk_topXX_selCln_acrossSmp,3)
        subset(chk_topXX_selCln_acrossSmp, nSmp_of_7!=7)
        
        ## sort by the [nSmp_of_7] column
        chk_topXX_selCln_acrossSmp <- chk_topXX_selCln_acrossSmp %>%
          arrange(-nSmp_of_7) %>%
          as.data.frame()
        head(chk_topXX_selCln_acrossSmp,3)
        tail(chk_topXX_selCln_acrossSmp,3)
        # ---
        
        ## add a column: [dc_topXX_selCln_comp_tx]
        dc_topXX_selCln_comp_tx$nSmp_of_7 <- chk_topXX_selCln_acrossSmp$nSmp_of_7[match(rownames(dc_topXX_selCln_comp_tx), rownames(chk_topXX_selCln_acrossSmp))]
        head(dc_topXX_selCln_comp_tx,3)
        tail(dc_topXX_selCln_comp_tx,3)
        ## sort by the [nSmp_of_7] column
        dc_topXX_selCln_comp_tx <- dc_topXX_selCln_comp_tx %>%
          arrange(-nSmp_of_7) %>%
          as.data.frame()
        head(dc_topXX_selCln_comp_tx,3)
        tail(dc_topXX_selCln_comp_tx,3)
        dim(dc_topXX_selCln_comp_tx) #120 8
        
        ## add row names - depending on the [ccall] variables
        if (ccall=="aa") {
          ## add a column: [chain2_TRB_aa]
          if (chain=="TRB") {
            dc_topXX_selCln_comp_tx$chain2_TRB_aa <- rownames(dc_topXX_selCln_comp_tx)
          } else if (chain=="TRA") {
            dc_topXX_selCln_comp_tx$chain1_TRA_aa <- rownames(dc_topXX_selCln_comp_tx)
          } else if (chain=="both") {
            dc_topXX_selCln_comp_tx$chain_both_aa <- rownames(dc_topXX_selCln_comp_tx)
          }
          head(dc_topXX_selCln_comp_tx,3)
          dim(dc_topXX_selCln_comp_tx) #120 9
          ## change the order of the columns
          # which(colnames(dc_topXX_selCln_comp_tx)=="chain2_TRB_aa") #9
          dc_topXX_selCln_comp_tx <- dc_topXX_selCln_comp_tx[,c(9,1:8)]
          head(dc_topXX_selCln_comp_tx,3)
          dim(dc_topXX_selCln_comp_tx) #120 9
          
        } else if (ccall=="strict") {
          ## add a column: [chain2_TRB_aa]
          if (chain=="TRB") {
            dc_topXX_selCln_comp_tx$chain2_TRB_strict <- rownames(dc_topXX_selCln_comp_tx)
          } else if (chain=="TRA") {
            dc_topXX_selCln_comp_tx$chain1_TRA_strict <- rownames(dc_topXX_selCln_comp_tx)
          } else if (chain=="both") {
            dc_topXX_selCln_comp_tx$chain_both_strict <- rownames(dc_topXX_selCln_comp_tx)
          }
          head(dc_topXX_selCln_comp_tx,3)
          dim(dc_topXX_selCln_comp_tx) #120 9
          ## change the order of the columns
          # which(colnames(dc_topXX_selCln_comp_tx)=="chain2_TRB_aa") #9
          dc_topXX_selCln_comp_tx <- dc_topXX_selCln_comp_tx[,c(9,1:8)]
          head(dc_topXX_selCln_comp_tx,3)
          dim(dc_topXX_selCln_comp_tx) #120 9
        } #closer: if (ccall=="aa") {
        
        
        
        ## save as TXT: [dc_topXX_selCln_comp_tx]
        fn_txt_dc_topXX_selCln_comp_tx <- paste0(dir_data_byChain,"clonalCompare_GCAR1_PT01_",cgrp,"_n7__",chain,"_top",nTopClone,"_selClones_",ccall,"__propPerSmp",".txt")
        if (!file.exists(fn_txt_dc_topXX_selCln_comp_tx)) {
          write.table(x = dc_topXX_selCln_comp_tx, file = fn_txt_dc_topXX_selCln_comp_tx, sep = "\t", col.names = T, row.names = F, quote = F)
          print(paste0("[DONE] save as TXT: ",fn_txt_dc_topXX_selCln_comp_tx))
        } else {
          print(paste0("[SKIP] TXT file exists: ",fn_txt_dc_topXX_selCln_comp_tx))
        }
      } #closer: if (isTRUE(mk_sum_tab)) {
      
    } #closer: for (chian in list_chain) {
  } #closer: for (cgrp in list_contigGrp) {
  # ---

} #closer: if (isTRUE(mk_plt_final)) {
# ===


# ===
## check cell-level cell type prediction (based on SCimilarity) ####

# ---
## step 1 | load the "updated" clone metadata: [cln_all_up] ####
fn_csv_cln_all_up <- paste0(homedir,"collaboration/","RCCI__GCAR1","/03b_result_VDJ/","clones_GCAR1_PT01_VDJ_all__nSmp",nSmp,"_updated",".csv")
cln_all_up <- read.table(file = fn_csv_cln_all_up, sep = ",", header = T, row.names = NULL, quote = NULL)
head(cln_all_up,3)
dim(cln_all_up) #81645 13->22

## step 1-1 | add columns: [chain1_strict], [chain2_strict] ##
## - note: need to add these columns to map with the top clone list (which was extracted using the "strict" mdoe)
cln_all_up <- cln_all_up %>%
  mutate(
    chain1_strict=paste(chain1_genes,chain1_nt,sep=";"),
    chain2_strict=paste(chain2_genes,chain2_nt,sep=";"),
    
  ) %>%
  as.data.frame()
head(cln_all_up,3)
dim(cln_all_up) #81645 24

## step 1-2 | add a column: [TRA_TRB_strict] ##
cln_all_up$TRA_TRB_strict <- paste(cln_all_up$chain1_strict, cln_all_up$chain2_strict, sep = "_")
head(cln_all_up,3)
dim(cln_all_up) #81645 25

## step 1-3 | add an identifier column: [sample__TRA_TRB_strict] ##
cln_all_up$sample__TRA_TRB_strict <- paste(cln_all_up$sample, cln_all_up$TRA_TRB_strict, sep = "__")
head(cln_all_up,3)
dim(cln_all_up) #81645 26

## step 1-4 | update the [cln_all_up] CSV file ##
write.table(x = cln_all_up, file = fn_csv_cln_all_up, sep = ",", col.names = T, row.names = F, quote = F)
print(paste0("[DONE] update the CSV: ",fn_csv_cln_all_up))

## step 1-5 | add cell-level annotation columns ##
## load the <obs.csv> file - generated from running [SCimilarity] ####
fn_obs_c <- paste0(dir_res_scm,"SCim_GCAR1_PT01/mrg7/adata-to-csv_constrain_res2.4_VDJonly","/obs.csv") #constrained prediction
obs_c <- read.csv(fn_obs_c, header = T, row.names = 1)
head(obs_c,3)
dim(obs_c) #77651 19(nBatch3; GEX-subsetVDJ; constrained)
length(unique(rownames(obs_c))) #77651(constrained; nBatch3; GEX-subsetVDJ)

## step 1-5-1 | add columns: [RNA_snn_res.2.4] ##
cln_all_up$RNA_snn_res.2.4 <- obs_c$RNA_snn_res.2.4[match(cln_all_up$cell_barcode, obs_c$cell_barcode)]
head(cln_all_up,3)
# ---

# ---
## step 2 | load the clones-of-interest ####
## option 3 | [top30] clones from [GCAR-neg] cells ##
if (nTopClone==30) {
  fn_top30_selCln_comp__GCARneg_relab <- paste0(homedir,"collaboration/RCCI__GCAR1/03b_result_VDJ/0_top30_both_clones__fromEachSmp/clonalCompare_GCAR1_PT01_GCARneg_n7__both_top30_selClones_strict_reLab.txt")
  top30_selCln_comp__GCARneg_relab <- read.table(file = fn_top30_selCln_comp__GCARneg_relab, sep = "\t", header = T, row.names = NULL, quote = NULL)
  dim(top30_selCln_comp__GCARneg_relab) #602 4(top30; GCAR-neg; both TRA+TRB chains; strict mode)
  head(top30_selCln_comp__GCARneg_relab,3)
  
  ## add an identifier column: [Sample__original.clones]
  top30_selCln_comp__GCARneg_relab$Sample__original.clones <- paste(top30_selCln_comp__GCARneg_relab$Sample, top30_selCln_comp__GCARneg_relab$original.clones, sep="__")
  dim(top30_selCln_comp__GCARneg_relab) #602 5(top30; GCAR-neg; both TRA+TRB chains; strict mode)
  head(top30_selCln_comp__GCARneg_relab,3)
  
  # ---
  ## add cell cluster-level details ####
  ## step 0 | load the cell cluster metadata
  fn_md_clst_res2.4 <- paste0(homedir,"collaboration/RCCI__GCAR1/03a_result_GEX/metadata_annoList_SeuratClusters_GEX-VDJ_res2.4_nBatch3_FINAL.txt")
  md_clst_res2.4 <- read.table(file = fn_md_clst_res2.4, sep = "\t", header = T, row.names = NULL, quote = NULL, check.names = F)
  head(md_clst_res2.4,3)
  
  ## step 1 | add a column: [RNA_snn_res.2.4] - cell cluster
  ## - note: this column includes the Seurat Cell Cluster details (resolution=2.4)
  top30_selCln_comp__GCARneg_relab$RNA_snn_res.2.4 <- cln_all_up$RNA_snn_res.2.4[match(top30_selCln_comp__GCARneg_relab$Sample__original.clones, cln_all_up$sample__TRA_TRB_strict)]
  dim(top30_selCln_comp__GCARneg_relab) #602 6(top30; GCAR-neg; both TRA+TRB chains; strict mode)
  head(top30_selCln_comp__GCARneg_relab,3)
  
  ## step 2 | add columns: [Broad_cell_class],[Specific_cell_class],[Timepoint_range],[GCAR1_status],[Cell_type],[Cell_type__simplified_(w/_Broad_cell_class)]
  ## - details for each cell cluster (based on the "semi-automated" cluster-level annotation)
  top30_selCln_comp__GCARneg_relab$Broad_cell_class <- md_clst_res2.4$Broad_cell_class[match(top30_selCln_comp__GCARneg_relab$RNA_snn_res.2.4, md_clst_res2.4$seurat_clusters_res2.4)]
  top30_selCln_comp__GCARneg_relab$Specific_cell_class <- md_clst_res2.4$Specific_cell_class[match(top30_selCln_comp__GCARneg_relab$RNA_snn_res.2.4, md_clst_res2.4$seurat_clusters_res2.4)]
  top30_selCln_comp__GCARneg_relab$Timepoint_range <- md_clst_res2.4$Timepoint_range[match(top30_selCln_comp__GCARneg_relab$RNA_snn_res.2.4, md_clst_res2.4$seurat_clusters_res2.4)]
  top30_selCln_comp__GCARneg_relab$GCAR1_status <- md_clst_res2.4$GCAR1_status[match(top30_selCln_comp__GCARneg_relab$RNA_snn_res.2.4, md_clst_res2.4$seurat_clusters_res2.4)]
  top30_selCln_comp__GCARneg_relab$Cell_type <- md_clst_res2.4$Cell_type[match(top30_selCln_comp__GCARneg_relab$RNA_snn_res.2.4, md_clst_res2.4$seurat_clusters_res2.4)]
  top30_selCln_comp__GCARneg_relab$`Cell_type__simplified_(w/_Broad_cell_class)` <- md_clst_res2.4$`Cell_type__simplified_(w/_Broad_cell_class)`[match(top30_selCln_comp__GCARneg_relab$RNA_snn_res.2.4, md_clst_res2.4$seurat_clusters_res2.4)]
  dim(top30_selCln_comp__GCARneg_relab) #602 12(top30; GCAR-neg; both TRA+TRB chains; strict mode)
  head(top30_selCln_comp__GCARneg_relab,3)
  
  
  ## save as TXT: [top30_selCln_comp__GCARneg_relab] - w/ cell cluster (CC) details
  fn_top30_selCln_comp__GCARneg_relab_wCC <- paste0(homedir,"collaboration/RCCI__GCAR1/03b_result_VDJ/0_top30_both_clones__fromEachSmp/clonalCompare_GCAR1_PT01_GCARneg_n7__both_top30_selClones_strict_reLab__wCC.txt")
  if (!file.exists(fn_top30_selCln_comp__GCARneg_relab_wCC)) {
    write.table(x = top30_selCln_comp__GCARneg_relab, file = fn_top30_selCln_comp__GCARneg_relab_wCC, sep = "\t", col.names = T, row.names = F)
    print(paste0("[DONE] save as TXT: ",fn_top30_selCln_comp__GCARneg_relab_wCC))
  } else {
    print(paste0("[SKIP] TXT file exists: ",fn_top30_selCln_comp__GCARneg_relab_wCC))
  } #closer: if (!file.exists(fn_top30_selCln_comp__GCARneg_relab_wCC)) {
} #closer: if (nTopClone==30) {

## [top30] GCAR-neg ##
#       clones   Proportion             Sample                                                             original.clones
# 1  Clone: 38 0.0007871537 Enriched_Apheresis NA;NA_TRBV27.TRBD1.TRBJ2-5.TRBC2;TGTGCCAGCACAACAACAGGGGGCCGCGAGACCCAGTACTTC
# 2  Clone: 18 0.0009445844 Enriched_Apheresis       NA;NA_TRBV30.NA.TRBJ1-5.TRBC1;TGTGCCTGGAGTGAACAGGGACGTCAGCCCCAGCATTTT
# 3 Clone: 103 0.0015743073 Enriched_Apheresis    NA;NA_TRBV30.TRBD1.TRBJ2-2.TRBC2;TGTGCCTGGAGTGTACGGTGGGCCCGGGAGCTGTTTTTT
# ---



# ---
## step 3 | load the comparison table of the "top30" clone list extracted using the "strict" mode
## [top30] GCAR-neg ##
if (nTopClone==30) {
  fn_top30_selCln_comp__GCARneg <- paste0(homedir,"collaboration/RCCI__GCAR1/03b_result_VDJ/0_top30_both_clones__fromEachSmp/clonalCompare_GCAR1_PT01_GCARneg_n7__both_top30_selClones_strict__propPerSmp.txt")
  top30_selCln_comp__GCARneg <- read.table(file = fn_top30_selCln_comp__GCARneg, sep = "\t", header = T, row.names = NULL, quote = NULL)
  dim(top30_selCln_comp__GCARneg) #135 9(top30; GCAR-neg; both TRA+TRB chains; strict mode)
  head(top30_selCln_comp__GCARneg,3)
  
  ## add a column: [original.clones]
  top30_selCln_comp__GCARneg$original.clones <- top30_selCln_comp__GCARneg_relab$original.clones[match(top30_selCln_comp__GCARneg$chain_both_strict, top30_selCln_comp__GCARneg_relab$clones)]
  dim(top30_selCln_comp__GCARneg) #135 10(top30; GCAR-neg; both TRA+TRB chains; strict mode)
  head(top30_selCln_comp__GCARneg,3)
} #closer: if (nTopClone==30) {

## [top30] GCAR-neg ##
#   chain_both_strict    Enriched_Apheresis      GCAR1_Product          D08         D15          D22          D28          D42 nSmp_of_7                                                                                                                            original.clones
# 1        Clone: 109           0.001416877       0.0006677796 0.0006242197 0.000179051 7.132668e-05 0.0001379501 0.0004579354         7         TRAV8-6.TRAJ34.TRAC;TGTGCTGCCGACCCGAACACCGACAAGCTCATCTTT_TRBV6-1.NA.TRBJ1-1.TRBC1;TGTGCCAGCAGTGAAAATGGCAGGGCGGGGGACACTGAAGCTTTCTTT
# 2         Clone: 11           0.007084383       0.0006677796 0.0023928423 0.001074306 7.132668e-05 0.0023451511 0.0013738061         7       TRAV8-6.TRAJ8.TRAC;TGTGCTGTGACCTTAGGCGGCTTTCAGAAACTTGTATTT_TRBV4-3.NA.TRBJ1-3.TRBC1;TGCGCCAGCAGCCACCAACCGCCGAGCTCTGGAAACACCATATATTTT
# 3        Clone: 113           0.001259446       0.0003338898 0.0009363296 0.000179051 1.426534e-04 0.0005518002 0.0002616774         7 TRAV14/DV4.TRAJ9.TRAC;TGTGCAATGAGAGAGGGACTGACTGGAGGCTTCAAAACTATCTTT_TRBV5-5.NA.TRBJ1-3.TRBC1;TGTGCCAGCAGCCTAAGGGGACAGTCTGGAAACACCATATATTTT
# ---


# ---
## step 3-1 | add an annotation column: [chk_topXXX_harvestAND]
## - note: this column is initially created from the [top30_selCln_comp__GCARneg] to identify which clonotypes are found to be recurrent at [Harvest] and other time points

## [top30] GCAR-neg ##
if (nTopClone==30) {
  top30_selCln_comp__GCARneg_relab$chk_top30_harvestAND <- top30_selCln_comp__GCARneg$chk_top30_harvestAND[match(top30_selCln_comp__GCARneg_relab$original.clones, top30_selCln_comp__GCARneg$original.clones)]
  head(top30_selCln_comp__GCARneg_relab,3)
  
  subset(top30_selCln_comp__GCARneg_relab, !is.na(chk_top30_harvestAND))
  dim(subset(top30_selCln_comp__GCARneg_relab, !is.na(chk_top30_harvestAND))) #26 5
  subset(top30_selCln_comp__GCARneg_relab, original.clones=="TRAV13-1.TRAJ47.TRAC;TGTGCAGCACCCGGATATGGAAACAAGCTGGTCTTT_TRBV6-5.NA.TRBJ2-7.TRBC2;TGTGCCAGCAGAGCTAGCTACGAGCAGTACTTC") #checkpoint
  subset(top30_selCln_comp__GCARneg_relab, original.clones=="TRAV14/DV4.TRAJ50.TRAC;TGTGCAATGAGAGAGGGCCTTAGGTGGAAAACCTCCTACGACAAGGTGATATTT_TRBV30.NA.TRBJ1-6.TRBC1;TGTGCCCTTAAGGGGAGCTCCTATAATTCACCCCTCCACTTT") #checkpoint
  
  ## step 2-2 | add an identifier column: [Sample__original.clones]
  top30_selCln_comp__GCARneg_relab$Sample__original.clones <- paste(top30_selCln_comp__GCARneg_relab$Sample, top30_selCln_comp__GCARneg_relab$original.clones, sep = "__")
  head(top30_selCln_comp__GCARneg_relab,3)
  
  ## step 2-3 | add a column: [cell_barcode]
  top30_selCln_comp__GCARneg_relab$cell_barcode <- cln_all_up$cell_barcode[match(top30_selCln_comp__GCARneg_relab$Sample__original.clones, cln_all_up$sample__TRA_TRB_strict)]
  head(top30_selCln_comp__GCARneg_relab,3)
  tail(top30_selCln_comp__GCARneg_relab,3)
  subset(top30_selCln_comp__GCARneg_relab, original.clones=="TRAV13-1.TRAJ47.TRAC;TGTGCAGCACCCGGATATGGAAACAAGCTGGTCTTT_TRBV6-5.NA.TRBJ2-7.TRBC2;TGTGCCAGCAGAGCTAGCTACGAGCAGTACTTC") #checkpoint
} #closer: if (nTopClone==30) {
# ---


# ---
## step 4 | load the cell-level cell type prediction - generated from running [SCimilarity] ####
fn_obs_c <- paste0(dir_res_scm,"SCim_GCAR1_PT01/mrg7/adata-to-csv_constrain_res2.4_VDJonly","/obs.csv") #constrained prediction
obs_c <- read.csv(fn_obs_c, header = T, row.names = 1)
head(obs_c,3)
dim(obs_c) #77651 19(nBatch3; GEX-subsetVDJ; constrained)
length(unique(rownames(obs_c))) #77651(constrained; nBatch3; GEX-subsetVDJ)
# ===


# ---
## step 6 | make an alluvial plot for the [Enriched] vs. [Harvest_GCARpos] vs. [Harvest_GCARneg] ####
if (nTopClone==30) {
  ## set the TXT file names to load
  fn_top30_selCln__all <- paste0(dir_res_scR,"2024-11-01_scRepertoire_GCAR1_PT01_VDJ/clonalCompare_GCAR1_PT01_","all","_n7__TRB_top30_selClones_aa.txt")
  fn_top30_selCln__GCARpos <- paste0(dir_res_scR,"2024-11-01_scRepertoire_GCAR1_PT01_VDJ/clonalCompare_GCAR1_PT01_","GCARpos","_n7__TRB_top30_selClones_aa.txt")
  fn_top30_selCln__GCARneg <- paste0(dir_res_scR,"2024-11-01_scRepertoire_GCAR1_PT01_VDJ/clonalCompare_GCAR1_PT01_","GCARneg","_n7__TRB_top30_selClones_aa.txt")
  
  ## load the top30 TRB chain clones (from each time point)
  ## [all] ##
  top30_selCln__all <- read.table(file = fn_top30_selCln__all, sep = "\t", header = T, row.names = NULL, quote = NULL)
  head(top30_selCln__all,3)
  dim(top30_selCln__all) #598 3
  top30_selCln__all$grp_contig <- "all"
  
  ## [GCARpos] ##
  top30_selCln__GCARpos <- read.table(file = fn_top30_selCln__GCARpos, sep = "\t", header = T, row.names = NULL, quote = NULL)
  head(top30_selCln__GCARpos,3)
  dim(top30_selCln__GCARpos) #249 3
  top30_selCln__GCARpos$grp_contig <- "GCARpos"
  
  ## [GCARneg] ##
  top30_selCln__GCARneg <- read.table(file = fn_top30_selCln__GCARneg, sep = "\t", header = T, row.names = NULL, quote = NULL)
  head(top30_selCln__GCARneg,3)
  dim(top30_selCln__GCARneg) #586 3
  top30_selCln__GCARneg$grp_contig <- "GCARneg"
  # ---
  
  ## bind rows of each data frame
  agg_top30_selCln <- bind_rows(top30_selCln__all,top30_selCln__GCARpos,top30_selCln__GCARneg)
  dim(agg_top30_selCln) #1433 4
  head(agg_top30_selCln,3)
  ## add a column: [sample_new]
  agg_top30_selCln <- agg_top30_selCln %>%
    mutate(
      sample_new=case_when(
        Sample=="Enriched_Apheresis" ~ "Enriched_Apheresis", #"Apheresis",
        Sample=="GCAR1_Product" ~ "GCAR1_Product", #"GCAR1",
        Sample!="Enriched_Apheresis" & Sample!="GCAR1_Product" ~ Sample
      )
    ) %>%
    as.data.frame()
  head(agg_top30_selCln,3)
  
  ## add a column: [selected_set]
  agg_top30_selCln <- agg_top30_selCln %>%
    mutate(
      selected_set=case_when(
        sample_new=="Enriched_Apheresis" & grp_contig=="all" ~ "Y_selected",
        sample_new=="GCAR1_Product" & grp_contig=="GCARpos" ~ "Y_selected",
        sample_new=="GCAR1_Product" & grp_contig=="GCARneg" ~ "Y_selected"
      ) 
    )
  head(agg_top30_selCln,3)
  dim(agg_top30_selCln) #1433 6
  table(agg_top30_selCln$selected_set)
  
  ## subset the aggregated data frame
  agg_top30_selCln__E_Hpn <- subset(agg_top30_selCln, selected_set=="Y_selected")
  dim(agg_top30_selCln__E_Hpn) #142 6
  head(agg_top30_selCln__E_Hpn,3)
  ## add a column: [grp_sample_new]
  agg_top30_selCln__E_Hpn$grp_sample_new <- paste(agg_top30_selCln__E_Hpn$sample_new, agg_top30_selCln__E_Hpn$grp_contig, sep = "_")
  head(agg_top30_selCln__E_Hpn,3)
  ## set the order of the [grp_sample_new]
  agg_top30_selCln__E_Hpn$grp_sample_new <- factor(agg_top30_selCln__E_Hpn$grp_sample_new,
                                                   levels = c("Harvest_GCARpos","Enriched_all","Harvest_GCARneg"))
  table(agg_top30_selCln__E_Hpn$grp_sample_new)
  
  ## save as TXT: [agg_top30_selCln__E_Hpn]
  fn_txt_agg_top30_selCln__E_Hpn <- paste0(mydir_res,"top30_TRB_clones_",prj,"_ALL+GCARPOS+NEG",".txt")
  if (!file.exists(fn_txt_agg_top30_selCln__E_Hpn)) {
    write.table(x = agg_top30_selCln__E_Hpn, file = fn_txt_agg_top30_selCln__E_Hpn, sep = "\t", col.names = T, row.names = F, quote = F)
    print(paste0("[DONE] save as TXT: ",fn_txt_agg_top30_selCln__E_Hpn))
  }
  # ---
  
  ## step 7 | aggregate all clones from the selected TRB chain clone list ####
  ## add a column: [grp_contig]
  tab_cln_comp__allClones_all$grp_contig <- "all"
  tab_cln_comp__allClones_GCARpos$grp_contig <- "GCARpos"
  tab_cln_comp__allClones_GCARneg$grp_contig <- "GCARneg"
  
  ## bind rows of each data frame
  agg_allClones_selCln <- bind_rows(tab_cln_comp__allClones_all, tab_cln_comp__allClones_GCARpos, tab_cln_comp__allClones_GCARneg)
  dim(agg_allClones_selCln) #95010 4
  head(agg_allClones_selCln,3)
  subset(agg_allClones_selCln, clones=="CASRASYEQYF")
} #closer: if (nTopClone==30) {
# ---



#### visualization ####
## colour palette ##
## TRB clone ##
list_clones_TRB <- sort(unique(agg_top30_selCln__E_Hpn$clones))
nCloneTRB <- length(list_clones_TRB)
clnColor <- setNames( c(kelly(n = 22),alphabet2(n = 26),cols25(n = 25),glasbey(n = 32),okabe(n = 8))[1:length(list_clones_TRB)], list_clones_TRB ) #108

# ---
## plot 3 | stacked bar plots w/ alluvial ####
## titles ##
mytitle <- paste0("[",prj,"] Top30 TRB chain clones")
mysubt <- paste0("nCloneTRB: ",nCloneTRB)
## labels ##
mylab_x <- ""
mylab_y <- "Proportion of clones"
lgd_fill <- "TRB clone"


ggsbar_EHpn <- 
  ggplot(data = agg_top30_selCln__E_Hpn, aes(x = grp_sample_new, y=Proportion, fill=clones)) +
  theme_bw() +
  geom_bar(
    stat = "identity",
    width = .25
  ) +
  ## white border around the bar plots
  geom_flow(aes(alluvium=clones), curve_type = "linear", alpha=0.25, color="white") +
  geom_col(width = .55, color = "white") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "right", #"right" #"bottom"
    strip.background = element_rect(fill = "transparent")
  ) +
  scale_fill_manual(values = clnColor) +
  scale_y_continuous(expand = c(0.001,0.001)) + #to start from zero (w/o additional space on the beginning of the y-axis)
  ylab(label = mylab_y) +
  xlab(label = mylab_x) +
  labs(
    title = mytitle,
    subtitle = mysubt,
    fill = lgd_fill
  )

## save as PDF: [ggsbar_EHpn]
fn_pdf_ggsbar_EHpn <- paste0(mydir_res,"3_sBar_cloneTRB_",prj,"_nSmp",length(unique(agg_top30_selCln__E_Hpn$sample_new)),"_nClone",nCloneTRB,".pdf")
ggsave(filename = fn_pdf_ggsbar_EHpn, plot = ggsbar_EHpn, device = "pdf", width = 20, height = 16, units = "in", dpi = 300)
print(paste0("[DONE] save as PDF: ",fn_pdf_ggsbar_EHpn))
# ---
# ===


## load the metadata from the [Seurat Object]
fn_md_mrg_sobj <- paste0(homedir,"collaboration/RCCI__GCAR1/03a_result_GEX/01_Seurat_basic/","metadata_mrg7_sobj_GCAR1_PT01_GEX-subsetVDJ_res2.4",".txt")
md_mrg_sobj <- read.table(file = fn_md_mrg_sobj, sep = "\t", header = T, row.names = 1)
dim(md_mrg_sobj) #224380 6 #188417 8(mrg6) #77651 12(GEX-subsetVDJ;nBatch3)
head(md_mrg_sobj,3)
length(unique(md_mrg_sobj$barcode)) #198427(mrg; prev) #170235(mrg6) #74305(mrg7;GEX-subsetVDJ;nBatch3)

## check whether any NA values are included in the [expr_GCAR] column
subset(md_mrg_sobj, is.na(expr_GCAR)) #0 #checkpoint
# ---


#### violin-box plot ####
if (isTRUE(final_plt_vbox)) {
  library(ggpubr)  
  # ---
  ## 
  if (isTRUE(plt_vbox_vdj)) {    
    ## step 1 | load the "selected" top30 clones (from both TRA+TRB chains of the TRB-filtered VDJ+ cells) ####
    ## - notes: the conditions used to extract these top30 clones
    ##          #1. TRB-filtered VDJ+ cells
    ##          #2. both chains (TRA + TRB)
    ##          #3. strict mode (VDJC gene + CDR nucleotide)
    
    fn_top30_selCln_both_fltTRB_GCARpos <- paste0(homedir,"collaboration/RCCI__GCAR1/","03b_result_VDJ","/0_top30_both_clones__fromEachSmp__TRB-filtered/","clonalCompare_GCAR1_PT01_GCARpos_n7__both_top30_selClones_strict.txt")
    fn_top30_selCln_both_fltTRB_GCARneg <- paste0(homedir,"collaboration/RCCI__GCAR1/","03b_result_VDJ","/0_top30_both_clones__fromEachSmp__TRB-filtered/","clonalCompare_GCAR1_PT01_GCARneg_n7__both_top30_selClones_strict.txt")
    
    ## [GCAR1-pos] ##
    top30_selCln_both_fltTRB_GCARpos <- read.table(file = fn_top30_selCln_both_fltTRB_GCARpos, sep = "\t", header = T, row.names = NULL, check.names = F)
    dim(top30_selCln_both_fltTRB_GCARpos) #233 4
    
    ## [GCAR1-neg] ##
    top30_selCln_both_fltTRB_GCARneg <- read.table(file = fn_top30_selCln_both_fltTRB_GCARneg, sep = "\t", header = T, row.names = NULL, check.names = F)
    dim(top30_selCln_both_fltTRB_GCARneg) #602 4
    # ---
    
    ## step 2 | extract clones identified at [D15] ####
    top30_selCln_both_fltTRB_GCARpos__D15 <- subset(top30_selCln_both_fltTRB_GCARpos, Sample=="D15")
    top30_selCln_both_fltTRB_GCARneg__D15 <- subset(top30_selCln_both_fltTRB_GCARneg, Sample=="D15")
    print(paste0("--- top30 clones at [D15]: ",nrow(top30_selCln_both_fltTRB_GCARpos__D15)," (GCARpos) | ",nrow(top30_selCln_both_fltTRB_GCARneg__D15), " (GCARneg) ---"))
    
    ## step 2a | add an identifier column - before binding rows of each data frame ##
    top30_selCln_both_fltTRB_GCARpos__D15$expr_GCAR <- "pos"
    head(top30_selCln_both_fltTRB_GCARpos__D15,3)
    top30_selCln_both_fltTRB_GCARneg__D15$expr_GCAR <- "neg"
    head(top30_selCln_both_fltTRB_GCARneg__D15,3)
    
    ## step 2b | modify the order of the columns ##
    top30_selCln_both_fltTRB_GCARpos__D15 <- subset(top30_selCln_both_fltTRB_GCARpos__D15, select=c(expr_GCAR,clones,Proportion,Sample,original.clones))
    head(top30_selCln_both_fltTRB_GCARpos__D15,3)
    top30_selCln_both_fltTRB_GCARneg__D15 <- subset(top30_selCln_both_fltTRB_GCARneg__D15, select=c(expr_GCAR,clones,Proportion,Sample,original.clones))
    head(top30_selCln_both_fltTRB_GCARneg__D15,3)
    
    ## step 3 | bind rows of each data frame ##
    top30_selCln_both_fltTRB_exprGCAR__D15 <- bind_rows(top30_selCln_both_fltTRB_GCARpos__D15,top30_selCln_both_fltTRB_GCARneg__D15)
    head(top30_selCln_both_fltTRB_exprGCAR__D15,3)
    tail(top30_selCln_both_fltTRB_exprGCAR__D15,3)
    dim(top30_selCln_both_fltTRB_exprGCAR__D15) #153 5
    
    ## step 3a | make a mapping table for the "shared" T cell clones ##
    comp_top30_selCln_both_fltTRB_exprGCAR__D15 <- data.frame(clones_top30=unique(top30_selCln_both_fltTRB_exprGCAR__D15$original.clones))
    dim(comp_top30_selCln_both_fltTRB_exprGCAR__D15) #135 1
    
    ## step 3b | add a column: [shared_clone_id]
    ## - this step is required to make a box plot with lines joining paired points
    comp_top30_selCln_both_fltTRB_exprGCAR__D15$shared_clone_id <- seq(1,nrow(comp_top30_selCln_both_fltTRB_exprGCAR__D15),by=1)
    head(comp_top30_selCln_both_fltTRB_exprGCAR__D15,3)
    tail(comp_top30_selCln_both_fltTRB_exprGCAR__D15,3)
    dim(comp_top30_selCln_both_fltTRB_exprGCAR__D15) #135 2
    
    ## step 3c | add the column to the existing merged data frame: [top30_selCln_both_fltTRB_exprGCAR__D15] ##
    top30_selCln_both_fltTRB_exprGCAR__D15$shared_clone_id <- comp_top30_selCln_both_fltTRB_exprGCAR__D15$shared_clone_id[match(top30_selCln_both_fltTRB_exprGCAR__D15$original.clones, comp_top30_selCln_both_fltTRB_exprGCAR__D15$clones_top30)]
    head(top30_selCln_both_fltTRB_exprGCAR__D15,3)
    tail(top30_selCln_both_fltTRB_exprGCAR__D15,3)
    table(top30_selCln_both_fltTRB_exprGCAR__D15$shared_clone_id) #check whether any shared T cell clones are observed
    
    
    ## step 3d | fix the "order" of variable in the [variable]
    as.data.frame(table(top30_selCln_both_fltTRB_exprGCAR__D15$expr_GCAR))
    # 1  neg   97
    # 2  pos   56
    top30_selCln_both_fltTRB_exprGCAR__D15$expr_GCAR <- factor(top30_selCln_both_fltTRB_exprGCAR__D15$expr_GCAR,
                                                               levels = c("pos","neg"))
    as.data.frame(table(top30_selCln_both_fltTRB_exprGCAR__D15$expr_GCAR))
    head(top30_selCln_both_fltTRB_exprGCAR__D15,3)
    tail(top30_selCln_both_fltTRB_exprGCAR__D15,3)
    
    # ---
    ## step 4 | violin-box plot - with lines joining paired points (main Figure 4j): [ggvbox_pair] ####
    ## - source: [https://www.geeksforgeeks.org/how-to-connect-data-points-on-boxplot-with-lines-in-r/]
    mytitle <- paste0("[",prj,"] Shared T cell clonotypes")
    mysubt <- paste0("- using 'cell-level' GCAR-pos/neg labels -- [pos = nCount >0]\n- conditions:\n1] TRB-filtered;\n2] both (TRA+TRB);\n3] strict (VDJC gene + CDR nucleotide)")
    mylab_x <- ""
    mylab_y <- "Clonotype Proportion at D15 (both TRA+TRB)"
    mylab_fill <- "GCAR1 expr"
    
    ggvbox_pair <- 
      ggplot(data = top30_selCln_both_fltTRB_exprGCAR__D15, aes(x=expr_GCAR, y=Proportion, fill=expr_GCAR)) +
      theme_bw() +
      geom_violin(alpha=.1, color=greyColor_light) +
      geom_boxplot(width=.2, alpha=.4) + 
      geom_line(aes(group=shared_clone_id), colour=col_gcar_d15) +
      scale_fill_manual(values = exprColor) +
      labs(
        title = mytitle, subtitle = mysubt
      ) +
      stat_compare_means(na.rm = T, label.sep = "\n", label.x = 1.35) + #label.y=.0058,
      scale_fill_manual(values = exprColor) +
      xlab(label = mylab_x) +
      ylab(label = mylab_y) +
      labs(fill=mylab_fill)
    
    ## save as PDF: [ggvbox_pair]
    mywidth <- 5.5
    myheight <- 8
    fn_pdf_ggvbox <- paste0(mydir_res,"4j_vbox-pair_propVDJ-T__",prj,"_nSmp1-onlyD15","_byGCARexpr",".pdf")
    ggsave(filename = fn_pdf_ggvbox, plot = ggvbox_pair, device = "pdf", width = mywidth, height = myheight, units = "in", dpi = 300)
    print(paste0("[DONE] save as PDF: ",fn_pdf_ggvbox))
  } #closer: if (isTRUE(plt_vbox_vdj)) {
  # ---
  
} #closer: if (isTRUE(final_plt_vbox)) {
# ---
# ===