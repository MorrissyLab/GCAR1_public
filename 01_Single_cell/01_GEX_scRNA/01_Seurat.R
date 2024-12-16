######################################
#### single cell RNA-seq analysis ####
#### - project: [GCAR1]           ####
######################################
# version: 2.23.5 (2024-12-15)
# by hyojinsong
today <- Sys.Date()
# ---

#### packages ####
mypackages <- c(
  "BiocManager", "ggplot2", "plyr", "dplyr", "reshape2", "RColorBrewer", "pals" # "readr", "parallel",
  ,"Seurat", "hdf5r", "sctransform" #, "glmGamPoi" #,"nlme" #for Seurat
  ,"ggalluvial"
  # ,"reticulate"
)
# "SeuratDisk",
# "cellranger" 
# "biomaRt"
# "ggside", "ggpubr", "ggtree", "aplot", #for complex plots
# "scater", "scran", "qpdf", 
# "tidyr"


# for (pkg in mypackages) {
#   if (pkg %in% rownames(installed.packages()==F)) {
#     install.packages(pkg)
#   }
# }

## [glmGamPoi] ##
# install.packages('BiocManager')
# BiocManager::install('glmGamPoi')

## [anndata] ##
# install.packages("anndata", force=T, dependencies=T, quiet=F, repos="https://mirror.csclub.uwaterloo.ca/CRAN/")
## [cli] ##
# install.packages("cli", force=T, dependencies=T, quiet=F, repos="https://mirror.csclub.uwaterloo.ca/CRAN/")
## [reticulate]
# install.packages("reticulate", force=T, dependencies=T, quiet=F, repos="https://mirror.csclub.uwaterloo.ca/CRAN/")

## [SeuratDisk] ##
# if (!requireNamespace("remotes", quietly = TRUE)) {
#   install.packages("remotes")
# }
# remotes::install_github("mojaveazure/seurat-disk")

## [hdf5r] ##
# install.packages("hdf5r")

## [pals] ##
# install.packages("pals", force=T, dependencies=T, quiet=F, repos="https://mirror.csclub.uwaterloo.ca/CRAN/")


# install.packages("nlme",dependencies = T,force=T)
# use_condaenv('seurat5') #template

lapply(mypackages, library, character.only = TRUE)
# ---


#### running switches ####
save_wspace <- F
load_wspace <- F

mk_mdata <- F #nSmp=7 (both GEX & VDJ)

## [Seurat] ##
run_seurat <- F #T #[master-running-switch] #important to turn ON when updating the resolution!
tcell_only <- T #F(for initial runs) #T(GEX-subsetVDJ) #[master-running-switch] #either all(GEX) or subset(VDJ) - but should run "all" for the initial run

sr_mode <- "mrg" #each | mrg
sr_dtype <- "raw" #raw | filtered
if (isTRUE(tcell_only)) {
  sr_mode <- paste(sr_dtype,"vdj",sep="_") #each | mrg(7)
} #closer: if (isTRUE(tcell_only)) {

## step 1: [sobj]
sr_01_mk_sobj <- F #T
sr_01a_up_sobj <- T #update the [sobj] metadata w/ GCAR-expressing cells
sr_mk_plt_qc <- F
sr_mk_plt_1 <- F
sr_mk_plt_1ab <- F


## step 2: [sobj_flt]
sr_02_flt <- T #T
## step 3: [sobj_norm]
sr_03_norm <- F #T
## step 4: [feature selection]
sr_04_fsel <- F #T
sr_mk_plt_2 <- F
## step 5: [scaling]
sr_05_scale <- F #T
## step 6: [SCTransform]
sr_06_sct <- F #skip

## step 7: [PCA]
sr_07_pca <- F #T
sr_mk_plt_3 <- F
## step 8: [clustering]
sr_08_clst <- F #T
res_clst <- 2.4 #0.8(default) #2.4(selected0

## step 9: [UMAP]
sr_09_umap <- F #T #important to turn OFF if not needed!
sr_mk_plt_4 <- F #UMAP plots
sr_09a_afterUMAP <- F
sr_mk_plt_5 <- F #UMAP-based GOI plots
## save workspace - after UMAP step
sr_save_wspace <- F #set TRUE if you want to keep the track

## step 10: [DE]
sr_10_de <- F
sr_10_de_FindAllMarkers <- F #T (both=pos+neg)
sr_10_de_FindMarkers <- F
sr_10_de_vis <- F

## step 11: [cell type assignment]
sr_11_ctype <- F
# ---

## merge Seurat Objects
sr_mrg_sobj <- F #T #mrg(7)

ext_sobj_mdata <- F #T #[master-running-switch]
ext_sobj_mdata_flt <- F #T

add_sobj_mdata <- F #[master-running-switch]
add_sobj_mdata_vdj <- F
# add_sobj_mdata_covid <- F #done
add_sobj_mdata_trb <- T


## convert the SeuratObject to h5ad
sceasy_cnvt_to_h5ad <- T #[master-running-switch] #separate for-loop #better to run separately from the Seurat section due to the "out-of-memory" error


if (!isTRUE(sceasy_cnvt_to_h5ad)) {
  lapply(mypackages, library, character.only = TRUE)
}
# ---

#### project ####
if (!isTRUE(tcell_only)) {
  prj <- "GCAR1_PT01_GEX"
} else {
  prj <- "GCAR1_PT01_GEX-subsetVDJ"
}

if (prj=="GCAR1_PT01_GEX_VDJ" | prj=="GCAR1_PT01_GEX" | prj=="GCAR1_PT01_VDJ" | prj=="GCAR1_PT01_GEX-subsetVDJ") {
  prjid <- "GCAR1_PT01"
  sp <- "hs38_GCAR"
}


#### directory ####
homedir <- "/Users/your_username/"
dir_data <- paste0(homedir,"data/",prjid,"/")
dir_res <- paste0(homedir,"result/")
dir_res_cr <- paste0(homedir,"result/cellranger/",prjid,"/")
## [Seurat] ##
dir_res_sr <- paste0(homedir,"result/Seurat/")
if (!dir.exists(dir_res_sr)) {
  dir.create(dir_res_sr, recursive = T)
  print(paste0("[DONE] create a directory: ",dir_res_sr))
} else {
  print(paste0("[SKIP] directory exists: ",dir_res_sr))
}
setwd(dir = dir_res_sr)


#### metadata ####
if (isTRUE(mk_mdata)) {
  ## initialize metadata
  mdata <- data.frame(
    smp_id=c("Enriched","Harvest","D08","D15","D22","D28","D42")
  )
  ## add columns:
  mdata <- mdata %>%
    mutate(
      smp_type=case_when(
        smp_id=="Enriched" ~ "Enriched", #"untransducedCAR"
        smp_id=="Harvest" ~ "Harvest", #"transduced GCAR1"
        smp_id=="D08"|smp_id=="D15"|smp_id=="D22"|smp_id=="D28"|smp_id=="D42" ~ "afterInfusion"
      ),
      smp_name=case_when(
        smp_id=="Enriched" ~ "CD4_CD8_enriched_product_untransducedCARs",
        smp_id=="Harvest" ~ "GCAR1_product_beforeInfusion",
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
  ## load the metadata - skip this step if [run_seurat] running switch is turned on
  if (!isTRUE(run_seurat)) {
    fn_mdata <- paste0(dir_data,"metadata_GCAR1_PT01_n7.txt")
    if (file.exists(fn_mdata)) {
      mdata <- read.table(file = fn_mdata, sep = "\t", header = T, row.names = NULL, quote = NULL)
      print(paste0("[DONE] load the metadata: ", fn_mdata))
    } else {
      print(paste0("[CHECK] metadata does not exist: ", fn_mdata))
    }
  } #closer: if (!isTRUE(run_seurat)) {
  
} #closer: if (isTRUE(mk_mdata)) {


################
#### Seurat ####
################
if (isTRUE(run_seurat)) {
  #### project ####
  ## set [prj] - based on "tcell_only" running switch ####
  prj_gex <- "GCAR1_PT01_GEX"
  if (!isTRUE(tcell_only)) {
    prj <- prj_gex
  } else {
    prj <- "GCAR1_PT01_GEX-subsetVDJ"
  }
  
  ## set dependent variables
  if (prj=="GCAR1_PT01_GEX_VDJ" | prj=="GCAR1_PT01_GEX") {
    prjid <- "GCAR1_PT01"
    sp <- "hs38_GCAR"
  }
  print(paste0("=== prj: [",prj,"] | prjid: [",prjid,"] | species: [",sp,"] ==="))
  
  #### metadata ####
  fn_mdata <- paste0(dir_data,"metadata_GCAR1_PT01_n7.txt")
  if (file.exists(fn_mdata)) {
    mdata <- read.table(file = fn_mdata, sep = "\t", header = T, row.names = NULL, quote = NULL)
    print(paste0("[DONE] load the metadata: ", fn_mdata))
    
  } else {
    print(paste0("[CHECK] metadata does not exist: ", fn_mdata))
  }
  
  ## define a list of sample IDs
  if (sr_mode=="each") {
    list_smpid <- unique(mdata$smp_id)
    list_smpid <- list_smpid
  } else if (sr_mode=="mrg") {
    list_smpid <- c("mrg7") #concatenated_nBatch3
  }
  
  
  #### directory ####
  dir_res_sr <- paste0(homedir,"result/Seurat/")
  if (!dir.exists(dir_res_sr)) {
    dir.create(dir_res_sr, recursive = T)
    print(paste0("[DONE] create a directory: ",dir_res_sr))
  }
  
  
  ## set threshold ####
  minCell <- 1 #3(default)
  minFeat <- 100 #200(default)
  
  # ---  
  for (smpid in list_smpid) {
    print(paste0("=== smp_id: [",smpid,"] ==="))
    
    ## input directory ##
    if (sr_mode=="each") {
      dir_df <- paste0(dir_res_cr, "runs_", sp, "_multi/")
    } #closer: if (sr_mode=="each") {
    
    ## output directory ##
    mydir_res <- paste0(dir_res_sr,today,"_Seurat","_",prj,"/")
    mydir_smp <- paste0(mydir_res,smpid,"/")
    if (!dir.exists(mydir_smp)) {
      dir.create(mydir_smp, recursive = T)
      print(paste0("[DONE] create a directory: ",mydir_smp))
    } else {
      print(paste0("[SKIP] directory exists: ",mydir_smp))
    }
    
    ## step 0 | load the 10x Cell Ranger output files (as an input for Seurat) ####
    if (isTRUE(sr_01_mk_sobj)) {
      print(paste0("--- STEP 0 | load the 10x Cell Ranger output files | smp_id: [",smpid,"]"," (",which(list_smpid==smpid)," /",length(list_smpid),")"," ---"))
      
      if (sr_mode=="each") {
        
        if (sr_dtype=="raw") {
          fn_h5_10x <- paste0(dir_df, smpid, "/outs/multi/count/","raw_feature_bc_matrix.h5")
          dir_df_smp <- paste0(dir_df, smpid, "/outs/multi/count/","raw_feature_bc_matrix/")
          
        } else if (sr_dtype=="filtered") {
          fn_h5_10x <- paste0(dir_res_cr,prjid,"/runs_",sp,"_multi/", smpid, "/outs/multi/count/","sample_filtered_feature_bc_matrix.h5")
          # h5_10x <- Read10X_h5(fn_h5_10x)
          dir_df_smp <- paste0(dir_df, smpid, "/outs/multi/count/","sample_filtered_feature_bc_matrix/")
        } #closer: if (sr_dtype=="raw") {
        print(dir_df_smp)
        
        ## set the option to make the SeuratObject in "v3" format
        ## - this function has to run before loading the "Seurat" package
        ver_seurat <- "v3"
        options(Seurat.object.assay.version = ver_seurat)
        print(paste0("[DONE] set Seurat format version to: [",ver_seurat,"]"))
        
        ## read 10x data 
        gex <- Read10X(data.dir = dir_df_smp)
        
        ## step 1 | create a SeuratObject ####
        ## - [min.cells]: Include features detected in at least this many cells. Will subset the counts matrix as well. To reintroduce excluded features, create a new object with a lower cutoff
        ## - [min.features]: Include cells where at least this many features are detected
        print(paste0("--- STEP 1 | create a SeuratObject | smp_id: [",smpid,"]"," (",which(list_smpid==smpid)," /",length(list_smpid),")"," ---"))
        
        # minCell <- 1 #3(default)
        # minFeat <- 100 #200(default)
        
        if (sr_dtype=="raw"|sr_dtype=="filtered") {
          sobj <- CreateSeuratObject(counts = gex, project = smpid, min.cells = minCell, min.features = minFeat)
        }
        
        ## check the Seurat format version
        print(sobj@version)
        # print(class(sobj[["RNA"]]))
        
        
        ## define the format of the mitochondrial genes ####
        # list_patt_mt <- c(
        #   "^MT-", #for [hs38]
        #   "^mt-", #for [mm10]
        #   "^Mt-",
        #   "^GRCh38-MT-",
        #   "^GRCh38-mt-"
        # )
        patt_mt <- "^MT-"
        PercentageFeatureSet(sobj, pattern = patt_mt) %>% head() #checkpoint
        ## add a column: [percent.mt]
        sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = patt_mt)
        
        ## save as RDS: [sobj]
        fn_rds_sobj <- paste0(dir_res_sr,"sobj_",prj,"_",smpid,"_",sr_dtype,".rds")
        if (!file.exists(fn_rds_sobj)) {
          saveRDS(object = sobj, file = fn_rds_sobj)
          print(paste0("[DONE] save as RDS: ",fn_rds_sobj))
        }
      } #closer: if (sr_mode=="each") {
      
      ## need to load the [sobj] RDS file if 'CreateSeuratObject' step is skipped:
    } else {
      print(paste0("[SKIP] 'sr_01_mk_sobj' step:"))
      
      if (!isTRUE(tcell_only)) {
        fn_rds_sobj <- paste0(dir_res_sr,"sobj_",prj_gex,"_",smpid,"_",sr_dtype,".rds")
      } else {
        fn_rds_sobj <- paste0(dir_res_sr,"sobj_",prj_gex,"_",smpid,"_",sr_dtype,"_vdj",".rds")
      }
      if (file.exists(fn_rds_sobj)) {
        sobj <- readRDS(fn_rds_sobj)
        print(paste0("[DONE] load the RDS: ",fn_rds_sobj))
        
        print(head(colnames(sobj))) #checkpoint
        print(tail(colnames(sobj))) #checkpoint
        
        ## check the unique samples included in the merged SeuratObject
        unique(sapply(X = strsplit(colnames(sobj), split = "_"), FUN = "[", 1)) #checkpoint
        ## check the number of cells from each sample
        table(sobj$orig.ident)
        
      } else {
        print(paste0("[CHECK] RDS file does not exist: ",fn_rds_sobj))
      }
      print(":")
    } #closer: if (isTRUE(sr_01_mk_sobj)) {
    
    ## add a new metadata column based on GCAR expression ####
    if (sr_mode=="each") {
      if (isTRUE(sr_01a_up_sobj)) {
        ## load the RDS: [sobj] - before "flt" step
        fn_rds_sobj <- paste0(dir_res_sr,"sobj_",prj,"_",smpid,"_",sr_dtype,".rds")
        if (file.exists(fn_rds_sobj)) {
          sobj <- readRDS(file = fn_rds_sobj)
          print(paste0("[DONE] load the RDS: ",fn_rds_sobj))
        }
        
        ## subset the SeuratObject if GCAR expression is above zero
        sobj_gcar <- subset(sobj, subset=GCAR>0)
        
        ## extract cell barcode names
        cb <- rownames(sobj_gcar@meta.data)
        
        ## add a column to the existing metadata of the SeuratObject
        sobj$barcode <- rownames(sobj@meta.data)
        sobj@meta.data <- sobj@meta.data %>%
          mutate(
            expr_GCAR=ifelse((sobj$barcode %in% cb),
                             yes="pos",
                             no="neg"
            )
          )
        print(head(sobj@meta.data)) #checkpoint
        print(tail(sobj@meta.data)) #checkpoint
        
        ## update the RDS: [sobj] - w/ GCAR expression
        fn_rds_sobj <- paste0(dir_res_sr,"sobj_",prj,"_",smpid,"_",sr_dtype,".rds")
        saveRDS(object = sobj, file = fn_rds_sobj)
        print(paste0("[DONE] update the RDS: ",fn_rds_sobj," - w/ GCAR expr"))
      } #closer: if (isTRUE(sr_01a_up_sobj)) {
    } #closer: if (sr_mode=="each") {
    
    ## check the distribution of each modality
    # sobj$orig.ident #checkpoint
    ## [nCount_RNA] ##
    # hist(sobj$nCount_RNA)
    print(summary(sobj$nCount_RNA))
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    # 109     441    1543    2778    4055   41443   #[Harvest] - minCell=3 & minFeat=100
    
    ## [nFeature_RNA] ##
    # hist(sobj$nFeature_RNA)
    print(summary(sobj$nFeature_RNA))
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    # 100     246     958    1192    1794    7780   #[Harvest] - minCell=3 & minFeat=100
    
    ## [percent.mt] ##
    # hist(sobj$percent.mt)
    print(summary(sobj$percent.mt))
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    # 0.000   2.886   5.193  11.047  12.955  90.589   #[Harvest] - minCell=3 & minFeat=100
    
    
    if (isTRUE(sr_mk_plt_qc)) {
      ## perform QC based on the number of detected molecules for each modality and mitochondrial percentage
      if (isTRUE(sr_mk_plt_1)) {
        ## plot 1: Visualize QC metrics as a violin plot ####
        fn_pdf_vln <- paste0(mydir_smp, "1_vlnplot_", smpid, ".pdf")
        if (!file.exists(fn_pdf_vln)) {
          if (sr_mode=="each") {
            pdf(file = fn_pdf_vln, width = 7, height = 8)
            print(
              VlnPlot(sobj, 
                      features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), 
                      ncol = 3,
                      log = TRUE, 
                      pt.size = .1, 
                      alpha = .4,
                      cols = "light blue" #"transparent"
              ) + 
                NoLegend()
            )
          } else if (sr_mode=="mrg") {
            pdf(file = fn_pdf_vln, width = 15, height = 8)
            print(
              VlnPlot(sobj, 
                      features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), 
                      ncol = 3,
                      log = TRUE, 
                      pt.size = .1, 
                      alpha = .4,
                      cols = "light blue", #"transparent"
                      split.by = "orig.ident",
                      group.by = "orig.ident"
              ) + 
                NoLegend()
            )
          }
          dev.off()
          print(paste0("[DONE] save as PDF: ",fn_pdf_vln))
        }
      } #closer: if (isTRUE(sr_mk_plt_1)) {
      
      if (isTRUE(sr_mk_plt_1ab)) {
        ## make a data frame for the detailed QC - using my code ####
        md_sobj <- as.data.frame(sobj@meta.data)
        dim(md_sobj) #18430 4 #nCell nColumn
        head(md_sobj,3)
        ## make a column: [cell_barcode]
        md_sobj$cell_barcode <- rownames(md_sobj)
        dim(md_sobj) #18430 5 #nCell nColumn
        head(md_sobj,3)
        ## melt the data frame
        m_md_sobj <- melt(data = md_sobj, 
                          id.vars = c("orig.ident","cell_barcode"), 
                          measure.vars = c("nCount_RNA","nFeature_RNA","percent.mt")
        )
        dim(m_md_sobj) #55290 4
        print(head(m_md_sobj,3)) #checkpoint
        print(tail(m_md_sobj,3)) #checkpoint
        
        ## summarize statistical values: [median],[q25],[q75]
        ## [nCount_RNA] ##
        med_nCntRNA <- round(median(subset(m_md_sobj, variable=="nCount_RNA")$value),2)
        q25_nCntRNA <- round(quantile(subset(m_md_sobj, variable=="nCount_RNA")$value, probs = 0.25, na.rm = T),2)
        q75_nCntRNA <- round(quantile(subset(m_md_sobj, variable=="nCount_RNA")$value, probs = 0.75, na.rm = T),2)
        ## [nFeature_RNA] ##
        med_nFeatRNA <- round(median(subset(m_md_sobj, variable=="nFeature_RNA")$value),2)
        q25_nFeatRNA <- round(quantile(subset(m_md_sobj, variable=="nFeature_RNA")$value, probs = 0.25, na.rm = T),2)
        q75_nFeatRNA <- round(quantile(subset(m_md_sobj, variable=="nFeature_RNA")$value, probs = 0.75, na.rm = T),2)
        ## [percent.mt] ##
        med_pctMT <- round(median(subset(m_md_sobj, variable=="percent.mt")$value),2)
        q25_pctMT <- round(quantile(subset(m_md_sobj, variable=="percent.mt")$value, probs = 0.25, na.rm = T),2)
        q75_pctMT <- round(quantile(subset(m_md_sobj, variable=="percent.mt")$value, probs = 0.75, na.rm = T),2)
        
        
        ## visualize: QC plots ####
        ## plot 1a: [ggvbox_qc] ####
        mytitle <- paste0("[",prjid,"] QC metrics - '",smpid,"'")
        mysubt <- paste0("sobj created w/ minCell=",minCell," & minFeat=",minFeat,
                         "\n{nCountRNA}: med=",med_nCntRNA," | q25=",q25_nCntRNA," | q75=",q75_nCntRNA,
                         "\n{nFeature_RNA}: med=",med_nFeatRNA," | q25=",q25_nFeatRNA," | q75=",q75_nFeatRNA,
                         "\n{percent.mt}: med=",med_pctMT," | q25=",q25_pctMT," | q75=",q75_pctMT
        )
        mywidth <- 5
        ggvbox_qc <- 
          ggplot(data = m_md_sobj, aes(x=variable,y=value)) +
          theme_bw() +
          theme(
            strip.background = element_rect(fill = "transparent")
          ) +
          geom_violin(alpha=.1) +
          geom_boxplot(width=.2, fill="light blue", alpha=.4) + 
          facet_wrap(~variable, scales = "free", nrow = 1) +
          xlab(label = "") +
          ylab(label = "") +
          theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
          labs(
            title = mytitle, subtitle = mysubt
          )
        if (sr_mode=="mrg") {
          mywidth <- 20
          ggvbox_qc <- ggvbox_qc +
            facet_wrap(~orig.ident, nrow = 1, scales = "fixed")
          ggvbox_qc_freeY <- ggvbox_qc +
            facet_wrap(~orig.ident, nrow = 1, scales = "free_y")
        }
        
        ## save as PDF: [ggvbox_qc]
        fn_pdf_ggvbox_qc <- paste0(mydir_smp, "1a_ggvbox_sobj-QC_",smpid,".pdf")
        if (!file.exists(fn_pdf_ggvbox_qc)) {
          ggsave(filename = fn_pdf_ggvbox_qc, plot = ggvbox_qc, device = "pdf", width = mywidth, height = 6.5, units = "in", dpi = 300)
          print(paste0("[DONE] save as PDF: ",fn_pdf_ggvbox_qc))
        }
        ## save as PDF: [ggvbox_qc_freeY]
        if (sr_mode=="mrg") {
          fn_pdf_ggvbox_qc <- paste0(mydir_smp, "1a_ggvbox_sobj-QC_",smpid,"_freeY",".pdf")
          if (!file.exists(fn_pdf_ggvbox_qc)) {
            ggsave(filename = fn_pdf_ggvbox_qc, plot = ggvbox_qc_freeY, device = "pdf", width = mywidth, height = 6.5, units = "in", dpi = 300)
            print(paste0("[DONE] save as PDF: ",fn_pdf_ggvbox_qc))
          }
        } #closer: if (sr_mode=="mrg") {
        
        
        ## plot 1b: [ggvbox_qc_jit] ####
        mytitle <- paste0("[",prjid,"] QC metrics - '",smpid,"'")
        mysubt <- paste0("sobj created w/ minCell=",minCell," & minFeat=",minFeat,
                         "\n{nCountRNA}: med=",med_nCntRNA," | q25=",q25_nCntRNA," | q75=",q75_nCntRNA,
                         "\n{nFeature_RNA}: med=",med_nFeatRNA," | q25=",q25_nFeatRNA," | q75=",q75_nFeatRNA,
                         "\n{percent.mt}: med=",med_pctMT," | q25=",q25_pctMT," | q75=",q75_pctMT
        )
        mywidth <- 5
        ggvbox_qc_jit <- 
          ggplot(data = m_md_sobj, aes(x=variable,y=value)) +
          theme_bw() +
          theme(
            strip.background = element_rect(fill = "transparent")
          ) +
          geom_violin(alpha=.1) +
          geom_boxplot(width=.2, fill="light blue", alpha=.4) + 
          geom_jitter(width=.2, alpha=.1) +
          facet_wrap(~variable, scales = "free", nrow = 1) +
          xlab(label = "") +
          ylab(label = "") +
          theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
          labs(
            title = mytitle, subtitle = mysubt
          )
        if (sr_mode=="mrg") {
          mywidth <- 20
          ggvbox_qc_jit <- ggvbox_qc_jit +
            facet_wrap(~orig.ident, nrow = 1, scales = "fixed")
          ggvbox_qc_jit_freeY <- ggvbox_qc_jit +
            facet_wrap(~orig.ident, nrow = 1, scales = "free_y")
        }
        
        ## save as PDF: [ggvbox_qc_jit]
        fn_pdf_ggvbox_qc_jit <- paste0(mydir_smp, "1b_ggvbox_sobj-QC_",smpid,"_wJitter.pdf")
        if (!file.exists(fn_pdf_ggvbox_qc_jit)) {
          ggsave(filename = fn_pdf_ggvbox_qc_jit, plot = ggvbox_qc_jit, device = "pdf", width = mywidth, height = 6.5, units = "in", dpi = 300)
          print(paste0("[DONE] save as PDF: ",fn_pdf_ggvbox_qc_jit))
        }
        ## save as PDF: [ggvbox_qc_freeY]
        if (sr_mode=="mrg") {
          fn_pdf_ggvbox_qc_jit <- paste0(mydir_smp, "1b_ggvbox_sobj-QC_",smpid,"_wJitter_freeY.pdf")
          if (!file.exists(fn_pdf_ggvbox_qc_jit)) {
            ggsave(filename = fn_pdf_ggvbox_qc_jit, plot = ggvbox_qc_jit_freeY, device = "pdf", width = mywidth, height = 6.5, units = "in", dpi = 300)
            print(paste0("[DONE] save as PDF: ",fn_pdf_ggvbox_qc_jit))
          }
        } #closer: if (sr_mode=="mrg") {
      } #closer: if (isTRUE(sr_mk_plt_1ab)) {
    } #closer: if (isTRUE(sr_mk_plt)) {
    
    
    ## note: From [step 2], the code will run each step for every sample included in the [list_smpid] for better efficiency
    ## step 2 | filter the SeuratObject ####
    if (isTRUE(sr_02_flt)) {
      print(paste0("--- STEP 2 | filter the SeuratObject | smp_id: [",smpid,"]"," (",which(list_smpid==smpid)," /",length(list_smpid),")"," ---"))
      
      ## load the RDS: [sobj]
      fn_rds_sobj <- paste0(dir_res_sr,"sobj_",prj,"_",smpid,"_",sr_dtype,".rds")
      if (file.exists(fn_rds_sobj)) {
        sobj <- readRDS(fn_rds_sobj)
        print(paste0("[DONE] load the RDS: ",fn_rds_sobj))
      } else {
        print(paste0("[CHECK] RDS file does not exist: ",fn_rds_sobj))
      }
      
      ## subset cells w/ VDJ data ####
      if (isTRUE(tcell_only)) {
        ## load the "source" metadata 
        fn_md_vdj <- paste0(dir_res,"scRepertoire/","metadata_clones_GCAR1_PT01_VDJ_nSmp7.txt")
        
        md_vdj <- read.table(file = fn_md_vdj, sep = "\t", header = T, row.names = NULL, quote = NULL)
        print(dim(md_vdj)) #81520 10
        print(head(md_vdj,3))
        print(tail(md_vdj,3))
        ## add columns: [cell_VDJ],[sample_new_VDJ]
        md_vdj$cell_VDJ <- "VDJ"
        md_vdj$sample_new_VDJ <- paste(md_vdj$sample_new, md_vdj$cell_VDJ, sep = "_")
        print(head(md_vdj,3))
        print(dim(md_vdj)) #81520 12
        
        ## select columns to add to the existing Seurat object
        md_vdj_sel <- subset(md_vdj, select=c(sample_new,cell_barcode,sample_new_VDJ,cell_VDJ))
        ## set row names
        rownames(md_vdj_sel) <- md_vdj$cell_barcode
        print(head(md_vdj_sel,3))
        print(dim(md_vdj_sel)) #81521 2
        # ---
        
        ## define the column name from the "source" Seurat Object
        # cname_clst <- "cell_VDJ"
        cname_clst <- colnames(md_vdj_sel)
        
        ## checkpoint:
        print("[BEFORE] add the metadata: ")
        print(head(sobj@meta.data))
        print(tail(sobj@meta.data))
        
        ## add the extracted metadata column to the newer version of the UMAP Seurat Object
        sobj <- AddMetaData(object = sobj,
                            metadata = md_vdj_sel,
                            col.name = cname_clst
        )
        
        ## checkpoint:
        print("[AFTER] add the metadata: ")
        print(head(sobj@meta.data))
        print(tail(sobj@meta.data))
        
        ## update the RDS file: [sobj]
        print(paste0("[START] update 'sobj' RDS | ", Sys.time()))
        saveRDS(object = sobj, file = fn_rds_sobj)
        print(paste0("[DONE] update the RDS: ",fn_rds_sobj))
        print(paste0("[END] update 'sobj' RDS | ", Sys.time()))
        
        ## save as a new RDS file: [sobj_vdj]
        sr_dtype <- "raw_vdj" #raw | filtered
        sobj_vdj <- subset(sobj, subset=cell_VDJ=="VDJ")
        print(dim(sobj_vdj))
        fn_rds_sobj_vdj <- paste0(dir_res_sr,"sobj_",prj,"_",smpid,"_",sr_dtype,".rds")
        print(paste0("[START] save as a new 'sobj' RDS | ", Sys.time()))
        saveRDS(object = sobj_vdj, file = fn_rds_sobj_vdj)
        print(paste0("[DONE] save as a new the RDS: ",fn_rds_sobj_vdj))
        print(paste0("[END] save as a new 'sobj' RDS | ", Sys.time()))
        
        print(paste0("[DONE] subset the 'sobj': ",fn_rds_sobj_vdj))
        
        # ---
        ## calculate the "median reads (i.e., counts) per sample"
        print(with(sobj_vdj@meta.data, tapply(nCount_RNA,orig.ident, median)))
        
        ## calculate the "median genes (i.e., features) per sample"
        print(with(sobj_vdj@meta.data, tapply(nFeature_RNA,orig.ident, median)))        
        # ---
        
      } #closer: if (isTRUE(tcell_only)) {
      
      
      ## set thresholds for each measure (based on the violin plots): 
      ## [nCount_RNA]: number of cells
      ## [nFeature_RNA]: number of genes
      ## [percent.mt]: select cells with the mitochondrial percentage lower than 25%
      if (smpid=="Harvest") {
        # thrCell <- 0
        # thrFeat <- 0
        thrMT <- 5
      } else if (smpid=="Enriched") {
        # thrCell <- 0
        # thrFeat <- 0
        thrMT <- 5
      } else if (smpid=="D08") {
        # thrCell <- 0
        # thrFeat <- 0
        thrMT <- 5
      } else if (smpid=="D15") {
        # thrCell <- 0
        # thrFeat <- 0
        thrMT <- 5
      } else if (smpid=="D22") {
        # thrCell <- 0
        # thrFeat <- 0
        thrMT <- 5
      } else if (smpid=="D28") {
        # thrCell <- 0
        # thrFeat <- 0
        thrMT <- 5
      } else if (smpid=="D42") {
        # thrCell <- 0
        # thrFeat <- 0
        thrMT <- 5
      } else if (smpid=="mrg" | smpid=="mrg7") {
        thrMT <- 5
      }
      
      ## set the Seurat Object for filtering step
      if (!isTRUE(tcell_only)) {
        sobj <- sobj
        sr_dtype <- sr_dtype #raw | filtered
      } else if (isTRUE(tcell_only)) {
        sobj <- sobj_vdj
        sr_dtype <- paste0(sr_dtype,"_vdj") #raw | filtered
      } #closer: if (!isTRUE(tcell_only)) {
      
      ## select cells based on the threshold set above
      sobj_flt <- subset(
        x = sobj,
        subset = 
          # nCount_RNA >thrCell & 
          # nFeatureRNA >thrFeat & 
          percent.mt <thrMT
      )
      
      ## save as RDS: [sobj_flt] - initial ####
      fn_rds_sobj_flt <- paste0(dir_res_sr,"sobj_flt_",prj,"_",smpid,"_",sr_dtype,".rds")
      saveRDS(object = sobj_flt, file = fn_rds_sobj_flt)
      print(paste0("[DONE] save as RDS: ",fn_rds_sobj_flt))
      
      ## make a comparison table - before vs. after the subsetting step
      ## from [STEP 1]:
      # minCell <- 1 #3(default)
      # minFeat <- 100 #200(default)
      # dim(sobj_flt) #nUMI nCell
      
      df_comp_sobj <- data.frame(
        sample=smpid,
        parameter_minCell=minCell,
        parameter_minFeat=minFeat,
        nUMI_sobj=nrow(sobj),
        nCell_sobj=ncol(sobj),
        threshold_pctMT=thrMT,
        nUMI_sobj_flt=nrow(sobj_flt),
        nCell_sobj_flt=ncol(sobj_flt)
      )
      ## add columns: [nUMI_excluded],[nCell_excluded]
      df_comp_sobj <- df_comp_sobj %>%
        mutate(
          nUMI_excluded=nUMI_sobj-nUMI_sobj_flt,
          nCell_excluded=nCell_sobj-nCell_sobj_flt
        ) %>% 
        as.data.frame()
      
      fn_txt_comp_sobj <- paste0(dir_res_sr,"compare_sobj_",prj,"_nSmp",length(list_smpid),".txt")
      if (!file.exists(fn_txt_comp_sobj)) {
        if (smpid==list_smpid[1]) {
          write.table(x = df_comp_sobj, file = fn_txt_comp_sobj, sep = "\t", col.names = T, row.names = F, quote = F, append = F)
        } else {
          write.table(x = df_comp_sobj, file = fn_txt_comp_sobj, sep = "\t", col.names = F, row.names = F, quote = F, append = T)
        }
        print(paste0("[DONE] save as TXT: ",fn_txt_comp_sobj," - ",smpid," (",which(list_smpid==smpid)," /",length(list_smpid),")"))
      } #closer: if (!file.exists(fn_txt_comp_sobj)) {
      
      print(":")
    } #closer: if (isTRUE(sr_02_flt)) {
    
    
    ## step 3 | normalize the SeuratObject - standard steps ####
    if (isTRUE(sr_03_norm)) {
      print(paste0("--- STEP 3 | normalize the SeuratObject - standard steps | smp_id: [",smpid,"]"," (",which(list_smpid==smpid)," /",length(list_smpid),")"," ---"))
      
      ## set the Seurat Object for normalization step
      if (!isTRUE(tcell_only)) {
        sr_dtype <- sr_dtype #raw | filtered
      } else if (isTRUE(tcell_only)) {
        sr_dtype <- paste0(sr_dtype,"_vdj") #raw | filtered
      }
      
      ## load the RDS: [sobj_flt]
      fn_rds_sobj_flt <- paste0(dir_res_sr,"sobj_flt_",prj,"_",smpid,"_",sr_dtype,".rds")
      if (file.exists(fn_rds_sobj_flt)) {
        sobj_flt <- readRDS(fn_rds_sobj_flt)
        print(paste0("[DONE] load the RDS: ",fn_rds_sobj_flt))
      } else {
        print(paste0("[CHECK] RDS file does not exist: ",fn_rds_sobj_flt))
      }
      
      ## re-set the assay slot
      DefaultAssay(sobj_flt) <- "RNA"
      # sobj_flt[["SCT"]] <- NULL
      
      ## set the Seurat format version
      assay_seurat <- "RNA"
      
      # After removing unwanted cells from the dataset, the next step is to normalize the data. 
      # By default, we employ a global-scaling normalization method “LogNormalize” that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result. 
      # In Seurat v5, Normalized values are stored in pbmc[["RNA"]]$data.
      
      print(paste0("[START] NormalizeData | ", Sys.time()))
      sobj_flt <- NormalizeData(sobj_flt, 
                                assay = assay_seurat, 
                                normalization.method = "LogNormalize", 
                                scale.factor = 10000
      )
      print(paste0("[END] NormalizeData | ", Sys.time()))
      
      ## update the [sobj_flt] RDS
      fn_rds_sobj_flt <- paste0(dir_res_sr,"sobj_flt_",prj,"_",smpid,"_",sr_dtype,".rds")
      sobj_flt <- saveRDS(object = sobj_flt, file = fn_rds_sobj_flt)
      print(paste0("[DONE] update the RDS: ",fn_rds_sobj_flt," - after 'NormalizeData' step"))
      
    } #closer: if (isTRUE(sr_03_norm)) {
    
    
    ## step 4 | Identification of highly variable features (feature selection) ####
    # We next calculate a subset of features that exhibit high cell-to-cell variation in the dataset (i.e, they are highly expressed in some cells, and lowly expressed in others). 
    # We and others have found that focusing on these genes in downstream analysis helps to highlight biological signal in single-cell datasets.
    if (isTRUE(sr_04_fsel)) {
      print(paste0("--- STEP 4 | feature selection | smp_id: [",smpid,"]"," (",which(list_smpid==smpid)," /",length(list_smpid),")"," ---"))
      
      ## set the Seurat Object for normalization step
      if (!isTRUE(tcell_only)) {
        sr_dtype <- sr_dtype #raw | filtered
      } else if (isTRUE(tcell_only)) {
        sr_dtype <- paste0(sr_dtype,"_vdj") #raw | filtered
      }
      
      ## load the RDS: [sobj_flt]
      fn_rds_sobj_flt <- paste0(dir_res_sr,"sobj_flt_",prj,"_",smpid,"_",sr_dtype,".rds")
      if (file.exists(fn_rds_sobj_flt)) {
        sobj_flt <- readRDS(fn_rds_sobj_flt)
        print(paste0("[DONE] load the RDS: ",fn_rds_sobj_flt))
      } else {
        print(paste0("[CHECK] RDS file does not exist: ",fn_rds_sobj_flt))
      }
      
      ## re-set the assay slot
      DefaultAssay(sobj_flt) <- "RNA"
      # sobj_flt[["SCT"]] <- NULL
      
      ## check the Seurat format version
      print(sobj_flt@version)
      # print(class(sobj_flt[["RNA"]]))
      ## convert the v5 assay to v3 assay
      ## - source: "https://satijalab.org/seurat/articles/seurat5_essential_commands.html"
      # sobj_flt[["RNA3"]] <- as(object = sobj_flt[["RNA"]], Class = "Assay")
      # print(class(sobj_flt[["RNA3"]]))
      
      ## set the Seurat format version
      assay_seurat <- "RNA"
      
      ## run the [FindVariableFeatures] function
      ## - "Basically you should not be using the integrated assay to find variable features"
      ##   (source: "https://github.com/satijalab/seurat/issues/2042")
      sobj_flt <- FindVariableFeatures(sobj_flt, 
                                       assay = assay_seurat, 
                                       selection.methcod = "vst", 
                                       nfeatures = 3000, 
                                       verbose=T)
      
      ## update the [sobj_flt] RDS
      fn_rds_sobj_flt <- paste0(dir_res_sr,"sobj_flt_",prj,"_",smpid,"_",sr_dtype,".rds")
      sobj_flt <- saveRDS(object = sobj_flt, file = fn_rds_sobj_flt)
      print(paste0("[DONE] update the RDS: ",fn_rds_sobj_flt," - after 'FindVariableFeatures' step"))
      
      
      if (isTRUE(sr_mk_plt_2)) {
        ## Identify the 10 most highly variable genes
        top10 <- head(VariableFeatures(sobj_flt, assay = assay_seurat), 10)
        print(top10)
        
        ## plot variable features with and without labels
        plot1 <- VariableFeaturePlot(sobj_flt, assay = assay_seurat, selection.method = "vst")
        plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
        # plot1 + plot2
        
        ## plot 2: [ggvbox_qc] ####
        fn_pdf_var <- paste0(mydir_smp,"2_scat_variableFeat_",smpid,".pdf")
        if (!file.exists(fn_pdf_var)) {
          pdf(file = fn_pdf_var, width = 16, height = 7)
          print(
            CombinePlots(plots = list(plot1, plot2), ncol = 2)
          )
          dev.off()
          print(paste0("[DONE] save as PDF: ",fn_pdf_var))
        }
      } #closer: if (isTRUE(sr_mk_plt_2)) {
      
      print(":")
    } #closer:  if (isTRUE(sr_04_fsel)) {
    
    
    ## step 5 | Scaling the data ####
    # Next, we apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA. 
    
    # The ScaleData() function:
    # - Shifts the expression of each gene, so that the mean expression across cells is 0
    # - Scales the expression of each gene, so that the variance across cells is 1
    #   -- This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
    # - The results of this are stored in pbmc[["RNA"]]$scale.data
    # - By default, only variable features are scaled.
    # - You can specify the features argument to scale additional features
    if (isTRUE(sr_05_scale)) {
      print(paste0("--- STEP 5 | scale the data | smp_id: [",smpid,"]"," (",which(list_smpid==smpid)," /",length(list_smpid),")"," ---"))
      
      ## set the Seurat Object for normalization step
      if (!isTRUE(tcell_only)) {
        sr_dtype <- sr_dtype #raw | filtered
      } else if (isTRUE(tcell_only)) {
        sr_dtype <- paste0(sr_dtype,"_vdj") #raw | filtered
      }
      
      ## load the RDS: [sobj_flt]
      fn_rds_sobj_flt <- paste0(dir_res_sr,"sobj_flt_",prj,"_",smpid,"_",sr_dtype,".rds")
      if (file.exists(fn_rds_sobj_flt)) {
        sobj_flt <- readRDS(fn_rds_sobj_flt)
        print(paste0("[DONE] load the RDS: ",fn_rds_sobj_flt))
      } else {
        print(paste0("[CHECK] RDS file does not exist: ",fn_rds_sobj_flt))
      }
      
      ## re-set the assay slot
      DefaultAssay(sobj_flt) <- "RNA"
      # sobj_flt[["SCT"]] <- NULL
      
      ## check the Seurat format version
      print(sobj_flt@version)
      # print(class(sobj_flt[["RNA"]]))
      ## convert the v5 assay to v3 assay
      ## - source: "https://satijalab.org/seurat/articles/seurat5_essential_commands.html"
      # sobj_flt[["RNA3"]] <- as(object = sobj_flt[["RNA"]], Class = "Assay")
      # print(class(sobj_flt[["RNA3"]]))
      
      ## set the Seurat format version
      assay_seurat <- "RNA"
      
      all.genes <- rownames(sobj_flt)
      
      print(paste0("[START] ScaleData | ", Sys.time()))
      sobj_flt <- ScaleData(sobj_flt, assay = assay_seurat, features = all.genes)
      print(paste0("[END] ScaleData | ", Sys.time()))
      
      ## update the [sobj_flt] RDS
      fn_rds_sobj_flt <- paste0(dir_res_sr,"sobj_flt_",prj,"_",smpid,"_",sr_dtype,".rds")
      sobj_flt <- saveRDS(object = sobj_flt, file = fn_rds_sobj_flt)
      print(paste0("[DONE] update the RDS: ",fn_rds_sobj_flt," - after 'ScaleData' step"))
      
      print(":")
    } #closer: if (isTRUE(sr_05_scale)) {
    
    
    ## step 6 | normalize the SeuratObject - using [sctransform] ####
    # - Biological heterogeneity in single-cell RNA-seq data is often confounded by technical factors including sequencing depth. 
    # - The number of molecules detected in each cell can vary significantly between cells, even within the same celltype. 
    # - Interpretation of scRNA-seq data requires effective pre-processing and normalization to remove this technical variability.
    
    # - In our manuscript we introduce a modeling framework for the normalization and variance stabilization of molecular count data from scRNA-seq experiments. 
    # - This procedure omits the need for heuristic steps including pseudocount addition or log-transformation and improves common downstream analytical tasks 
    #   such as variable gene selection, dimensional reduction, and differential expression. We named this method [sctransform].
    
    # The use of SCTransform replaces the need to run [NormalizeData], [FindVariableFeatures], or [ScaleData] (described below.)
    
    if (isTRUE(sr_06_sct)) {
      print(paste0("--- STEP 6 | normalize the SeuratObject - using [sctransform] | smp_id: [",smpid,"]"," (",which(list_smpid==smpid)," /",length(list_smpid),")"," ---"))
      
      ## set the Seurat Object for normalization step
      if (!isTRUE(tcell_only)) {
        sr_dtype <- sr_dtype #raw | filtered
      } else if (isTRUE(tcell_only)) {
        sr_dtype <- paste0(sr_dtype,"_vdj") #raw | filtered
      }
      
      ## load the RDS: [sobj_flt]
      fn_rds_sobj_flt <- paste0(dir_res_sr,"sobj_flt_",prj,"_",smpid,"_",sr_dtype,".rds")
      if (file.exists(fn_rds_sobj_flt)) {
        sobj_flt <- readRDS(fn_rds_sobj_flt)
        print(paste0("[DONE] load the RDS: ",fn_rds_sobj_flt))
      } else {
        print(paste0("[CHECK] RDS file does not exist: ",fn_rds_sobj_flt))
      }
      
      ## step 6-1: Apply sctransform normalization #skip
      ## store mitochondrial percentage in object meta data
      # patt_mt <- "^MT-"
      # sobj_flt <- PercentageFeatureSet(sobj_flt, pattern = patt_mt, col.name = "percent.mt")
      
      ## run sctransform
      print(paste0("[START] SCTransform | ", Sys.time()))
      sobj_flt <- SCTransform(sobj_flt, vars.to.regress = "percent.mt", verbose = T)
      print(paste0("[END] SCTransform | ", Sys.time()))
      
      ## update the [sobj_flt] RDS
      fn_rds_sobj_flt <- paste0(dir_res_sr,"sobj_flt_",prj,"_",smpid,"_",sr_dtype,"_res",res_clst,".rds")
      sobj_flt <- saveRDS(object = sobj_flt, file = fn_rds_sobj_flt)
      print(paste0("[DONE] update the RDS: ",fn_rds_sobj_flt," - after 'SCTransform' step"))
      
      
      ## step 6-2: Perform dimensionality reduction by PCA and UMAP embedding #skip
      # These are now standard steps in the Seurat workflow for visualization and clustering
      # sobj_flt <- RunPCA(sobj_flt, verbose = T)
      # sobj_flt <- RunUMAP(sobj_flt, dims = 1:30, verbose = T)
      # 
      # sobj_flt <- FindNeighbors(sobj_flt, dims = 1:30, verbose = T)
      # sobj_flt <- FindClusters(sobj_flt, verbose = T)
      # DimPlot(sobj_flt, label = TRUE)
      
      ## check the normalized values stored for [sctransform]
      # sobj_flt[["SCT"]]$scale.data
      
      ## Users can individually annotate clusters based on canonical markers. 
      ## However, the sctransform normalization reveals sharper biological distinctions compared to the standard Seurat workflow, in a few ways:
      
      # - Clear separation of at least 3 "CD8 T cell" populations (naive, memory, effector), based on CD8A, GZMK, CCL5, CCR7 expression
      # - Clear separation of three "CD4 T cell" populations (naive, memory, IFN-activated) based on S100A4, CCR7, IL32, and ISG15
      # - Additional developmental sub-structure in "B cell" cluster, based on TCL1A, FCER2
      # - Additional separation of "NK cells" into CD56dim vs. CD56bright clusters, based on XCL1 and FCGR3A
      
      ## These are now standard steps in the Seurat workflow for visualization and clustering
      ## Visualize canonical marker genes as violin plots.
      # VlnPlot(
      #   sobj_flt, 
      #   features = c(
      #     "CD8A", "GZMK", "CCL5", "CCR7", #CD8 T cells
      #     "S100A4", "IL32", "ISG15", #CD4 T cells
      #     "TCL1A", "FCER2", #B cells
      #     "NCAM1", "XCL1", "FCGR3A" #NK cells
      #     # "ANXA1", "CD3D"
      #   ),
      #   pt.size = 0.2, ncol = 4
      # )
      print(":")
    } #closer: if (isTRUE(sr_06_sct)) {
    
    ## step 7 | Perform linear dimensional reduction ####
    # Next we perform PCA on the scaled data. By default, only the previously determined variable features are used as input, 
    #  but can be defined using features argument if you wish to choose a different subset 
    #  (if you do want to use a custom subset of features, make sure you pass these to "ScaleData" first).
    if (isTRUE(sr_07_pca)) {
      
      ## step 6a: run PCA ####
      # For the first principal components, Seurat outputs a list of genes with the most positive and negative loadings, 
      #  representing modules of genes that exhibit either correlation (or anti-correlation) across single-cells in the dataset.
      print(paste0("--- STEP 7 | linear dimensional reduction (PCA) | smp_id: [",smpid,"]"," (",which(list_smpid==smpid)," /",length(list_smpid),")"," ---"))
      
      ## set the Seurat Object for normalization step
      if (!isTRUE(tcell_only)) {
        sr_dtype <- sr_dtype #raw | filtered
      } else if (isTRUE(tcell_only)) {
        sr_dtype <- paste0(sr_dtype,"_vdj") #raw | filtered
      }
      
      ## load the RDS: [sobj_flt]
      fn_rds_sobj_flt <- paste0(dir_res_sr,"sobj_flt_",prj,"_",smpid,"_",sr_dtype,".rds")
      if (file.exists(fn_rds_sobj_flt)) {
        sobj_flt <- readRDS(fn_rds_sobj_flt)
        print(paste0("[DONE] load the RDS: ",fn_rds_sobj_flt))
      } else {
        print(paste0("[CHECK] RDS file does not exist: ",fn_rds_sobj_flt))
      }
      
      ## run the [RunPCA] function
      print(paste0("[START] RunPCA | ", Sys.time()))
      sobj_flt <- RunPCA(sobj_flt, features = VariableFeatures(object = sobj_flt))
      print(paste0("[END] RunPCA | ", Sys.time()))
      
      ## update the [sobj_flt] RDS
      fn_rds_sobj_flt <- paste0(dir_res_sr,"sobj_flt_",prj,"_",smpid,"_",sr_dtype,"_res",res_clst,".rds")
      sobj_flt <- saveRDS(object = sobj_flt, file = fn_rds_sobj_flt)
      print(paste0("[DONE] update the RDS: ",fn_rds_sobj_flt," - after 'RunPCA' step"))
      
      
      ## Seurat provides several useful ways of visualizing both cells and features that define the PCA, including VizDimReduction(), DimPlot(), and DimHeatmap()
      ## Examine and visualize PCA results a few different ways
      print(sobj_flt[["pca"]], dims = 1:5, nfeatures = 5)
      
      if (isTRUE(sr_mk_plt_3)) {
        ## plot 3a: PC dot plots (top-scored genes in each PC) ####
        pcnum_from <- 1
        pcnum_to <- 5
        nPC <- (pcnum_to-pcnum_from)+1
        nTopGene <- 50
        fn_pdf_pc_dot <- paste0(mydir_smp, "3a_dot_top",nTopGene,"Genes_PC",pcnum_from,"-to-PC",pcnum_to,"_",smpid,".pdf")
        if (!file.exists(fn_pdf_pc_dot)) {
          pdf(file = fn_pdf_pc_dot, width = 16.5, height = 9)
          print(
            VizDimLoadings(sobj_flt, 
                           dims = pcnum_from:pcnum_to, #1:5
                           ncol = nPC,
                           nfeatures = nTopGene, 
                           balanced = T, 
                           combine = T, #T(plots in one page) #F(plots are in the individual pages)
                           reduction = "pca"
            )
          )
          dev.off()
          print(paste0("[DONE] save as PDF: ",fn_pdf_pc_dot))
        }
        
        ## plot 3b: PC scatter plot ####
        pcnum_from <- 1
        pcnum_to <- 2
        fn_pdf_pc_scat <- paste0(mydir_smp, "3b_scat_PC",pcnum_from,"-vs-PC",pcnum_to,"_",smpid,".pdf")
        if (!file.exists(fn_pdf_pc_scat)) {
          pdf(file = fn_pdf_pc_scat, width = 9, height = 9)
          print(
            DimPlot(sobj_flt, 
                    reduction = "pca", 
                    dims = c(pcnum_from,pcnum_to),
            ) + 
              NoLegend()
          )
          dev.off()
          print(paste0("[DONE] save as PDF: ",fn_pdf_pc_scat))
        }
        
        ## plot 3c: PC heatmap ####
        pcnum_from <- 1
        pcnum_to <- 5
        fn_pdf_pc_hmp <- paste0(mydir_smp, "3c_hmp_PC",pcnum_from,"-to-PC",pcnum_to,"_",smpid,".pdf")
        if (!file.exists(fn_pdf_pc_hmp)) {
          pdf(file = fn_pdf_pc_hmp, width = 16.5, height = 9)
          print(
            DimHeatmap(sobj_flt, 
                       dims = pcnum_from:pcnum_to, #1:5
                       cells = 500, 
                       balanced = TRUE, 
                       combine = T,
                       ncol = 3
            )
          )
          dev.off()
          print(paste0("[DONE] save as PDF: ",fn_pdf_pc_hmp))
        }
        
        ## step 6b: Determine the ‘dimensionality’ of the dataset ####
        # To overcome the extensive technical noise in any single feature for scRNA-seq data, 
        #  Seurat clusters cells based on their PCA scores, with each PC essentially representing a ‘metafeature’ that combines information across a correlated feature set. 
        # The top principal components therefore represent a robust compression of the dataset. 
        # However, how many components should we choose to include? 10? 20? 100?
        
        ## plot 3d: elbow plot ####
        fn_pdf_elb <- paste0(mydir_smp, "3d_elbowPlot_",smpid,".pdf")
        if (!file.exists(fn_pdf_elb)) {
          pdf(file = fn_pdf_elb, width = 8, height = 5)
          print(
            ElbowPlot(sobj_flt)
          )
          dev.off()
          print(paste0("[DONE] save as PDF: ",fn_pdf_elb))
        }
      } #closer: if (isTRUE(sr_mk_plt_3)) {
      
      print(":")
    } #closer: if (isTRUE(sr_07_pca)) {
    
    
    ## step 8 | Cluster the cells ####
    if (isTRUE(sr_08_clst)) {
      
      print(paste0("--- STEP 8 | cell clsustering | smp_id: [",smpid,"]"," (",which(list_smpid==smpid)," /",length(list_smpid),")"," ---"))
      
      ## set the Seurat Object for normalization step
      if (!isTRUE(tcell_only)) {
        sr_dtype <- sr_dtype #raw | filtered
      } else if (isTRUE(tcell_only)) {
        sr_dtype <- paste0(sr_dtype,"_vdj") #raw | filtered
      }
      
      ## load the RDS: [sobj_flt]
      fn_rds_sobj_flt <- paste0(dir_res_sr,"sobj_flt_",prj,"_",smpid,"_",sr_dtype,".rds")
      if (file.exists(fn_rds_sobj_flt)) {
        sobj_flt <- readRDS(fn_rds_sobj_flt)
        print(paste0("[DONE] load the RDS: ",fn_rds_sobj_flt))
      } else {
        print(paste0("[CHECK] RDS file does not exist: ",fn_rds_sobj_flt))
      }
      
      sobj_flt <- FindNeighbors(sobj_flt, dims = 1:10)
      sobj_flt <- FindClusters(sobj_flt, resolution = res_clst)
      print(head(sobj_flt@meta.data))
      print(tail(sobj_flt@meta.data))
      
      ## update the [sobj_flt] RDS
      fn_rds_sobj_flt <- paste0(dir_res_sr,"sobj_flt_",prj,"_",smpid,"_",sr_dtype,".rds")
      sobj_flt <- saveRDS(object = sobj_flt, file = fn_rds_sobj_flt)
      print(paste0("[DONE] update the RDS: ",fn_rds_sobj_flt," - after 'FindClusters' step"))
      
      
      ## check [cluster IDs] of the first 5 cells
      # print(head(Idents(sobj_flt), 5))
      
      print(":")
    } #closer: if (isTRUE(sr_08_clst)) {
    
    
    ## step 9 | Run non-linear dimensional reduction (UMAP/tSNE) ####
    # Seurat offers several non-linear dimensional reduction techniques, such as tSNE and UMAP, to visualize and explore these datasets. 
    
    # The goal of these algorithms is to learn underlying structure in the dataset, in order to place similar cells together in low-dimensional space. 
    # Therefore, cells that are grouped together within graph-based clusters determined above should co-localize on these dimension reduction plots.
    if (isTRUE(sr_09_umap)) {
      print(paste0("--- STEP 9 | non-linear dimensional reduction (UMAP/tSNE) | smp_id: [",smpid,"]"," (",which(list_smpid==smpid)," /",length(list_smpid),")"," ---"))
      
      ## set the Seurat Object for normalization step
      if (!isTRUE(tcell_only)) {
        sr_dtype <- sr_dtype #raw | filtered
      } else if (isTRUE(tcell_only)) {
        sr_dtype <- paste0(sr_dtype,"_vdj") #raw | filtered
      }
      
      ## load the RDS: [sobj_flt]
      fn_rds_sobj_flt <- paste0(dir_res_sr,"sobj_flt_",prj,"_",smpid,"_",sr_dtype,".rds")
      if (file.exists(fn_rds_sobj_flt)) {
        sobj_flt <- readRDS(fn_rds_sobj_flt)
        print(paste0("[DONE] load the RDS: ",fn_rds_sobj_flt))
      } else {
        print(paste0("[CHECK] RDS file does not exist: ",fn_rds_sobj_flt))
      }
      
      ## set seed
      seed <- 42 #42(default)
      ## run [RunUMAP]
      sobj_flt <- RunUMAP(sobj_flt, dims = 1:10, seed.use = seed)
      
      ## note that you can set `label = TRUE` or use the LabelClusters function to help label
      ## individual clusters
      # DimPlot(sobj_flt, reduction = "umap", label = T)
      
      ## re-set [prj] ####
      if (!isTRUE(tcell_only)) {
        prj <- "GCAR1_PT01_GEX"
      } else {
        prj <- "GCAR1_PT01_GEX-subsetVDJ"
      }
      
      ## save as RDS: [sobj_flt] - after UMAP ####
      fn_rds_sobj_flt_umap <- paste0(dir_res_sr,"sobj_flt_umap_",prj,"_",smpid,"_",sr_dtype,"_res",res_clst,".rds")
      # if (!file.exists(fn_rds_sobj_flt_umap)) {
      saveRDS(object = sobj_flt, file = fn_rds_sobj_flt_umap)
      print(paste0("[DONE] save as RDS:: ",fn_rds_sobj_flt_umap))
      
    } else {
      print(paste0("[SKIP] 'sr_09_umap' step:"))
      
    } #closer: if (isTRUE(sr_09_umap)) {
    
    
    ## visualize using the [sobj_flt_umap]
    if (isTRUE(sr_mk_plt_4)) {
      
      ## colour palette ##
      if (sr_mode=="mrg") {
        # based on [Dark2] (n=8)
        col_gcar_e <- "#666666" #dark_grey
        col_gcar_h <- "#bf5b17" #brown
        col_gcar_d08 <- "#a6d854" #light_green
        col_gcar_d15 <- "#e7298a" #pink
        col_gcar_d22 <- "#e6ab02" #yellow
        col_gcar_d28 <- "#1b9e77" #dark_green
        col_gcar_d42 <- "#386cb0" #blue
        col_smp <- setNames(c(col_gcar_e,col_gcar_h,col_gcar_d08,col_gcar_d15,col_gcar_d22,col_gcar_d28,col_gcar_d42), c("Enriched","Harvest","D08","D15","D22","D28","D42"))
      }
      
      ## plot 4a: rnaUMAP plot ####
      fn_pdf_umap <- paste0(mydir_smp, "4a_rnaUMAP_", smpid, ".pdf")
      if (!file.exists(fn_pdf_umap)) {
        pdf(file = fn_pdf_umap, width = 15, height = 13)
        print(
          DimPlot(object = sobj_flt_umap, 
                  reduction = "umap", 
                  label = T, repel = T,
                  label.size = 6,
                  seed = 42 #make the seed same as the value used in the [RunUMAP]
          )
        )
        dev.off()
        print(paste0("[DONE] save as PDF: ",fn_pdf_umap))
      }
      
      ## plot 4a: rnaUMAP plot (group by each sample) ####
      if (sr_mode=="mrg") {
        fn_pdf_umap <- paste0(mydir_smp, "4a_rnaUMAP_", smpid,"_bySmp",".pdf")
        if (!file.exists(fn_pdf_umap)) {
          pdf(file = fn_pdf_umap, width = 15, height = 13)
          print(
            DimPlot(object = sobj_flt_umap, 
                    reduction = "umap", 
                    label = T, repel = T,
                    label.size = 6,
                    seed = 42, #make the seed same as the value used in the [RunUMAP]
                    group.by = "orig.ident",
                    cols = col_smp
            )
          )
          dev.off()
          print(paste0("[DONE] save as PDF: ",fn_pdf_umap))
        }
      } #closer: if (sr_mode=="mrg") {
      
      ## colour palette ##
      nClst <- length(unique(sobj_flt_umap$seurat_clusters))
      if (nClst<=25) {
        clstColor <- setNames(cols25(n = 25)[1:nClst], sort(unique(as.factor(sobj_flt_umap$seurat_clusters))))
      } else if (nClst>25) {
        clstColor <- setNames(c(cols25(n = 25),kelly(n=22)[3:22])[1:nClst], sort(unique(as.factor(sobj_flt_umap$seurat_clusters))))
      }
      # print(paste0("--- [clstColor]:\n",clstColor)) #checkpoint
      
      ## plot 4b: rnaUMAP plot - w/ cluster numbers ####
      fn_pdf_umap <- paste0(mydir_smp, "4b_rnaUMAP_", smpid, "_wClustNum", ".pdf")
      if (!file.exists(fn_pdf_umap)) {
        pdf(file = fn_pdf_umap, width = 15, height = 13)
        print(
          DimPlot(object = sobj_flt_umap, reduction = "umap", 
                  label = T, repel = T, 
                  cols = clstColor, 
                  seed = 42, #make the seed same as the value used in the [RunUMAP]
                  label.size = 8
          )
        )
        dev.off()
        print(paste0("[DONE] save as PDF: ",fn_pdf_umap))
      }
      
      ## plot 4b: rnaUMAP plot - w/ cluster numbers (group by each sample) ####
      if (sr_mode=="mrg") {
        fn_pdf_umap <- paste0(mydir_smp, "4b_rnaUMAP_", smpid, "_wClustNum","_bySmp", ".pdf")
        if (!file.exists(fn_pdf_umap)) {
          pdf(file = fn_pdf_umap, width = 28, height = 16)
          print(
            DimPlot(object = sobj_flt_umap, reduction = "umap", 
                    label = T, repel = T, 
                    cols = clstColor, 
                    seed = 42, #make the seed same as the value used in the [RunUMAP]
                    label.size = 5,
                    split.by = "orig.ident", 
                    ncol = 4
            )
          )
          dev.off()
          print(paste0("[DONE] save as PDF: ",fn_pdf_umap))
        }
      }
      
      
      ## plot 4c: rnaUMAP plot - w/ label box ####
      fn_pdf_umap <- paste0(mydir_smp, "4c_rnaUMAP_", smpid, "_wClustNum-BoxLabel", ".pdf")
      if (!file.exists(fn_pdf_umap)) {
        pdf(file = fn_pdf_umap, width = 15, height = 13)
        print(
          DimPlot(object = sobj_flt_umap, reduction = "umap", 
                  label = T, repel = T, 
                  cols = clstColor, 
                  seed = 42, #make the seed same as the value used in the [RunUMAP]
                  label.box = T
          )
        )
        dev.off()
        print(paste0("[DONE] save as PDF: ",fn_pdf_umap))
      }        
    } #closer: if (isTRUE(sr_mk_plt_4)) {
    
    ## save workspace ####
    if (isTRUE(sr_save_wspace)) {
      print(paste0("--- save as workspace - after UMAP/tSNE | smp_id: [",smpid,"]"," (",which(list_smpid==smpid)," /",length(list_smpid),")"," ---"))
      
      if (isTRUE(sr_09_umap)) {
        fn_wspace <- paste0(dir_res_sr,"Seurat_afterUMAP_", prj,"_",smpid,"_res",res_clst,"_",today,".RData")        
        print(paste0("[START] ", Sys.time()))
        gc()
        save.image(file = fn_wspace) # elapsed time: xx min. approx.
        print(paste0("[END] ", Sys.time()))
        print(paste0("[DONE] save as RData: ",fn_wspace))
        gc()
        print(fn_wspace)
      } #closer: if (isTRUE(sr_09_umap)) {
      
    } #closer: if (isTRUE(sr_save_wspace)) {
    # } #closer: if (isTRUE(sr_09_umap)) {
    print(":")
  } #closer: for (smpid in list_smpid) {
  
  ## step 9a | after UMAP ####
  if (isTRUE(sr_09a_afterUMAP)) {
    for (smpid in list_smpid) {
      print(paste0("--- STEP 9a | after UMAP | smp_id: [",smpid,"]"," (",which(list_smpid==smpid)," /",length(list_smpid),")"," ---"))
      
      ## load the RDS: [sobj_flt_umap]
      fn_rds_sobj_flt_umap <- paste0(dir_res_sr,"sobj_flt_umap_",prj,"_",smpid,"_",sr_dtype,"_res",res_clst,".rds")
      if (file.exists(fn_rds_sobj_flt_umap)) {
        sobj_flt_umap <- readRDS(fn_rds_sobj_flt_umap)
        print(paste0("[DONE] load the RDS: ",fn_rds_sobj_flt_umap))
        
        ## modify the order of the [orig.ident] variable
        print(unique(sobj_flt_umap$orig.ident))
        sobj_flt_umap$orig.ident <- factor(x = sobj_flt_umap$orig.ident, 
                                           levels = c(
                                             "Enriched","Harvest",
                                             "D08","D15","D22","D28","D42"
                                           )
        )
        print(unique(sobj_flt_umap$orig.ident))
        
      } else {
        print(paste0("[CHECK] RDS file does not exist: ",fn_rds_sobj_flt_umap))
      }
      
      ## plot 5 | visualize using the [sobj_flt_umap]
      if (isTRUE(sr_mk_plt_5)) {
        ## set seed
        seed <- 42 #42(default)
        set.seed(seed)
        
        ## define a list of genes-of-interest (GOIs)
        list_goi <- c("GPNMB","GCAR",
                      "CD8A","CD8B",
                      "MYC", #MYC(aka c-Myc)
                      "CD3D","CD3E","CD3G", #"CD3W"
                      "PDCD1","CD274" #CD274(aka PDL1)
        )
        
        ## plot dimensions ##
        ## set the "number of columns" of the [FeaturePlot]
        nGOI <- length(list_goi)
        if (nGOI<=4) { #1-4
          ncol_umap <- nGOI
        } else if (nGOI>4 & nGOI<7) { #5-6
          ncol_umap <- 3
        } else if (nGOI>6 & nGOI<9) { #7-8
          ncol_umap <- 4
        } else if (nGOI>8 & nGOI<11) { #9-10
          ncol_umap <- 5
        } else {
          print("[CHECK] 'nGOI' variable - to set the number of columns for the UMAP plots")
        }
        ## set the width & height of the plot
        if (nGOI<=4) {
          mywidth_umap <- ncol_umap*8.5
          myheight_umap <- 8.5
        } else {
          mywidth_umap <- ncol_umap*8.5
          myheight_umap <- round(nGOI/ncol_umap)*8.5
          # myheight_umap <- round(nGOI/7)*8.5
        }
        
        ## colour palette ##
        colors_ft2 <- c("#bdbdbd33", "#BC312A") #customized (light_grey | dark_red)
        
        ## plot 5: rnaUMAP plot - for GOI (quick version) ####
        fn_pdf_ft <- paste0(mydir_smp, "5_featurePlot_", smpid, "_nGOI",nGOI, ".pdf")
        # if (!file.exists(fn_pdf_ft)) {
        pdf(file = fn_pdf_ft, width = mywidth_umap, height = myheight_umap)
        print(
          FeaturePlot(sobj_flt_umap,
                      features = list_goi,
                      reduction = "umap",
                      order = T,
                      ncol = ncol_umap,
                      cols = colors_ft2
          )
        )
        dev.off()
        print(paste0("[DONE] save as PDF: ",fn_pdf_ft))
        # }
        
        
        #### gene lists ####
        ## geneList 0 | [GOI14] ####
        ## GOI ##
        list_goi <- c(
          "GPNMB","GCAR",
          "ASPSCR1","TFE3",
          "CD3D","CD3E","CD3G",
          "CD8A","CD8B",
          "CD4",
          "CD24","CD276",
          ## macrophage-related genes (from the chordoma kick-off meeting (Feb 23, 2024) ##
          "CD68","MRC1","CD163"
        )
        ngoi <- length(list_goi)
        print(paste0("--- nGOI: ",ngoi," ---"))
        
        ## geneList 1 | [mg_sc_syns] ####
        ## load the sarcoma-derived cell type marker genes ####
        ## - source: [Supplementary Table S6] from ModingLab's Nature Cancer 2024 paper (PMID: 38429415 | DOI: 10.1038/s43018-024-00743-y)
        fn_mg_sc_syns <- paste0(homedir,"data/0_gene_list/human_sarcoma/geneList_Moding2024-TableS6_scRNA_Broad-SYNS_CellType9.txt")
        mg_sc_syns <- read.table(file = fn_mg_sc_syns, sep = "\t", header = T, row.names = NULL, quote = NULL)
        dim(mg_sc_syns) #53 2
        head(mg_sc_syns,3)
        as.data.frame(table(mg_sc_syns$cell_type))
        # 1               B Cells    2
        # 2           CD4 T Cells    8
        # 3           CD8 T Cells    6
        # 4       Dendritic Cells    4
        # 5     Endothelial Cells    4
        # 6           Fibroblasts   14
        # 7            Mast Cells    5
        # 8 Monocytes/Macrophages    6
        # 9              NK Cells    4
        subset(mg_sc_syns, gene=="DUSP6") #checkpoint #Mast Cells
        
        mg_sc_syns_b <- subset(mg_sc_syns, cell_type=="B Cells")$gene
        mg_sc_syns_tcd4 <- subset(mg_sc_syns, cell_type=="CD4 T Cells")$gene
        mg_sc_syns_tcd8 <- subset(mg_sc_syns, cell_type=="CD8 T Cells")$gene
        mg_sc_syns_dc <- subset(mg_sc_syns, cell_type=="Dendritic Cells")$gene
        mg_sc_syns_endo <- subset(mg_sc_syns, cell_type=="Endothelial Cells")$gene
        mg_sc_syns_fibro <- subset(mg_sc_syns, cell_type=="Fibroblasts")$gene
        mg_sc_syns_mast <- subset(mg_sc_syns, cell_type=="Mast Cells")$gene
        mg_sc_syns_mono_mac <- subset(mg_sc_syns, cell_type=="Monocytes/Macrophages")$gene
        mg_sc_syns_nk <- subset(mg_sc_syns, cell_type=="NK Cells")$gene
        
        list_ctype_mg_sc_syns <- c("B","T_CD4","T_CD8",
                                   "DC","Endo","Fibro",
                                   "Mast","Mono_Mac",
                                   "NK"
        )
        length(list_ctype_mg_sc_syns) #9
        
        # ---
        ## geneList 2 | [mg_sc_sts2] ####
        ## load the sarcoma-derived cell type marker genes ####
        ## - source: [Supplementary Table S8] from ModingLab's Nature Cancer 2024 paper (PMID: 38429415 | DOI: 10.1038/s43018-024-00743-y)
        fn_mg_sc_sts2 <- paste0(homedir,"data/0_gene_list/human_sarcoma/geneList_Moding2024-TableS8_scRNA_Stanford-STS2_CellState13CellType5.txt")
        mg_sc_sts2 <- read.table(file = fn_mg_sc_sts2, sep = "\t", header = T, row.names = NULL, quote = NULL)
        dim(mg_sc_sts2) #53 2
        head(mg_sc_sts2,3)
        ## add a column: [manually_curated]
        ## - source: highlighted in bold font - according to the [Supplementary Table 8]
        mg_sc_sts2 <- mg_sc_sts2 %>%
          mutate(
            manually_curated=case_when(
              gene=="GZMA"|gene=="GZMK"|gene=="PDCD1"|gene=="TIGIT"|
                gene=="LILRB1"|gene=="TRAF2"|
                gene=="ARHGAP9"|gene=="MALAT1"|
                gene=="TRIM29"|gene=="SCUBE2"|
                gene=="COX7C"|gene=="NDUFB8"|
                gene=="NRP1"|gene=="G0S2"|
                gene=="CDC20"|gene=="ATP6V1D"|
                gene=="SMAD3"|gene=="EPB41L3"|
                gene=="RGS16"|gene=="TGFBI"|
                gene=="CLEC5A"|gene=="SPP1"|
                gene=="MRC1"|gene=="CD163"|gene=="IRF5"|gene=="MYO1F"|
                gene=="TPSAB1"|gene=="CD63"|gene=="ATF3"|gene=="TMSB4X" ~ "Y"
            )
          ) %>%
          as.data.frame()
        table(mg_sc_sts2$manually_curated) #33
        
        as.data.frame(table(mg_sc_sts2$cell_type))
        # 1           CD8 T cells   60 (nCellState=3)
        # 2 Epithelial-like cells   80 (nCellState=4)
        # 3 Fibroblast-like cells   40 (nCellState=2)
        # 4            Mast cells   40 (nCellState=2)
        # 5 Monocytes/macrophages   40 (nCellState=2)
        as.data.frame(table(mg_sc_sts2$cell_state_id))
        # 1            S01 CD8 T cells   20
        # 2  S01 Epithelial-like cells   20
        # 3  S01 Fibroblast-like cells   20
        # 4             S01 Mast cells   20
        # 5  S01 Monocytes/macrophages   20
        # 6            S02 CD8 T cells   20
        # 7  S02 Epithelial-like cells   20
        # 8  S02 Fibroblast-like cells   20
        # 9             S02 Mast cells   20
        # 10           S03 CD8 T cells   20
        # 11 S03 Epithelial-like cells   20
        # 12 S03 Monocytes/macrophages   20
        # 13 S04 Epithelial-like cells   20
        as.data.frame(table(mg_sc_sts2$cell_state_type,mg_sc_sts2$cell_type))
        # 1   S01           CD8 T cells   20
        # 2   S02           CD8 T cells   20
        # 3   S03           CD8 T cells   20
        # 5   S01 Epithelial-like cells   20
        # 6   S02 Epithelial-like cells   20
        # 7   S03 Epithelial-like cells   20
        # 8   S04 Epithelial-like cells   20
        # 9   S01 Fibroblast-like cells   20
        # 10  S02 Fibroblast-like cells   20
        # 13  S01            Mast cells   20
        # 14  S02            Mast cells   20
        # 17  S01 Monocytes/macrophages   20
        # 19  S03 Monocytes/macrophages   20
        as.data.frame(table(mg_sc_sts2$ecotype))
        # 1  SE1   60
        # 2  SE2  120
        # 3  SE3   40
        subset(mg_sc_sts2, gene=="CD109") #checkpoint #Monocytes/macrophages
        
        # mg_sc_sts2_se1 <- subset(mg_sc_sts2, ecotype=="SE1")$gene
        # mg_sc_sts2_se2 <- subset(mg_sc_sts2, ecotype=="SE2")$gene
        # mg_sc_sts2_se3 <- subset(mg_sc_sts2, ecotype=="SE3")$gene
        
        mg_sc_sts2_tcd8 <- subset(mg_sc_sts2, cell_type=="CD8 T cells")$gene
        mg_sc_sts2_tcd8_s01 <- subset(mg_sc_sts2, cell_type=="CD8 T cells" & cell_state_type=="S01")$gene
        mg_sc_sts2_tcd8_s02 <- subset(mg_sc_sts2, cell_type=="CD8 T cells" & cell_state_type=="S02")$gene
        mg_sc_sts2_tcd8_s03 <- subset(mg_sc_sts2, cell_type=="CD8 T cells" & cell_state_type=="S03")$gene
        
        mg_sc_sts2_epi <- subset(mg_sc_sts2, cell_type=="Epithelial-like cells")$gene
        mg_sc_sts2_epi_s01 <- subset(mg_sc_sts2, cell_type=="Epithelial-like cells" & cell_state_type=="S01")$gene
        mg_sc_sts2_epi_s02 <- subset(mg_sc_sts2, cell_type=="Epithelial-like cells" & cell_state_type=="S02")$gene
        mg_sc_sts2_epi_s03 <- subset(mg_sc_sts2, cell_type=="Epithelial-like cells" & cell_state_type=="S03")$gene
        mg_sc_sts2_epi_s04 <- subset(mg_sc_sts2, cell_type=="Epithelial-like cells" & cell_state_type=="S04")$gene
        
        mg_sc_sts2_fibro <- subset(mg_sc_sts2, cell_type=="Fibroblasts")$gene
        mg_sc_sts2_fibro_s01 <- subset(mg_sc_sts2, cell_type=="Fibroblasts" & cell_state_type=="S01")$gene
        mg_sc_sts2_fibro_s02 <- subset(mg_sc_sts2, cell_type=="Fibroblasts" & cell_state_type=="S02")$gene
        
        
        mg_sc_sts2_mast <- subset(mg_sc_sts2, cell_type=="Mast Cells")$gene
        mg_sc_sts2_mast_s01 <- subset(mg_sc_sts2, cell_type=="Mast Cells" & cell_state_type=="S01")$gene
        mg_sc_sts2_mast_s02 <- subset(mg_sc_sts2, cell_type=="Mast Cells" & cell_state_type=="S02")$gene
        
        mg_sc_sts2_mono_mac <- subset(mg_sc_sts2, cell_type=="Monocytes/Macrophages")$gene 
        mg_sc_sts2_mono_mac_s01 <- subset(mg_sc_sts2, cell_type=="Monocytes/Macrophages" & cell_state_type=="S01")$gene
        mg_sc_sts2_mono_mac_s03 <- subset(mg_sc_sts2, cell_type=="Monocytes/Macrophages" & cell_state_type=="S03")$gene
        
        ## define a list of the unique [cell_type]
        unique(mg_sc_sts2$cell_type)
        list_ctype_mg_sc_sts2 <- c("T_CD8",
                                   "Epi",
                                   "Fibro",                          
                                   "Mono_Mac",
                                   "Mast"
        )
        length(list_ctype_mg_sc_sts2) #5
        # ---
        
        ## geneList 3 | [mg_b_sts5] ####
        ## load the Sarcoma Ecotype-level cell state marker genes ####
        ## - source: [Supplementary Table S4] from ModingLab's Nature Cancer 2024 paper (PMID: 38429415 | DOI: 10.1038/s43018-024-00743-y)
        fn_mg_b_sts5 <- paste0(homedir,"data/0_gene_list/human_sarcoma/geneList_Moding2024-TableS4_bRNA_Stanford-STS5_CellState23CellType9.txt")
        mg_b_sts5 <- read.table(file = fn_mg_b_sts5, sep = "\t", header = T, row.names = NULL, quote = NULL)
        dim(mg_b_sts5) #3927 7
        head(mg_b_sts5,3)
        as.data.frame(table(mg_b_sts5$cell_type))
        # 1               B cells  297
        # 2           CD8 T cells  419
        # 3       Dendritic cells  277
        # 4     Endothelial cells  787
        # 5 Epithelial-like cells  903
        # 6 Fibroblast-like cells  597
        # 7            Mast cells   79
        # 8 Monocytes/macrophages  528
        # 9                  PMNs   40
        as.data.frame(table(mg_b_sts5$annotation))
        # 1                              Activated   73
        # 2                       Activated memory  260
        # 3                      Classical myeloid  173
        # 4                               Effector   13
        # 5                   Estrogen responsive   418
        # 6                    Exhausted cytotoxic  382
        # 7                         KRAS activated   91
        # 8                 KRAS activated/hypoxic  306
        # 9                           M1/M2 hybrid  340
        # 10             M2-like immunosuppressive  188
        # 11                                Mature  110
        # 12                    Mature neutrophils   19
        # 13                                Memory   52
        # 14                  MYC/MTORC1 activated  193
        # 15                                 Naive    9
        # 16                           Normal-like  516
        # 17 Oxidative phosphorylation upregulated  201
        # 18                        Pro-angiogenic  271
        # 19          Pro-inflammatory neutrophils   21
        # 20                   TNFA/KRAS activated  291
        subset(as.data.frame(table(mg_b_sts5$CS_SE,mg_b_sts5$annotation)),Freq>0)
        # 7              S01_Mast cells                             Activated   73
        # 24                S01_B cells                      Activated memory  260
        # 58        S02_Dendritic cells                     Classical myeloid  173
        # 80            S02_CD8 T cells                              Effector   13
        # 97  S01_Epithelial-like cells                  Estrogen responsive   418
        # 117           S01_CD8 T cells                   Exhausted cytotoxic  382
        # 159 S03_Epithelial-like cells                        KRAS activated   91
        # 176 S02_Fibroblast-like cells                KRAS activated/hypoxic  306
        # 206 S03_Monocytes/macrophages                          M1/M2 hybrid  340
        # 215 S01_Monocytes/macrophages             M2-like immunosuppressive   62
        # 224 S02_Monocytes/macrophages             M2-like immunosuppressive  126
        # 233       S01_Dendritic cells                                Mature  104
        # 246            S02_Mast cells                                Mature    6
        # 262                  S01_PMNs                    Mature neutrophils   19
        # 286               S02_B cells                                Memory   28
        # 296           S03_CD8 T cells                                Memory   24
        # 322 S04_Epithelial-like cells                  MYC/MTORC1 activated  193
        # 341               S03_B cells                                 Naive    9
        # 358     S02_Endothelial cells                           Normal-like  516
        # 382 S02_Epithelial-like cells Oxidative phosphorylation upregulated  201
        # 395     S01_Endothelial cells                        Pro-angiogenic  271
        # 432                  S02_PMNs          Pro-inflammatory neutrophils   21
        # 443 S01_Fibroblast-like cells                   TNFA/KRAS activated  291
        
        mg_b_sts5_texh <- subset(mg_b_sts5, annotation=="Exhausted cytotoxic")$gene #382
        
        ## geneList 4 | [mg_b_sts5_mrk] ####
        ## - source: [Supplementary Table S10] from ModingLab's Nature Cancer 2024 paper (PMID: 38429415 | DOI: 10.1038/s43018-024-00743-y)
        fn_mg_b_sts5_mrk <- paste0(homedir,"data/0_gene_list/human_sarcoma/geneList_Moding2024-TableS10_SarcEcotype-CellState23CellType9.txt")
        mg_b_sts5_mrk <- read.table(file = fn_mg_b_sts5_mrk, sep = "\t", header = T, row.names = NULL, quote = NULL)
        dim(mg_b_sts5_mrk) #52 7
        head(mg_b_sts5_mrk,3)
        ## fix the typos from the original publication: [GZMA],[GZMK]
        # mg_b_sts5_mrk$key_marker_gene[mg_b_sts5_mrk$key_marker_gene=="GRZMA"] <- "GZMA" #fixed
        # mg_b_sts5_mrk$key_marker_gene[mg_b_sts5_mrk$key_marker_gene=="GRZMK"] <- "GZMK" #fixed
        as.data.frame(table(mg_b_sts5_mrk$cell_type))
        # 1               B cells    7
        # 2           CD8 T cells    8
        # 3       Dendritic cells    4
        # 4     Endothelial cells    5
        # 5 Epithelial-like cells    8
        # 6 Fibroblast-like cells    4
        # 7            Mast cells    4
        # 8 Monocytes/macrophages    8
        # 9                  PMNs    4
        ## check the number of the [key_marker_gene] ####
        length(mg_b_sts5_mrk$key_marker_gene) #52
        length(unique(mg_b_sts5_mrk$key_marker_gene)) #46
        # ---
        # ===
        
        ## plot 5a: UMAP plots - using the marker gene lists ####
        list_mgset <- c(
          # "scRNA_sts2", #done #nUniqueGene=250 #nUniqueGene(manuallyCurated)=30
          ## "bRNA_sts5", #nUniqueGene=3507
          "bRNA_sts5_mrk", #nUniqueGene=46
          "scRNA_syns", #nUniqueGene=49
          "asps"
        )
        
        make_fplot <- T
        make_fplot_byCtype <- T
        
        ## make plots using the marker gene lists
        for (mgset in list_mgset) {
          print(paste0("-- mgset: [",mgset,"] ---"))
          
          if (mgset=="scRNA_sts2") {
            mg_sarc <- subset(mg_sc_sts2,manually_curated=="Y")
            mgeneset <- unique(mg_sarc$gene)
          } else if (mgset=="scRNA_sts2_all") {
            mg_sarc <- mg_sc_sts2
            mgeneset <- unique(mg_sarc$gene)
          } else if (mgset=="bRNA_sts5") {
            mg_sarc <- mg_b_sts5
            mgeneset <- unique(mg_sarc$gene)
          } else if (mgset=="bRNA_sts5_mrk") {
            mg_sarc <- mg_b_sts5_mrk
            # mgeneset <- unique(mg_b_sts5_mrk$key_marker_gene)
            which(colnames(mg_sarc)=="key_marker_gene") #6
            colnames(mg_sarc)[6] <- "gene"
            mgeneset <- unique(mg_sarc$gene)
          } else if (mgset=="scRNA_syns") {
            mg_sarc <- mg_sc_syns
            mgeneset <- unique(mg_sc_syns$gene)
          } else if (mgset=="asps") {
            mgeneset <- c(
              "ASPSCR1","TFE3",
              "GPNMB","GCAR"
            )
          }
          
          n_mgene <- length(mgeneset)
          print(paste0("--- mgset: [",mgset,"] | nGene: ",n_mgene," ---"))
          
          
          ## make feature plots using the marker gene lists ####
          if (isTRUE(make_fplot)) {
            ## feature plots: the expression levels of the marker genes ##
            # n_mgene <- length(mgeneset)
            if (n_mgene<7) {
              ncol_umap <- n_mgene
              mywidth_umap <- ncol_umap*4
              myheight_umap <- 4
            } else {
              ncol_umap <- 7
              mywidth_umap <- ncol_umap*4
              myheight_umap <- round(n_mgene/7)*4
              
              if (mgset=="scRNA_sts2") {
                myheight_umap <- myheight_umap+4
              }
            }
            
            if (sr_mode=="mrg") {
              myheight2_umap <- myheight_umap*2
            }
            
            ## featurePlot: by all marker genes for each cell type - orderT ##
            set.seed(seed)
            fn_pdf_ft <- paste0(mydir_smp, "5a_featurePlot_", smpid,"_",mgset,"_n", n_mgene, "_orderT.pdf")
            if (!file.exists(fn_pdf_ft)) {
              pdf(file = fn_pdf_ft, width = mywidth_umap, height = myheight_umap)
              print(
                FeaturePlot(sobj_flt_umap, 
                            features = mgeneset, 
                            reduction = "umap", 
                            order = T, 
                            ncol = ncol_umap, 
                            cols = colors_ft2
                )
              )
              dev.off()
              print(paste0("[DONE] save as PDF: ",fn_pdf_ft))
            }
            
            ## featurePlot: by all marker genes for each cell type - orderT (split by each sample) ##
            if (sr_mode=="mrg") {
              if (n_mgene <10) {
                fn_pdf_ft <- paste0(mydir_smp, "5a_featurePlot_", smpid,"_",mgset,"_n", n_mgene, "_orderT","_bySmp",".pdf")
                if (!file.exists(fn_pdf_ft)) {
                  pdf(file = fn_pdf_ft, width = mywidth_umap, height = myheight2_umap)
                  print(
                    FeaturePlot(sobj_flt_umap, 
                                features = mgeneset, 
                                reduction = "umap", 
                                order = T, 
                                ncol = ncol_umap, 
                                cols = colors_ft2,
                                split.by = "orig.ident"
                    )
                  )
                  dev.off()
                  print(paste0("[DONE] save as PDF: ",fn_pdf_ft))
                }
              } else {
                print(paste0("[SKIP] FeaturePlot - due to large gene list (nCType=",n_mgene,")"))
              } #closer: if (n_mgene_byCType <10) {
            } #closer: if (sr_mode=="mrg") {
            
          } #closer: if (isTRUE(make_fplot)) {
          
          
          ## make feature plots by cell type - using the marker genes ####
          if (isTRUE(make_fplot_byCtype)) {
            if (mgset=="scRNA_sts2") {
              list_ctype_mg <- c(
                "T_CD8", 
                "Epi", 
                "Fibro", 
                "Mono_Mac", 
                "Mast" 
              )
              # mgeneset <- unique(mg_sc_sts2$gene)
            } else if (mgset=="scRNA_syns") {
              list_ctype_mg <- c(
                # "B" #, #note: no expression of: [CD79A],[MS4A1],[GZMB]
                "T_CD4",#note; no expression of: [CTLA4]
                "T_CD8",
                "DC",
                "Endo",
                "Fibro",
                "Mast",
                "Mono_Mac",
                "NK"
              )
            } else if (mgset=="bRNA_sts5_mrk") {
              list_ctype_mg <- c(
                "B",
                "T_CD8",
                "DC",
                "Endo",
                "Epi",
                "Fibro",
                "Mast",
                "Mono_Mac",
                "PMN"
              )
            } else if (mgset=="asps") {
              list_ctype_mg <- c(
                "ASPS"
              )
            }
            
            for (ctype in list_ctype_mg) {
              print(paste0("--- ctype: [",ctype,"] ---"))
              
              if (ctype=="B") {
                mg_sarc_byCType <- unique(subset(mg_sarc, cell_type=="B Cells"|cell_type=="B cells")$gene)
              } else if (ctype=="T_CD4") {
                mg_sarc_byCType <- unique(subset(mg_sarc, cell_type=="CD4 T Cells"|cell_type=="CD4 T cells")$gene)
              } else if (ctype=="T_CD8") {
                mg_sarc_byCType <- unique(subset(mg_sarc, cell_type=="CD8 T Cells"|cell_type=="CD8 T cells")$gene)
              } else if (ctype=="DC") {
                mg_sarc_byCType <- unique(subset(mg_sarc, cell_type=="Dendritic Cells"|cell_type=="Dendritic cells")$gene)
              } else if (ctype=="Endo") {
                mg_sarc_byCType <- unique(subset(mg_sarc, cell_type=="Endothelial Cells"|cell_type=="Endothelial cells")$gene)
              } else if (ctype=="Fibro") {
                mg_sarc_byCType <- unique(subset(mg_sarc, cell_type=="Fibroblasts"|cell_type=="Fibroblast-like cells")$gene)
              } else if (ctype=="Mast") {
                mg_sarc_byCType <- unique(subset(mg_sarc, cell_type=="Mast Cells"|cell_type=="Mast cells")$gene)
              } else if (ctype=="Mono_Mac") {
                mg_sarc_byCType <- unique(subset(mg_sarc, cell_type=="Monocytes/Macrophages"|cell_type=="Monocytes/macrophages")$gene)
              } else if (ctype=="NK") {
                mg_sarc_byCType <- unique(subset(mg_sarc, cell_type=="NK Cells"|cell_type=="NK cells")$gene)
              } else if (ctype=="Epi") {
                mg_sarc_byCType <- unique(subset(mg_sarc, cell_type=="Epithelial-like cells")$gene)
              } else if (ctype=="PMN") {
                mg_sarc_byCType <- unique(subset(mg_sarc, cell_type=="PMNs")$gene)
              } else if (ctype=="ASPS") {
                mg_sarc_byCType <- c("ASPSCR1","TFE3","GPNMB","GCAR")
              }
              
              n_mgene_byCType <- length(mg_sarc_byCType)
              
              ## colour palette ##
              colors_ft2 <- c("#bdbdbd33", "#BC312A") #customized (light_grey | dark_red)
              
              if (length(mg_sarc_byCType)<7) {
                ncol_umap <- length(mg_sarc_byCType)
                mywidth_umap <- ncol_umap*8.5
                myheight_umap <- 8.5
              } else {
                ncol_umap <- 7
                mywidth_umap <- ncol_umap*8.5
                myheight_umap <- round(n_mgene_byCType/7)*8.5
                
                if (mgset=="scRNA_sts2" & ctype=="T_CD8") {
                  myheight_umap <- myheight_umap*2
                } else if (mgset=="scRNA_syns" & ctype=="T_CD4") {
                  myheight_umap <- myheight_umap*2
                }
              }
              
              if (sr_mode=="mrg") {
                myheight2_umap <- myheight_umap*2
              }
              
              ## featurePlot: by marker gene list for each cell type - orderT ##
              fn_pdf_ft <- paste0(mydir_smp, "5b_featurePlot_", smpid,"_",mgset,"-", ctype,"_n",n_mgene_byCType, "_orderT.pdf")
              if (!file.exists(fn_pdf_ft)) {
                pdf(file = fn_pdf_ft, width = mywidth_umap, height = myheight_umap)
                print(
                  FeaturePlot(sobj_flt_umap, 
                              features = mg_sarc_byCType, 
                              reduction = "umap", 
                              order = T, 
                              ncol = ncol_umap, 
                              cols = colors_ft2
                  )
                )
                dev.off()
                print(paste0("[DONE] save as PDF: ",fn_pdf_ft))
              }
              
              ## featurePlot: by marker gene list for each cell type - orderT (split by each sample) ##
              if (sr_mode=="mrg") {
                if (n_mgene_byCType <10) {
                  fn_pdf_ft <- paste0(mydir_smp, "5b_featurePlot_", smpid,"_",mgset,"-", ctype,"_n",n_mgene_byCType, "_orderT","_bySmp",".pdf")
                  if (!file.exists(fn_pdf_ft)) {
                    pdf(file = fn_pdf_ft, width = mywidth_umap, height = myheight2_umap)
                    print(
                      FeaturePlot(sobj_flt_umap, 
                                  features = mg_sarc_byCType, 
                                  reduction = "umap", 
                                  order = T, 
                                  ncol = ncol_umap, 
                                  cols = colors_ft2,
                                  split.by = "orig.ident"
                      )
                    )
                    dev.off()
                    print(paste0("[DONE] save as PDF: ",fn_pdf_ft))
                  }
                } else {
                  print(paste0("[SKIP] FeaturePlot - due to large gene list (nCType=",n_mgene_byCType,")"))
                } #closer: if (n_mgene_byCType <10) {
              } #closer: if (sr_mode=="mrg") {
              
            } #closer:  for (ctype in list_ctype_mg) {
          } #closer:  if (isTRUE(make_fplot_byGrp)) {
        } #closer:  for (mgset in list_mgset) {
        
        
      } #closer: if (isTRUE(sr_mk_plt_5)) {
    } #closer: for (smpid in list_smpid) { 
    print(":")
  } #closer: if (isTRUE(sr_09a_afterUMAP)) {
  
  ## step 10 | Finding differentially expressed features (cluster biomarkers) ####
  ## - source: "https://satijalab.org/seurat/articles/de_vignette"
  if (isTRUE(sr_10_de)) {
    
    ## install & load the "recommended" package for more efficient calculation of the "Wilcoxon Rank Sum Test" #skip-for-now (installed using "conda")
    #   # For a (much!) faster implementation of the Wilcoxon Rank Sum Test,
    #   # (default method for FindMarkers) please install the presto package
    #   # --------------------------------------------
    
    install.packages('devtools')
    devtools::install_github('immunogenomics/presto', dependencies = T, force = T)
    
    #   # --------------------------------------------
    #   # After installation of presto, Seurat will automatically use the more 
    #   # efficient implementation (no further action necessary).    
    library(presto)
    
    
    for (smpid in list_smpid) {
      print(paste0("--- STEP 10 | find DE features (cluster biomarkers) | smp_id: [",smpid,"]"," (",which(list_smpid==smpid)," /",length(list_smpid),")"," ---"))
      
      ## set the Seurat Object for filtering step
      if (!isTRUE(tcell_only)) {
        sr_dtype <- sr_dtype #raw | filtered
      } else if (isTRUE(tcell_only)) {
        sr_dtype <- paste0(sr_dtype,"_vdj") #raw | filtered
      }
      
      ## load the RDS: [sobj_flt_umap] - after UMAP      
      fn_rds_sobj_flt_umap <- paste0(dir_res_sr,"sobj_flt_umap_",prj,"_",smpid,"_",sr_dtype,"_res",res_clst,".rds")
      if (file.exists(fn_rds_sobj_flt_umap)) {
        sobj_flt_umap <- readRDS(fn_rds_sobj_flt_umap)
        print(paste0("[DONE] load the RDS: ",fn_rds_sobj_flt_umap))
      } else {
        print(paste0("[CHECK] RDS file does not exist: ",fn_rds_sobj_flt_umap))
      }
      
      ## check the number of the Seurat Clusters
      print(sort(unique(sobj_flt_umap@meta.data$seurat_clusters)))
      nClst <- length(unique(sobj_flt_umap@meta.data$seurat_clusters))
      print(nClst)
      
      ## example ##
      # ## find all markers of [cluster 2]
      # cluster2.markers <- FindMarkers(sobj_flt_umap, ident.1 = 2)
      # head(cluster2.markers, n = 5)
      
      ## find all markers distinguishing [cluster 5] from [clusters 0 and 3]
      # cluster5.markers <- FindMarkers(sobj_flt_umap, ident.1 = 5, ident.2 = c(0, 3))
      # head(cluster5.markers, n = 5)
      
      ## Seurat - FindAllMarkers ####
      if (isTRUE(sr_10_de_FindAllMarkers)) {
        
        ## find markers for every cluster compared to all remaining cells, report only the positive ones
        print(paste0("[BEGIN] run 'FindAllMarkers': ", Sys.time()))
        sobj_flt_umap.markers <- FindAllMarkers(sobj_flt_umap, only.pos = FALSE)
        print(paste0("[DONE] run 'FindAllMarkers': ", Sys.time()))
        
        print(paste0("[BEGIN] filter 'avg_log2FC >1' by each cluster: ", Sys.time()))
        sobj_flt_umap.markers %>%
          group_by(cluster) %>%
          dplyr::filter(avg_log2FC > 1)
        print(paste0("[DONE] filter 'avg_log2FC >1' by each cluster: ", Sys.time()))
        
        ## convert the [FindAllMarkers] output as data frame
        df_markers_byClst <- as.data.frame(sobj_flt_umap.markers)
        ## change the order of columns
        print(colnames(df_markers_byClst)) #checkpoint
        df_markers_byClst <- subset(df_markers_byClst, select=c(cluster,gene,p_val, 
                                                                avg_log2FC, pct.1, pct.2, p_val_adj
        ))
        print(colnames(df_markers_byClst)) #checkpoint
        
        ## save as TXT: [df_markers_byClst]
        fn_txt_markers_byClst <- paste0(dir_res_sr,"markers_byClst_",prj,"_",sr_dtype,"_res",res_clst,"_nClst",nClst,".txt")
        print(paste0("[BEGIN] save 'df_markers_byClst' as TXT: ", Sys.time()))
        write.table(x = df_markers_byClst, file = fn_txt_markers_byClst, sep = "\t", col.names = T, row.names = F, quote = F)
        print(paste0("[DONE] save 'df_markers_byClst' as TXT: ", Sys.time()))
        
        print(paste0("[DONE] save as TXT: ",fn_txt_markers_byClst))
        
      } #closer: if (isTRUE(sr_10_de_FindAllMarkers)) {
      
      
      ## Seurat has several tests for differential expression which can be set with the test.use parameter (see our DE vignette for details). 
      ## For example, the ROC test returns the ‘classification power’ for any individual marker (ranging from 0 - random, to 1 - perfect).
      if (isTRUE(sr_10_de_FindMarkers)) {
        cluster0.markers <- FindMarkers(sobj_flt_umap, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
      } #closer: if (isTRUE(sr_10_de_FindMarkers)) {
      
      
      if (isTRUE(sr_10_de_vis)) {
        ## violin plot - using [data]
        VlnPlot(sobj_flt_umap, features = c("MS4A1", "CD79A"))
        
        ## you can plot raw counts as well
        ## violin plot - using [counts]
        VlnPlot(sobj_flt_umap, 
                features = c("NKG7", "PF4"), 
                slot = "counts", log = TRUE
        )
        
        ## UMAP-based plot:
        FeaturePlot(sobj_flt_umap, 
                    features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A")
        )
      } #closer: if (isTRUE(sr_10_de_FindAllMarkers)) {
      
      ## DoHeatmap() generates an expression heatmap for given cells and features. In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.
      sobj_flt_umap.markers %>%
        group_by(cluster) %>%
        dplyr::filter(avg_log2FC > 1) %>%
        slice_head(n = 10) %>%
        ungroup() -> top10
      DoHeatmap(sobj_flt_umap, features = top10$gene) + NoLegend()
    } #closer: for (smpid in list_smpid) {
    print(":")
  } #closer: if (isTRUE(sr_10_de)) {
  
  
  ## step 11 | Assigning cell type identity to clusters ####
  # Fortunately in the case of the example dataset, we can use canonical markers to easily match the unbiased clustering to known cell types:
  
  #   Cluster ID	Markers	Cell Type
  # 0	IL7R, CCR7	Naive CD4+ T
  # 1	CD14, LYZ	CD14+ Mono
  # 2	IL7R, S100A4	Memory CD4+
  # 3	MS4A1	B
  # 4	CD8A	CD8+ T
  # 5	FCGR3A, MS4A7	FCGR3A+ Mono
  # 6	GNLY, NKG7	NK
  # 7	FCER1A, CST3	DC
  # 8	PPBP	Platelet
  if (isTRUE(sr_11_ctype)) {
    for (smpid in list_smpid) {
      print(paste0("--- STEP 11 | assign cell type identity to clusters | smp_id: [",smpid,"]"," (",which(list_smpid==smpid)," /",length(list_smpid),")"," ---"))
      
      ## load the RDS: [sobj_flt_umap]
      fn_rds_sobj_flt_umap <- paste0(dir_res_sr,"sobj_flt_umap_",prj,"_",smpid,"_",sr_dtype,".rds")
      if (file.exists(fn_rds_sobj_flt_umap)) {
        sobj_flt_umap <- readRDS(fn_rds_sobj_flt_umap)
        print(paste0("[DONE] load the RDS: ",fn_rds_sobj_flt_umap))
      } else {
        print(paste0("[CHECK] RDS file does not exist: ",fn_rds_sobj_flt_umap))
      }
      
      new.cluster.ids <- c(
        "Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
        "NK", "DC", "Platelet"
      )
      names(new.cluster.ids) <- levels(sobj_flt_umap)
      ## rename cell clusters based on the marker genes
      sobj_flt_umap <- RenameIdents(sobj_flt_umap, new.cluster.ids)
      
      ## UMAP-based plot:
      print(
        DimPlot(sobj_flt_umap, 
                reduction = "umap", label = TRUE, pt.size = 0.5) + 
          NoLegend()
      )
      
      ## UMAP-based plot - w/ details:
      # library(ggplot2)
      plot <- 
        DimPlot(sobj_flt_umap, 
                reduction = "umap", 
                label = TRUE, label.size = 4.5
        ) + 
        xlab("UMAP 1") + 
        ylab("UMAP 2") +
        theme(axis.title = element_text(size = 18), 
              legend.text = element_text(size = 18)
        ) + 
        guides(colour = guide_legend(override.aes = list(size = 10)))
      ## save as PDF:
      fn_pdf <- paste0("UMAP_",prjid,"_",smpid,".pdf")
      ggsave(filename = fn_pdf, height = 7, width = 12, plot = plot, quality = 50)
      
      ## save as RDS: [sobj_flt_umap] - w/ new cell type identities ####
      fn_rds_sobj_flt_umap_ctype <- paste0(dir_res_sr,"sobj_flt_umap_",prj,"_",smpid,"_",sr_dtype,"_wCellType",".rds")
      saveRDS(object = sobj_flt_umap, file = fn_rds_sobj_flt_umap_ctype)
      print(paste0("[DONE] save as RDS:: ",fn_rds_sobj_flt_umap_ctype))
      
    } #closer: for (smpid in list_smpid) {
    print(":")
  } #closer: if (isTRUE(sr_11_ctype)) {
  
} #closer: if (isTRUE(run_seurat)) {


#### Seurat: merge Seurat Objects ####
if (isTRUE(sr_mrg_sobj)) {
  # Install devtools from CRAN
  # install.packages("devtools", force=T, dependencies=T, quiet=F, repos="https://mirror.csclub.uwaterloo.ca/CRAN/")
  # Or the development version from GitHub:
  # devtools::install_github("r-lib/devtools")
  # library(devtools)
  
  # devtools::install_github('satijalab/seurat-data')
  # library(SeuratData)
  
  if (prj=="GCAR1_PT01_GEX_VDJ" | prj=="GCAR1_PT01_GEX" | prj=="GCAR1_PT01_VDJ") {
    prjid <- "GCAR1_PT01"
    # sp <- "hs38_GCAR"
    
    ## set a file path for each Seurat Object
    fn_sobj1 <- paste0(dir_res_sr,"sobj_GCAR1_PT01_GEX_Enriched_raw.rds")
    fn_sobj2 <- paste0(dir_res_sr,"sobj_GCAR1_PT01_GEX_Harvest_raw.rds")
    fn_sobj3 <- paste0(dir_res_sr,"sobj_GCAR1_PT01_GEX_D08_raw.rds")
    fn_sobj4 <- paste0(dir_res_sr,"sobj_GCAR1_PT01_GEX_D15_raw.rds")
    fn_sobj5 <- paste0(dir_res_sr,"sobj_GCAR1_PT01_GEX_D22_raw.rds")
    fn_sobj6 <- paste0(dir_res_sr,"sobj_GCAR1_PT01_GEX_D28_raw.rds")
    fn_sobj7 <- paste0(dir_res_sr,"sobj_GCAR1_PT01_GEX_D42_raw.rds")
    
    ## load individual Seurat Objects to merge
    sobj1 <- readRDS(fn_sobj1)
    print(paste0("[DONE] read the RDS: ",fn_sobj1))
    sobj2 <- readRDS(fn_sobj2)
    print(paste0("[DONE] read the RDS: ",fn_sobj2))
    sobj3 <- readRDS(fn_sobj3)
    print(paste0("[DONE] read the RDS: ",fn_sobj3))
    sobj4 <- readRDS(fn_sobj4)
    print(paste0("[DONE] read the RDS: ",fn_sobj4))
    sobj5 <- readRDS(fn_sobj5)
    print(paste0("[DONE] read the RDS: ",fn_sobj5))
    sobj6 <- readRDS(fn_sobj6)
    print(paste0("[DONE] read the RDS: ",fn_sobj6))
    sobj7 <- readRDS(fn_sobj7) #exclude-for-now
    print(paste0("[DONE] read the RDS: ",fn_sobj7))
    
    ## make a list of Seurat Objects
    ## [entire list] ##
    list_sobj7 <- c(sobj1,sobj2,
                    sobj3,sobj4,sobj5,sobj6,sobj7
    )
    
    
    ## make a list of [smpid]
    print(head(mdata,3))
    list_smpid <- unique(mdata$smp_id)
    print(list_smpid)
    
    ## count the number of samples included in the [list_smpid]
    nSmp <- length(list_smpid)
    
    ## merge Seurat Objects into one Seurat Object
    if (nSmp==7) {
      mrg_sobj <- merge(sobj1, y=c(sobj2,sobj3,sobj4,sobj5,sobj6,sobj7), 
                        add.cell.ids=c("Enriched","Harvest","D08","D15","D22","D28","D42"),
                        project=prj
      )
      # mrg_sobj <- merge(sobj1, y=list_sobj6, 
      #                   add.cell.ids=list_smpid,
      #                   project=prj
      # )
      
    } else if (nSmp==6) {
      mrg_sobj <- merge(sobj1, y=c(sobj2,sobj3,sobj4,sobj5,sobj6), 
                        add.cell.ids=c("Enriched","Harvest","D08","D15","D22","D28"),
                        project=prj
      )
    }
    print(head(colnames(mrg_sobj))) #checkpoint
    print(tail(colnames(mrg_sobj))) #checkpoint
    
    ## check the unique samples included in the merged SeuratObject
    unique(sapply(X = strsplit(colnames(mrg_sobj), split = "_"), FUN = "[", 1)) #checkpoint
    ## check the number of cells from each sample
    table(mrg_sobj$orig.ident)
    
    ## save the merged Seurat Object as RDS
    # "sobj_GCAR1_PT01_GEX_mrg_raw.rds"
    fn_rds_mrg_sobj <- paste0(dir_res_sr,"sobj_",prj,"_","mrg",nSmp,"_",sr_dtype,".rds")
    saveRDS(object = mrg_sobj, file = fn_rds_mrg_sobj)
    print(paste0("[DONE] save 'mrg_sobj' as RDS: ",fn_rds_mrg_sobj))
  } #closer: if (prj=="GCAR1_PT01_GEX_VDJ" | prj=="GCAR1_PT01_GEX" | prj=="GCAR1_PT01_VDJ") {
  
} #closer: if (isTRUE(sr_mrg_sobj)) {


# ===
#### metadata from Seurat Object ####
## extract the metadata from the "merged" Seurat Object ####
if (isTRUE(ext_sobj_mdata)) {
  if (isTRUE(ext_sobj_mdata_flt)) {
    list_smpid <- unique(mdata$smp_id)
  } else {
    smpid <- "mrg7" #"mrg" #"mrg7"
    list_smpid <- c(smpid)
  } #closer: if (isTRUE(ext_sobj_mdata_flt)) {
  print(paste0("--- list_smpid: "))
  print(list_smpid)
  print(":")
  
  if (!isTRUE(tcell_only)) {
    sr_dtype <- sr_dtype #raw | filtered
  } else if (isTRUE(tcell_only)) {
    sr_dtype <- paste0(sr_dtype,"_vdj") #raw | filtered
  }
  
  for (smpid in list_smpid) {
    print(paste0("--- smpid: [ ",smpid," ] ---"))
    
    if (isTRUE(ext_sobj_mdata_flt)) {
      fn_rds_mrg_sobj <- paste0(dir_res_sr,"0_rds__","GCAR1_PT01_GEX/","sobj_flt_",prj,"_",smpid,"_",sr_dtype,".rds")
    } else {
      fn_rds_mrg_sobj <- paste0(dir_res_sr,"0_rds__","GCAR1_PT01_GEX/","sobj_flt_umap_",prj,"_",smpid,"_",sr_dtype,"_res",res_clst,".rds")
    } #closer: if (isTRUE(ext_sobj_mdata_flt)) {
    
    if (file.exists(fn_rds_mrg_sobj)) {
      mrg_sobj <- readRDS(fn_rds_mrg_sobj)
      print(paste0("[DONE] load the RDS: ",fn_rds_mrg_sobj))
    } else {
      print(paste0("[CHECK] RDS file does not exist: ",fn_rds_mrg_sobj))
    }
    
    ## extract the metadata as "data frame" format
    md_mrg_sobj <- as.data.frame(mrg_sobj@meta.data)
    print(dim(md_mrg_sobj))
    
    ## save as TXT: [md_mrg_sobj]
    if (isTRUE(ext_sobj_mdata_flt)) {
      fn_txt_md <- paste0(dir_res_sr,"metadata_",smpid,"_sobj_flt_",prj,".txt")
    } else {
      fn_txt_md <- paste0(dir_res_sr,"metadata_",smpid,"_sobj_",prj,"_res",res_clst,".txt")
    } #closer: if (isTRUE(ext_sobj_mdata_flt)) {
    write.table(x = md_mrg_sobj, file = fn_txt_md, sep = "\t", col.names = T, row.names = T, quote = F)
    if (file.exists(fn_txt_md)) {
      print(paste0("[DONE] save as TXT: ",fn_txt_md))
    } else {
      print(paste0("[CHECK] TXT did not create successfully: ",fn_txt_md))
    }
    
  } #closer: for (smpid in list_smpid) {
  
} #closer: if (isTRUE(ext_sobj_mdata)) {
# ===


# ===
#### add metadata to the Seurat Object ####
if (isTRUE(add_sobj_mdata)) {
  ## add metadata from [VDJ] data ####
  if (isTRUE(add_sobj_mdata_vdj)) {
    ## load the "source" metadata 
    fn_md_vdj <- paste0(dir_res,"scRepertoire/","metadata_clones_GCAR1_PT01_VDJ_nSmp7.txt")
    md_vdj <- read.table(file = fn_md_vdj, sep = "\t", header = T, row.names = NULL, quote = NULL)
    print(dim(md_vdj)) #81520 10
    print(head(md_vdj,3))
    print(tail(md_vdj,3))
    ## add columns: [cell_VDJ],[sample_new_VDJ]
    md_vdj$cell_VDJ <- "VDJ"
    md_vdj$sample_new_VDJ <- paste(md_vdj$sample_new, md_vdj$cell_VDJ, sep = "_")
    print(head(md_vdj,3))
    print(dim(md_vdj)) #81520 12
    
    ## select columns to add to the existing Seurat object
    md_vdj_sel <- subset(md_vdj, select=c(sample_new,cell_barcode,sample_new_VDJ,cell_VDJ))
    ## set row names
    rownames(md_vdj_sel) <- md_vdj$cell_barcode
    print(head(md_vdj_sel,3))
    print(dim(md_vdj_sel)) #81521 2
    # ---
    
    ## load the Seurat UMAP object
    fn_rds_mrg_sobj_umap2.4 <- paste0(dir_res_sr,"sobj_flt_umap_GCAR1_PT01_GEX_mrg_raw_res2.4.rds")
    if (file.exists(fn_rds_mrg_sobj_umap2.4)) {
      mrg_sobj_umap2.4 <- readRDS(fn_rds_mrg_sobj_umap2.4)
      print(paste0("[DONE] load the RDS: ",fn_rds_mrg_sobj_umap2.4))
    } else {
      print(paste0("[CHECK] RDS file does not exist: ",fn_rds_mrg_sobj_umap2.4))
    }
    
    ## define the column name from the "source" Seurat Object
    # cname_clst <- "cell_VDJ"
    cname_clst <- colnames(md_vdj_sel)
    
    ## checkpoint:
    print("[BEFORE] add the metadata: ")
    print(head(mrg_sobj_umap2.4@meta.data))
    print(tail(mrg_sobj_umap2.4@meta.data))
    
    ## add the extracted metadata column to the newer version of the UMAP Seurat Object
    mrg_sobj_umap2.4 <- AddMetaData(object = mrg_sobj_umap2.4,
                                    metadata = md_vdj_sel,
                                    col.name = cname_clst
    )
    
    ## checkpoint:
    print("[AFTER] add the metadata: ")
    print(head(mrg_sobj_umap2.4@meta.data))
    print(tail(mrg_sobj_umap2.4@meta.data))
    
    ## update the RDS file: [mrg_sobj_umap2.4]
    print(paste0("[START] update 'sobj' RDS | ", Sys.time()))
    saveRDS(object = mrg_sobj_umap2.4, file = fn_rds_mrg_sobj_umap2.4)
    print(paste0("[DONE] update the RDS: ",fn_rds_mrg_sobj_umap2.4))
    print(paste0("[END] update 'sobj' RDS | ", Sys.time()))
    
  } #closer: if (isTRUE(add_sobj_mdata_vdj)) {
  # ---
  
  ## add metadata from [VD]] data ####
  if (isTRUE(add_sobj_mdata_trb)) {
    ## load the "source" metadata 
    fn_md_trb <- paste0(homedir,"collaboration/","RCCI__GCAR1","/03b_result_VDJ/","clones_GCAR1_PT01_VDJ_all__nSmp7","_updated",".csv")
    md_trb <- read.table(file = fn_md_trb, sep = ",", header = T, row.names = NULL, quote = NULL)
    dim(md_trb) #81645 13->22
    print(head(md_trb,3))
    print(tail(md_trb,3))
    
    # ---
    ## column-of-interest: [combLen_chain1_chain2_gene]
    # ## step 2 | add a column: [select_combLen_chain1_chain2_gene] - based on the subset ####
    # # head(cln_all_up,1)
    # cln_all_up <- cln_all_up %>%
    #   mutate(
    #     select_combLen_chain1_chain2_gene=case_when(
    #       combLen_chain1_chain2_gene=="1--0" | combLen_chain1_chain2_gene=="1--2" | combLen_chain1_chain2_gene=="2--2" ~ "N_exclude",
    #       combLen_chain1_chain2_gene!="1--0" & combLen_chain1_chain2_gene!="1--2" & combLen_chain1_chain2_gene!="2--2" ~ "Y_select"
    #     )
    #   ) %>%
    #   as.data.frame()
    # head(cln_all_up,3)
    # dim(cln_all_up) #81645 21
    # ---
    
    # ## add columns: [cell_VDJ],[sample_new_VDJ]
    # md_trb$cell_VDJ <- "VDJ"
    # md_trb$sample_new_VDJ <- paste(md_trb$sample_new, md_trb$cell_VDJ, sep = "_")
    # print(head(md_trb,3))
    # print(dim(md_trb)) #81520 12
    
    ## select columns to add to the existing Seurat object
    md_trb_sel <- subset(md_trb, select=c(sample_new,cell_barcode,combLen_chain1_chain2_gene,select_combLen_chain1_chain2_gene))
    ## set row names
    rownames(md_trb_sel) <- md_trb$cell_barcode
    print(head(md_trb_sel,3))
    print(dim(md_trb_sel)) #81521 2
    # ---
    
    ## load the Seurat UMAP object
    fn_rds_mrg_sobj_umap2.4 <- paste0(dir_res_sr,"0_rds__GCAR1_PT01_GEX/","sobj_flt_umap_GCAR1_PT01_GEX-subsetVDJ_mrg7_raw_vdj_res2.4.rds") #GEX data w/ VDJ+ cells only
    if (file.exists(fn_rds_mrg_sobj_umap2.4)) {
      mrg_sobj_umap2.4 <- readRDS(fn_rds_mrg_sobj_umap2.4)
      print(paste0("[DONE] load the RDS: ",fn_rds_mrg_sobj_umap2.4))
    } else {
      print(paste0("[CHECK] RDS file does not exist: ",fn_rds_mrg_sobj_umap2.4))
    }
    
    ## define the column name from the "source" Seurat Object
    # cname_clst <- "cell_VDJ"
    cname_clst <- colnames(md_trb_sel)
    
    ## checkpoint:
    print("[BEFORE] add the metadata: ")
    print(head(mrg_sobj_umap2.4@meta.data))
    print(tail(mrg_sobj_umap2.4@meta.data))
    
    ## add the extracted metadata column to the newer version of the UMAP Seurat Object
    mrg_sobj_umap2.4 <- AddMetaData(object = mrg_sobj_umap2.4,
                                    metadata = md_trb_sel,
                                    col.name = cname_clst
    )
    
    ## remove the "single-category" variable metadata column from the existing metadata
    ## ( as this column is not necessary when converting from [sobj] to [h5ad] )
    ## - "Dropping single category variables:[cell_VDJ]"
    mrg_sobj_umap2.4$cell_VDJ <- NULL
    mrg_sobj_umap2.4[['cell_VDJ']] <- NULL
    
    ## checkpoint:
    print("[AFTER] add the metadata: ")
    print(head(mrg_sobj_umap2.4@meta.data))
    print(tail(mrg_sobj_umap2.4@meta.data))
    
    ## update the RDS file: [mrg_sobj_umap2.4]
    print(paste0("[START] update 'sobj' RDS | ", Sys.time()))
    saveRDS(object = mrg_sobj_umap2.4, file = fn_rds_mrg_sobj_umap2.4)
    print(paste0("[DONE] update the RDS: ",fn_rds_mrg_sobj_umap2.4))
    print(paste0("[END] update 'sobj' RDS | ", Sys.time()))
    
  } #closer: if (isTRUE(add_sobj_mdata_trb)) {
  # ---
  
} #closer: if (isTRUE(add_sobj_mdata)) {
# ===


# ===
########################################
## convert the SeuratObject to h5ad ####
########################################
if (isTRUE(sceasy_cnvt_to_h5ad)) {
  
  ## set [sr_mode] variable
  sr_mode <- "mrg"
  
  #### packages ####
  mypackages_sceasy <- c(
    "BiocManager"
    ,"sceasy"
    ,"LoomExperiment", "SingleCellExperiment"
    ,"GenomicRanges"
  )
  
  ## [sceasy] ##
  devtools::install_github("cellgeni/sceasy", dependencies = T, force = T)
  ## dependent packages ##
  BiocManager::install(c("LoomExperiment", "SingleCellExperiment"), force = T)
  BiocManager::install("GenomicRanges", force = T)
  
  lapply(mypackages_sceasy, library, character.only = TRUE)
  
  # use_condaenv('EnvironmentName') #template
  use_condaenv('sceasy') #template
  
  ## Optionally, if you plan to convert between loom and anndata, please also ensure that the loompy package is installed:
  ## - source: [https://github.com/cellgeni/sceasy]
  # loompy <- reticulate::import('loompy')
  # ---
  
  ## define a list of sample IDs
  if (sr_mode=="each") {
    list_smpid <- unique(mdata$smp_id)
    # list_smpid <- list_smpid[5] #[D22]
    # list_smpid <- list_smpid[8] #[D42]
  } else if (sr_mode=="mrg") {
    list_smpid <- c("mrg7") #concatenated_nBatch3
  }
  print(sr_mode) #checkpoint
  print(paste0("--- [list_smpid]:"))
  print(list_smpid)
  
  for (smpid in list_smpid) {
    ## input ##
    ## load the RDS: [sobj] - after UMAP (or the [sobj] created at the initial step by each sample) ####
    if (smpid=="mrg" smpid=="mrg7") {
      fn_rds_sobj <- paste0(dir_res_sr,"0_rds__GCAR1_PT01_GEX/","sobj_flt_umap_GCAR1_PT01_GEX-subsetVDJ_mrg7_raw_vdj_res2.4.rds") #GEX data w/ VDJ+ cells only
      
    } else if (smpid!="mrg" & smpid!="mrg7") {
      fn_rds_sobj <- paste0(dir_res_sr,"0_rds__",prj,"/","sobj_flt_",prj,"_",smpid,"_",sr_dtype,".rds") #flt
    }
    
    if (file.exists(fn_rds_sobj)) {
      sobj <- readRDS(fn_rds_sobj)
      print(paste0("[DONE] load the RDS: ",fn_rds_sobj))
    } else {
      print(paste0("[CHECK] RDS file does not exist: ",fn_rds_sobj))
    }
    # ---
    
    ## output ##
    print(sr_dtype) #check whether [sr_dtype] variable includes "vdj"
    ## set the H5AD "output" file name
    if (smpid=="mrg" | smpid=="mrg7") {
      fn_h5ad <- paste0(prjid,"_",smpid,"_afterUMAP","_",sr_dtype,"_res",res_clst,".h5ad")
      
    } else {
      # "sobj_GCAR1_PT01_GEX_Product_raw.rds"
      fn_h5ad <- paste0(prjid,"_",smpid,"_fltSobj","_",sr_dtype,".h5ad")
    }
    # path_h5ad <- paste0(dir_res_sr,"0_h5ad__",prj,"/")
    path_h5ad <- paste0(dir_res_sr,"0_h5ad__GCAR1_PT01_GEX/")
    print(paste0("--- path_h5ad: < ",path_h5ad," > ---"))
    fn_h5ad_full <- paste0(path_h5ad,fn_h5ad)
    print(fn_h5ad_full)
    # ---
    
    ## run [sceasy::convertFormat] ####
    if (!file.exists(fn_h5ad_full)) {
      print(paste0("[BEGIN] run 'convertFormat': ", Sys.time()))
      ## convert: Seurat to AnnData ####
      # sceasy::convertFormat(seurat_object, from="seurat", to="anndata", outFile='filename.h5ad') #template
      
      sceasy::convertFormat(sobj, 
                            from="seurat", to="anndata", 
                            outFile=fn_h5ad_full
      )
      
      print(paste0("[DONE] run 'convertFormat': ", Sys.time()))
      
    } else {
      print(paste0("[SKIP] h5ad file exists: ",fn_h5ad_full))
    } #closer: if (!file.exists(fn_h5ad_full)) {
    
    ## check whether the <h5ad> file is successfully created
    if (file.exists(fn_h5ad_full)) {
      print(paste0("[DONE] convert SeuratObject to h5ad: ", fn_h5ad_full))
    } else {
      print(paste0("[CHECK] h5ad file does not exist: ", fn_h5ad_full))
    }
  } #closer: for (smpid in list_smpid) {
} #closer: if (isTRUE(sceasy_cnvt_to_h5ad)) {
# ===