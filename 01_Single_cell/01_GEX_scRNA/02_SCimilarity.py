#!/usr/bin/python3
import sys
import os

################################
#### SCimilarity            ####
#### - cell type annotation ####
################################
# version: 2.45.1 (2024-12-15)
# by hyojinsong

from datetime import datetime
now=datetime.now()
today=now.strftime("%Y-%m-%d")
print(today)

#### running switches ####
## Seurat input data type used to create [sobj] ##
sr_dtype="raw" #"raw" #"filtered"
## Seurat UMAP clustering resolution ##
res_clst=2.4 #0.8(default) #2.4(selected)
res_clst_str=str(res_clst)
cname_clst="seurat_clusters"

TCELL_ONLY=1
if (TCELL_ONLY==1):
    sr_dtype=sr_dtype+"_vdj" #raw | filtered
else:
    sr_dtype=sr_dtype #raw | filtered

## doublet cell prediction ##
PRED_DUB=1 #use [Scrublet] h5ad #[master-running-switch]
RUN_SCR=1 #run [Scrublet]
## Scrublet paramters
rate_exp_doub=0.05 #0.05(default)
rate_exp_doub_str=str(rate_exp_doub)
SAVE_H5AD_SCRUB=0
SAVE_CSV_SCRUB=1

## cell type annotation ##
ANNO_CTYPE=1 #[master-running-switch]

COMPUTE_EMBED=1
MK_PLT_1=0 #1
MK_PLT_1a=1 #1
MK_PLT_1b=1 #UMAP(bySeuratCellCluster)
MK_PLT_1c=0 #1 #UMAP(markerGeneExpr)
MARKER_NATURECANCER=0 #plot 1c
MARKER_FLOWPANEL=0 #plot 1c
MARKER_TCELL=1 #plot 1c
MK_PLT_1d=1 #GCAR-expressing cells

ANNO_UNCONSTRAIN=1
EXT_CTYPE=0
MK_PLT_2=0
MK_PLT_2a=0
MK_PLT_2b=0
MK_PLT_2c=0
SAVE_H5AD_UNCONSTRAIN=0
SAVE_CSV_UNCONSTRAIN=1

ANNO_CONSTRAIN=1
MK_PLT_3=0
MK_PLT_3a=0
MK_PLT_3b=0
MK_PLT_3c=0
SAVE_H5AD_CONSTRAIN=0
SAVE_CSV_CONSTRAIN=1

ANNO_QC=0
MK_PLT_4=1 #note: [min_dist] required
MK_PLT_4a=0 #note: [min_dist] not included 
MK_PLT_4b=1 #note: [min_dist] requires the "constrained" approach
# ---

## marker gene expression plot ##
VIS_EXPR=0 #1 #[master-running-switch]

## dot plot ## - selected
# EXPR_PLT="GCARpos" #"all" "VDJ" "GCARpos" "GCARneg"
EXPR_MK_DOT=1
EXPR_MK_PLT_1=1 #1
EXPR_MK_PLT_1a=1 #1 #w/Dendrograms
EXPR_MK_PLT_1b=1 #transposed(swap_axes)
EXPR_MK_PLT_1c=1 #w/o categories
EXPR_MK_PLT_1d=1 #needs troubleshooting
EXPR_MK_PLT_1e=1 #byCluster
EXPR_MK_PLT_1f=1 #byCluster (w/Dendro)w

## scoreUMAP plots ##
EXPR_MK_SCORE=1 #1
SCR_PLT_1=1
SCR_PLT_1a=1

## heatmap ##
EXPR_MK_HMP_SC=1 #1 #[master-running-switch]
# ---

## cluster-level annotation ##
VIS_CLST=0 #[master-running-switch]

VIS_CLST_CCLASS=0 #plot by [cell class]
VIS_CLST_TPOINT=0 #plot by [timepoint]
VIS_CLST_CTYPE=1 #plot by [cell type]

CLST_MK_PLT_1=0
CLST_MK_PLT_1a=1
CLST_MK_PLT_2=0 #subplots (1 x 4)
CLST_MK_PLT_3=0 #subplots (4 x 8)
# ---


#### varirables ####
PRJNAME=sys.argv[1]
SMP=sys.argv[2]
if (VIS_EXPR==1):
    EXPR_MARKER_SET=sys.argv[3] #"CAR-T"
    EXPR_PLT=sys.argv[4] #"all" "VDJ" "GCARpos" "GCARneg" "VDJ_GCARpos"
    print("--- PRJNAME: ["+PRJNAME+"] | SMP: ["+SMP+"] - EXPR_MARKER: ["+EXPR_MARKER_SET+"] | PLOT: ["+EXPR_PLT+"] ---")
else:
    print("--- PRJNAME: ["+PRJNAME+"] | SMP: ["+SMP+"] ---")
## set the number of samples
if (SMP=="mrg" or SMP=="mrg7"):
    nSmp=7
    nSmp_str=str(nSmp)


## STEP 0 | Required software and data ####
# If the model hasn’t been downloaded please uncomment and run the two command below
## using [curl]
# !curl -L -o /models/model_v1.1.tar.gz \
#   https://zenodo.org/records/10685499/files/model_v1.1.tar.gz?download=1
# !tar -xzvf /models/model_v1.1.tar.gz
## using [wget]
# mkdir -p models & wget -O models/model_v1.1.tar.gz https://zenodo.org/records/10685499/files/model_v1.1.tar.gz?download=1

# If the data hasn’t been downloaded please uncomment and run the two command below
## using [curl]
# !curl -L -o "/data/GSE136831_subsample.h5ad" \
#   https://zenodo.org/record/8242083/files/GSE136831_subsample.h5ad?download=1
## using [wget]
# mkdir -p data & wget -O data/GSE136831_subsample.h5ad https://zenodo.org/record/8242083/files/GSE136831_subsample.h5ad?download=1

#### modules ####
import scanpy as sc
import anndata as ad

if (RUN_SCR==1):
    import scrublet as scr

from matplotlib import pyplot as plt

import numpy as np
# seed=42
# np.random.seed(seed)
import pandas as pd

import seaborn as sns
## [seaborn] palettes ##
# https://seaborn.pydata.org/tutorial/color_palettes.html
## [xkcd_rgb] colour palette ##
# https://xkcd.com/color/rgb/

# import colorcet as cc
import glasbey

## save multiple pages into one PDF file - skip for now
# from matplotlib.backends.backend_pdf import PdfPages

import PyPDF2 #merge PDF files
from pdf2image import convert_from_path #convert PDF-to-PNG files
# ---

## project ##
PRJNAME="GCAR1_PT01_GEX"
if (PRJNAME=="GCAR1_PT01_GEX"):
    PRJ="GCAR1_PT01"

## directory ##
HOMEDIR="/Users/your_username/"
DATADIR=HOMEDIR+"result/Seurat/0_h5ad__"+PRJNAME+"/"
print("--- DATADIR: "+DATADIR+" ---")

RESDIR=HOMEDIR+"result/SCimilarity/"
if (PRJ=="GCAR1_PT01"):
    RESDIR_PRJ=RESDIR+"SCim_"+PRJ+"/"+SMP+"/"
else:
    print("[CHECK] 'PRJ' variable")
isExist=os.path.exists(RESDIR_PRJ)
if not isExist:
    os.makedirs(RESDIR_PRJ)
    print("[DONE] create a directory: ", RESDIR_PRJ, sep="")
else:
    print("[SKIP] directory exists: ", RESDIR_PRJ, sep="")

if (VIS_EXPR==1):
    RESDIR_PRJ_EXPR=RESDIR_PRJ+"vis_expr"+"_res"+res_clst_str+"/"+EXPR_MARKER_SET+"_"+EXPR_PLT+"/"
    isExist=os.path.exists(RESDIR_PRJ_EXPR)
    if not isExist:
        os.makedirs(RESDIR_PRJ_EXPR)
        print("[DONE] create a directory: ", RESDIR_PRJ_EXPR, sep="")
    else:
        print("[SKIP] directory exists: ", RESDIR_PRJ_EXPR, sep="")

    if (EXPR_MK_SCORE==1):
        RESDIR_PRJ_SCR=RESDIR_PRJ+"vis_expr"+"_res"+res_clst_str+"_score"+"/"+EXPR_MARKER_SET+"_"+EXPR_PLT+"/"
        isExist=os.path.exists(RESDIR_PRJ_SCR)
        if not isExist:
            os.makedirs(RESDIR_PRJ_SCR)
            print("[DONE] create a directory: ", RESDIR_PRJ_SCR, sep="")
        else:
            print("[SKIP] directory exists: ", RESDIR_PRJ_SCR, sep="")    
    
if (VIS_CLST==1):
    RESDIR_PRJ_CLST=RESDIR_PRJ+"vis_clst"+"_res"+res_clst_str+"/"
    isExist=os.path.exists(RESDIR_PRJ_CLST)
    if not isExist:
        os.makedirs(RESDIR_PRJ_CLST)
        print("[DONE] create a directory: ", RESDIR_PRJ_CLST, sep="")
    else:
        print("[SKIP] directory exists: ", RESDIR_PRJ_CLST, sep="")
# ---

## set plot dimension
plt.rcParams["figure.figsize"] = [5, 5] #figure size in inches

import warnings

warnings.filterwarnings("ignore")

## set dot size - for UMAP plots
if (SMP!="mrg"):
    size_dot=20
elif (SMP=="mrg" or SMP=="mrg7"):
    size_dot=10
print("- dot size: ["+str(size_dot)+"]")

# ---
## STEP 0 | Predict doublet cells - using Scrublet ####
if (PRED_DUB==1 and RUN_SCR==1):    
    if (PRJ=="GCAR1_PT01"):
        if (TCELL_ONLY!=1):
            PATH_DATA=DATADIR+PRJ+"_"+SMP+"_afterUMAP"+"_raw"+"_res"+res_clst_str+".h5ad"
        elif (TCELL_ONLY==1):
            PATH_DATA=DATADIR+PRJ+"_"+SMP+"_afterUMAP"+"_raw_vdj"+"_res"+res_clst_str+".h5ad"
        if os.path.exists(PATH_DATA):
            print("[START] load the h5ad file: "+PATH_DATA)
        else:
            print("[CHECK] h5ad file does not exist: "+PATH_DATA)
    adata=sc.read(PATH_DATA)
    print(adata.shape) #checkpoint
    # (16214, 25099) #Harvest
    # (23163, 25469) #D08
    print("[END] load the h5ad file: "+PATH_DATA)
    
    ## [Scrublet] running environment requirement ##
    ## - source: [https://pmc.ncbi.nlm.nih.gov/articles/PMC6625319/]
    # "we timed our implementation for different numbers of cells, using a computer with a 2.1 GHz processor and 48 GB of memory, 
    # though this amount of memory was only necessary for the largest datasets."
    
    ## Scrublet (within Scanpy's built-in function) ##
    ## - source: [https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pp.scrublet.html]
    
    # Predict doublets using Scrublet [Wolock et al., 2019] - [https://doi.org/10.1016/j.cels.2018.11.005]

    # Predict cell doublets using a nearest-neighbor classifier of observed transcriptomes and simulated doublets. 

    # Works best if the input is a raw (unnormalized) counts matrix from a single sample or a collection of similar samples from the same experiment. 
    # This function is a wrapper around functions that pre-process using Scanpy and directly call functions of Scrublet(). 

    # You may also undertake your own preprocessing, simulate doublets with scrublet_simulate_doublets(), and run the core scrublet function scrublet() with adata_sim set.
    # ---
    
    ## [Scanpy:Scrublet] parameters
    # [sim_doublet_ratio]: float (default: 2.0)
    # Number of doublets to simulate relative to the number of observed transcriptomes.

    # [expected_doublet_rate]: float (default: 0.05)
    # Where adata_sim not suplied, the estimated doublet rate for the experiment.

    # [stdev_doublet_rate]: float (default: 0.02)
    # Where adata_sim not suplied, uncertainty in the expected doublet rate.
    # ---
    
    ## set Scrublet paramters
    # rate_exp_doub=0.05 #0.05(default)
    # rate_exp_doub_str=str(rate_exp_doub)
    
    ## run Scrublet - using Scanpy's built-in function:
    adata_scr=sc.pp.scrublet(
        adata,
        adata_sim=None, 
        batch_key=None, 
        sim_doublet_ratio=2.0, #2.0(default)
        expected_doublet_rate=rate_exp_doub, #0.05(default)
        stdev_doublet_rate=0.02, 
        synthetic_doublet_umi_subsampling=1.0, 
        knn_dist_metric='euclidean', 
        normalize_variance=True, 
        log_transform=False, 
        mean_center=True, 
        n_prin_comps=30, 
        use_approx_neighbors=None, 
        get_doublet_neighbor_parents=False, 
        n_neighbors=None, 
        threshold=None, 
        verbose=True, 
        copy=True, #False(default) #True(return a copy of the input adata with Scrublet results added. Otherwise, Scrublet results are added in place.)
        random_state=0
        )
    
    ## returns: [sc.pp.scrublet()]
    # if copy=True it returns or else adds fields to adata. Those fields:

    # .obs['doublet_score']: 
    # Doublet scores for each observed transcriptome

    # .obs['predicted_doublet']
    # Boolean indicating predicted doublet status

    # .uns['scrublet']['doublet_scores_sim']
    # Doublet scores for each simulated doublet transcriptome

    # .uns['scrublet']['doublet_parents']
    # Pairs of .obs_names used to generate each simulated doublet transcriptome

    # .uns['scrublet']['parameters']
    # Dictionary of Scrublet parameters
    
    ## save the doublet-predicted anndata as h5ad file
    if (SAVE_H5AD_SCRUB==1):
        # pip install hdf5plugin
        # conda activate hdf5
        import hdf5plugin
        
        if (TCELL_ONLY==0):
            FN_H5AD_NEW=DATADIR+PRJ+"_"+SMP+"_afterUMAP"+"_raw_"+"_res"+res_clst_str+"_Scrub_expRate"+rate_exp_doub_str+".h5ad"
        elif (TCELL_ONLY==1):
            FN_H5AD_NEW=DATADIR+PRJ+"_"+SMP+"_afterUMAP"+"_raw_vdj_"+"_res"+res_clst_str+"_Scrub_expRate"+rate_exp_doub_str+".h5ad"
        
        if not (os.path.exists(FN_H5AD_NEW)):
            adata_scr.write_h5ad(
                filename=FN_H5AD_NEW
            )
            ## checkpoint ##
            if (os.path.exists(FN_H5AD_NEW)):
                print("[DONE] save the doublet-predicted anndata as h5ad: "+FN_H5AD_NEW)
            else:
                print("[CHECK] doublet prediction-labelled h5ad file did not create successfully: "+FN_H5AD_NEW)
    
    ## save the specific column of the unconstrained prediction-labelled anndata as CSV file
    if (SAVE_CSV_SCRUB==1):
        ## usage ##
        # dirname = Name of the directory to which to export.
        # skip_data = Skip the data matrix X.
        
        if (TCELL_ONLY==0):
            OUTDIR_CSV=RESDIR_PRJ+"adata_scr-to-csv_expRate_"+rate_exp_doub_str+"_res"+res_clst_str+"/"
        elif (TCELL_ONLY==1):
            OUTDIR_CSV=RESDIR_PRJ+"adata_scr-to-csv_expRate_"+rate_exp_doub_str+"_res"+res_clst_str+"_VDJonly"+"/"
        
        isExist=os.path.exists(OUTDIR_CSV)
        if not isExist:
                os.makedirs(OUTDIR_CSV)
                print("[DONE] create a directory: ", OUTDIR_CSV, sep="")
        else:
            print("[SKIP] directory exists: ", OUTDIR_CSV, sep="")

        adata[adata.obs.predicted_doublet.isin()].write_csvs(dirname=OUTDIR_CSV, skip_data=True, sep=',')

        FN_CSV=OUTDIR_CSV+"obs.csv"
        ## checkpoint ##
        if (os.path.exists(FN_CSV)):
            print("[DONE] save the doublet-predicted anndata as CSV: "+FN_CSV)
        else:
            print("[CHECK] doublet prediction-labelled CSV file did not create successfully: "+FN_CSV)
# ===



# ---
## STEP 1 | Prepare for SCimilarity: Import and normalize data ####
if (ANNO_CTYPE==1):
    print("--- STEP 1 | prepare for SCimilarity ---")
    from scimilarity.utils import lognorm_counts, align_dataset
    if (ANNO_CTYPE==1):
        from scimilarity import CellAnnotation

    ## Import SCimilarity - Cell annotation object ####
    ## Instantiate the CellAnnotation object
    ## Set model_path to the location of the uncompressed model
    ## model (v1.1) ##
    model_path = HOMEDIR+"tool/scimilarity/models/model_v1.1"


if (ANNO_CTYPE==1):
    ## step 1a: load the input H5AD file ####
    if (PRJ=="GCAR1_PT01"):
        if (PRED_DUB!=1):
            if (TCELL_ONLY!=1):
                PATH_DATA=DATADIR+PRJ+"_"+SMP+"_afterUMAP"+"_raw"+"_res"+res_clst_str+".h5ad"
            elif (TCELL_ONLY==1):
                PATH_DATA=DATADIR+PRJ+"_"+SMP+"_afterUMAP"+"_raw_vdj"+"_res"+res_clst_str+".h5ad"
        elif (PRED_DUB==1):
            if (TCELL_ONLY!=1):
                PATH_DATA=DATADIR+PRJ+"_"+SMP+"_afterUMAP"+"_raw"+"_res"+res_clst_str+"_Scrub_expRate"+rate_exp_doub_str+".h5ad"
            elif (TCELL_ONLY==1):
                PATH_DATA=DATADIR+PRJ+"_"+SMP+"_afterUMAP"+"_raw_vdj"+"_res"+res_clst_str+"_Scrub_expRate"+rate_exp_doub_str+".h5ad"

            if os.path.exists(PATH_DATA):
                print("[START] load the h5ad file: "+PATH_DATA)
            else:
                print("[CHECK] h5ad file does not exist: "+PATH_DATA)
    adata=sc.read(PATH_DATA)
    print(adata.shape) #checkpoint
    # (16214, 25099) #Harvest
    # (23163, 25469) #D08
    print("[END] load the h5ad file: "+PATH_DATA)

    ## add the "pre-assigned (initial)" cell type annotations to start the prediction
    if (PRJ=="GCAR1_PT01"):
        if (SMP!="mrg"):
            ct_init=np.random.choice(["peripheral_blood_cell"], size=(adata.n_obs,))
            adata.obs["celltype_raw"] = pd.Categorical(ct_init) 
        elif (SMP=="mrg" or SMP=="mrg7"):
            adata.obs["celltype_raw"] = pd.Categorical(adata.obs["orig.ident"])
    adata.obs
    print(adata.T) #checkpoint    
# ---

    ## step 1b: SCimilarity pre-processing ####
    # Note, SCimilarity was trained with high data dropout to increase robustness to differences in gene lists.
    ca = CellAnnotation(model_path=model_path) #requires high memory
    
    ## Match feature space with SCimilarity models
    # "SCimilarity’s gene expression ordering is fixed. New data should be reorderd to match that, so that it is consistent with how the model was trained. 
    # Genes that are not present in the new data will be zero filled to comply to the expected structure. 
    # Genes that are not present in SCimilarity’s gene ordering will be filtered out."

    # Note, SCimilarity was trained with high data dropout to increase robustness to differences in gene lists.
    adata = align_dataset(adata, ca.gene_order)
    print("[DONE] run the 'align_dataset' function")
    
    ## Normalize data consistent with SCimilarity
    # "It is important to match Scimilarity’s normalization so that the data matches the lognorm tp10k procedure used during model training."
    adata = lognorm_counts(adata)
    print("[DONE] run the 'lognorm_counts' function")


## colour palettes ##
pal_smp={
    "Enriched_Apheresis":"#666666", #dark_grey
    "Enriched":"#666666", #dark_grey
    "GCAR1_Product":"#bf5b17", #brown
    "Harvest":"#bf5b17", #brown
    "D08":"#a6d854", #light_green
    "D15":"#e7298a", #pink
    "D22":"#e6ab02", #yellow
    "D28":"#1b9e77", #dark_green
    "D42":"#386cb0" #blue
    }

pal_smp_cat={
    "Enriched_Apheresis":"#666666", #dark_grey
    "Enriched":"#666666", #dark_grey
    "GCAR1_Product":"#bf5b17", #brown
    "Harvest":"#bf5b17", #brown
    "D08":"#a6d854", #light_green
    "D15":"#e7298a", #pink
    "D22":"#e6ab02", #yellow
    "D28":"#1b9e77", #dark_green
    "D42":"#386cb0", #blue
    "the rest":"lightgrey"
    }
pal_smp_cat_enr={"Enriched":"#666666", "the rest":"lightgrey"}
pal_smp_cat_hvt={"Harvest":"#bf5b17", "the rest":"lightgrey"}
pal_smp_cat_d08={"D08":"#a6d854", "the rest":"lightgrey"}
pal_smp_cat_d15={"D15":"#666666", "the rest":"lightgrey"}
pal_smp_cat_d22={"D22":"#e6ab02", "the rest":"lightgrey"}
pal_smp_cat_d28={"D28":"#1b9e77", "the rest":"lightgrey"}
pal_smp_cat_d42={"D42":"#386cb0", "the rest":"lightgrey"}


pal_dist=sns.color_palette("crest_r", as_cmap=True)
pal_expr_num=sns.color_palette("crest", as_cmap=True)
pal_expr_cat={"neg":"lightgrey", "pos":"#bd0026"}
pal_expr_cat_pos={"neg":"lightgrey", "pos":"#be0119"} #"scarlet=#be0119"
pal_expr_cat_neg={"neg":"#1e488f", "pos":"lightgrey"} #"cobalt"="#1e488f"


pal_expr_pos=sns.blend_palette(["lightgrey", sns.xkcd_rgb["scarlet"]], as_cmap=True)
pal_expr_neg=sns.blend_palette(["lightgrey", sns.xkcd_rgb["cobalt"]], as_cmap=True)

pal_expr_dot=sns.color_palette("YlOrRd", as_cmap=True)
pal_expr_vln=sns.color_palette("YlGnBu", as_cmap=True)
# ---


# ---
#### cell type annotation ####
## tutorial:
## [https://genentech.github.io/scimilarity/notebooks/cell_annotation_tutorial.html]

if (ANNO_CTYPE==1):    
    ## STEP 2 | Compute embeddings ####
    if (COMPUTE_EMBED==1):
        print("--- STEP 2 | compute embeddings ---")
        ## Using the already trained models, SCimilarity can embed your new dataset.
        adata.obsm["X_scimilarity"] = ca.get_embeddings(adata.X) #[adata.X] must include the "count" data

        ## Compute visualization of embeddings ####
        # Use UMAP to visualize SCimilarity embeddings
        sc.pp.neighbors(adata, use_rep="X_scimilarity")
        sc.tl.umap(adata)

        ## Visualize author annotations on the SCimilarity embedding ####       
        ## colour palettes - for Seurat cell clusters ##
        ## make "query" metadata columns for each Seurat cell cluster    
        ## - "As [sets] don’t contain any duplicate items then printing the length of the set will give us the total number of unique items."    
        list_clst=sorted(list(set(adata.obs[cname_clst])))
        # list_clst=sorted(list(set(adata.obs.seurat_clusters)))
        
        ## count the number of uniuqe Seurat cell clusters
        nClst=len(list_clst)
        nClst_str=str(nClst)

        colors_clst=glasbey.extend_palette("tab10", palette_size=nClst, colorblind_safe=True, cvd_severity=100) #default:lightness_bounds=(10,90)
        pal_clst=sns.set_palette(sns.color_palette(colors_clst))

        ## plot 1: UMAP plot - raw cell types ####
        if (MK_PLT_1==1):
            if (SMP!="mrg"):
                pal_smp=None
            elif (SMP=="mrg" or SMP=="mrg7"):
                pal_smp=pal_smp

            FN_PLOT_UMAP=RESDIR_PRJ+"1_UMAP_cellType_raw_"+PRJ+"_"+SMP+".png"
            isExist=os.path.exists(FN_PLOT_UMAP)
            if not isExist:
                mytitle="GCAR_PT01 (nSmp="+nSmp_str+")"
                umap=sc.pl.umap(
                    adata, 
                    color="celltype_raw", 
                    palette=pal_smp,
                    legend_fontsize=5,
                    s=size_dot,
                    title=mytitle
                    )
                plt.savefig(FN_PLOT_UMAP, orientation='landscape', bbox_inches='tight', pad_inches=0.1)
                if (os.path.exists(FN_PLOT_UMAP)):
                    print("[DONE] save as PDF: "+FN_PLOT_UMAP)
                

            FN_PLOT_UMAP=RESDIR_PRJ+"1_UMAP_cellType_raw_"+PRJ+"_"+SMP+"_byClst"+"_res"+res_clst_str+"_nClst"+nClst_str+".png" #".pdf"
            isExist=os.path.exists(FN_PLOT_UMAP)
            if not isExist:
                if (res_clst==0.8):
                    mytitle="Seurat cell clusters | res="+res_clst_str+" (default)"
                else:
                    mytitle="Seurat cell clusters | res="+res_clst_str
                umap=sc.pl.umap(
                    adata, 
                    color=cname_clst, 
                    palette=pal_clst,
                    legend_fontsize=5,
                    s=size_dot,
                    title=mytitle
                    )
                plt.savefig(FN_PLOT_UMAP, orientation='landscape', bbox_inches='tight', pad_inches=0.1)
                if (os.path.exists(FN_PLOT_UMAP)):
                    print("[DONE] save as PDF: "+FN_PLOT_UMAP)

        if (MK_PLT_1a==1):
            if (SMP=="mrg" or SMP=="mrg7"):
                
                ## make 7 UMAP subplots - w/ different colour palette and/or cmap parameters ##
                ncol=4
                nrow=2
                figsize=5
                wspace=0.5
                # Adapt figure size based on number of rows and columns and added space between them
                # (e.g. wspace between columns)
                fig,axs=plt.subplots(
                    nrows=nrow, ncols=ncol,
                    figsize=(ncol * figsize + (ncol - 1) * wspace * figsize, nrow * figsize)
                )
                plt.subplots_adjust(wspace=wspace)

                FN_PLOT_UMAP=RESDIR_PRJ+"1a_UMAP_cellType_raw_"+PRJ+"_"+SMP+"_bySmp_wAgg"+"_res"+res_clst_str+"_nClst"+nClst_str+".png"
                isExist=os.path.exists(FN_PLOT_UMAP)
                if not isExist:     
                    mytitle="GCAR_PT01 (nSmp="+nSmp_str+")"
                    umap0=sc.pl.umap(adata, color="sample_new", palette=pal_smp, show=False, legend_fontsize=5, s=size_dot, ax=axs[0,0], title=mytitle)
                    
                    umap1=sc.pl.umap(adata, color="sample_new", groups=["Enriched"], palette=pal_smp_cat, show=False, legend_fontsize=5, s=size_dot, ax=axs[0,1], title="Enriched") #"Product/Enriched"
                    umap2=sc.pl.umap(adata, color="sample_new", groups=["Harvest"], palette=pal_smp_cat, show=False, legend_fontsize=5, s=size_dot, ax=axs[0,2], title="Harvest") #"Tcells/Harvest"
                    umap3=sc.pl.umap(adata, color="sample_new", groups=["D08"], palette=pal_smp_cat, show=False, legend_fontsize=5, s=size_dot, ax=axs[0,3], title="D08")
                    
                    umap4=sc.pl.umap(adata, color="sample_new", groups=["D15"], palette=pal_smp_cat, show=False, legend_fontsize=5, s=size_dot, ax=axs[1,0], title="D15")
                    umap5=sc.pl.umap(adata, color="sample_new", groups=["D22"], palette=pal_smp_cat, show=False, legend_fontsize=5, s=size_dot, ax=axs[1,1], title="D22")
                    umap6=sc.pl.umap(adata, color="sample_new", groups=["D28"], palette=pal_smp_cat, show=False, legend_fontsize=5, s=size_dot, ax=axs[1,2], title="D28")
                    umap7=sc.pl.umap(adata, color="sample_new", groups=["D42"], palette=pal_smp_cat, show=False, legend_fontsize=5, s=size_dot, ax=axs[1,3], title="D42")
                    
                    ## edit UMAP plot legend labels
                    legend_texts = umap1.get_legend().get_texts()
                    for legend_text in legend_texts:
                        if legend_text.get_text() == "NA":
                            legend_text.set_text("the rest")
                    
                    legend_texts = umap2.get_legend().get_texts()
                    for legend_text in legend_texts:
                        if legend_text.get_text() == "NA":
                            legend_text.set_text("the rest")

                    legend_texts = umap3.get_legend().get_texts()
                    for legend_text in legend_texts:
                        if legend_text.get_text() == "NA":
                            legend_text.set_text("the rest")

                    legend_texts = umap4.get_legend().get_texts()
                    for legend_text in legend_texts:
                        if legend_text.get_text() == "NA":
                            legend_text.set_text("the rest")

                    legend_texts = umap5.get_legend().get_texts()
                    for legend_text in legend_texts:
                        if legend_text.get_text() == "NA":
                            legend_text.set_text("the rest")

                    legend_texts = umap6.get_legend().get_texts()
                    for legend_text in legend_texts:
                        if legend_text.get_text() == "NA":
                            legend_text.set_text("the rest")

                    legend_texts = umap7.get_legend().get_texts()
                    for legend_text in legend_texts:
                        if legend_text.get_text() == "NA":
                            legend_text.set_text("the rest")

                    plt.savefig(FN_PLOT_UMAP, orientation='landscape', bbox_inches='tight', pad_inches=0.1)
                    if (os.path.exists(FN_PLOT_UMAP)):
                        print("[DONE] save as PDF: "+FN_PLOT_UMAP)


        if (MK_PLT_1b==1):
            if (SMP=="mrg" or SMP=="mrg7"):
                # list_smpid=set(adata.obs["orig.ident"])          
                
                ## make 7 UMAP subplots - w/ different colour palette and/or cmap parameters ##
                ncol=4
                nrow=2
                figsize=5
                wspace=0.5
                # Adapt figure size based on number of rows and columns and added space between them
                # (e.g. wspace between columns)
                fig,axs=plt.subplots(
                    nrows=nrow, ncols=ncol,
                    figsize=(ncol * figsize + (ncol - 1) * wspace * figsize, nrow * figsize)
                )
                plt.subplots_adjust(wspace=wspace)

                FN_PLOT_UMAP=RESDIR_PRJ+"1b_UMAP_cellType_raw_"+PRJ+"_"+SMP+"_bySmp_wClst"+"_res"+res_clst_str+"_nClst"+nClst_str+".png"
                isExist=os.path.exists(FN_PLOT_UMAP)
                if not isExist:             
                    if (res_clst==0.8):
                        mytitle="Seurat cell clusters | res="+res_clst_str+" (default)"
                    else:
                        mytitle="Seurat cell clusters | res="+res_clst_str
                    umap0=sc.pl.umap(adata, color=cname_clst, palette=pal_clst, show=False, legend_fontsize=5, s=size_dot, ax=axs[0,0], title=mytitle)
                    
                    umap1=sc.pl.umap(adata, color="sample_new", groups=["Enriched"], palette=pal_smp_cat, show=False, legend_fontsize=5, s=size_dot, ax=axs[0,1], title="Enriched") #"Product/Enriched"
                    umap2=sc.pl.umap(adata, color="sample_new", groups=["Harvest"], palette=pal_smp_cat, show=False, legend_fontsize=5, s=size_dot, ax=axs[0,2], title="Harvest") #"Tcells/Harvest"
                    umap3=sc.pl.umap(adata, color="sample_new", groups=["D08"], palette=pal_smp_cat, show=False, legend_fontsize=5, s=size_dot, ax=axs[0,3], title="D08")
                    
                    umap4=sc.pl.umap(adata, color="sample_new", groups=["D15"], palette=pal_smp_cat, show=False, legend_fontsize=5, s=size_dot, ax=axs[1,0], title="D15")
                    umap5=sc.pl.umap(adata, color="sample_new", groups=["D22"], palette=pal_smp_cat, show=False, legend_fontsize=5, s=size_dot, ax=axs[1,1], title="D22")
                    umap6=sc.pl.umap(adata, color="sample_new", groups=["D28"], palette=pal_smp_cat, show=False, legend_fontsize=5, s=size_dot, ax=axs[1,2], title="D28")
                    umap7=sc.pl.umap(adata, color="sample_new", groups=["D42"], palette=pal_smp_cat, show=False, legend_fontsize=5, s=size_dot, ax=axs[1,3], title="D42")
                    
                    ## edit UMAP plot legend labels
                    legend_texts = umap1.get_legend().get_texts()
                    for legend_text in legend_texts:
                        if legend_text.get_text() == "NA":
                            legend_text.set_text("the rest")
                    
                    legend_texts = umap2.get_legend().get_texts()
                    for legend_text in legend_texts:
                        if legend_text.get_text() == "NA":
                            legend_text.set_text("the rest")

                    legend_texts = umap3.get_legend().get_texts()
                    for legend_text in legend_texts:
                        if legend_text.get_text() == "NA":
                            legend_text.set_text("the rest")

                    legend_texts = umap4.get_legend().get_texts()
                    for legend_text in legend_texts:
                        if legend_text.get_text() == "NA":
                            legend_text.set_text("the rest")

                    legend_texts = umap5.get_legend().get_texts()
                    for legend_text in legend_texts:
                        if legend_text.get_text() == "NA":
                            legend_text.set_text("the rest")

                    legend_texts = umap6.get_legend().get_texts()
                    for legend_text in legend_texts:
                        if legend_text.get_text() == "NA":
                            legend_text.set_text("the rest")

                    legend_texts = umap7.get_legend().get_texts()
                    for legend_text in legend_texts:
                        if legend_text.get_text() == "NA":
                            legend_text.set_text("the rest")

                    plt.savefig(FN_PLOT_UMAP, orientation='landscape', bbox_inches='tight', pad_inches=0.1)
                    if (os.path.exists(FN_PLOT_UMAP)):
                        print("[DONE] save as PDF: "+FN_PLOT_UMAP)
        
        if (MK_PLT_1c==1):
            ## make a "dictionary" of marker genes for each marker gene list
            
            ## CAR-T immunophenotyping markers ##
            if (MARKER_TCELL==1):
                ## - source: Figure 4 - Anderson ND et al. Nature Medicine 2023 [PMID:37407840]
                marker_gene_list_name="CAR-T"
                marker_genes={
                    ## Stem-like memory ##
                    "stem-like_memory": ["HNRNPLL","CCR7","SELL","LEF1","IL7R","TCF7"],
                    ## Activation, cytotoxicity, effector function ##
                    "activation_cytotoxicity_effector": ["CD28","CD27","GZMK","GZMA","PRF1","CCL4","GZMM","CCL5","NKG7","GZMB","GNLY","LYAR","CXCR3","CXCR4","TXNIP"],
                    ## Pre-exhaustion, exhaustion ##
                    "pre-exhaustion_exhaustion": ["TIGIT","ENTPD1","HAVCR2","CTLA4","LAG3","PDCD1","EOMES","TOX","PMCH","BATF","PRDM1","CXCR1"],
                    ## Resident ##
                    "resident": ["MCM5","TNFRSF9","ITGAE","CD69","NR4A2"]
                }
                
                list_mrkctype=[
                    "memory",
                    "effector",
                    "exhaustion",
                    "resident"
                ]
                
            ## for the "individual" data ##
            if (SMP!="mrg"):
                for mct in list_mrkctype:
                    print("--- marker cell type: ["+mct+"] ---")
                    
                    FN_PLOT_UMAP=RESDIR_PRJ+"1c_UMAP_cellType_raw_"+PRJ+"_"+SMP+"_ASPS+"+mct+".png" #changed to PNG(default) due to large file size issues
                    isExist=os.path.exists(FN_PLOT_UMAP)
                    if not isExist:
                        ## subset the dictionary by each "marker cell type":
                        
                        ## CAR-T immunophenotyping markers ##
                        if (MARKER_TCELL==1):
                            if mct=="memory":
                                mct_name="stem-like memory"
                                marker_genes_mct={key: value for key, value in marker_genes.items() if key in {
                                    "stem-like_memory"
                                }}
                            elif mct=="effector":
                                mct_name="Activation, cytotoxicity, effector function"
                                marker_genes_mct={key: value for key, value in marker_genes.items() if key in {
                                    "activation_cytotoxicity_effector"
                                }}
                            elif mct=="exhaustion":
                                mct_name="Pre-exhaustion, exhaustion"
                                marker_genes_mct={key: value for key, value in marker_genes.items() if key in {
                                    "pre-exhaustion_exhaustion"
                                }}
                            elif mct=="Resident":
                                mct_name="Pre-exhaustion, exhaustion"
                                marker_genes_mct={key: value for key, value in marker_genes.items() if key in {
                                    "resident"
                                }}

                            
                        print(marker_genes_mct) #checkpoint
                        
                        ## count the number of "cell types"
                        nMCT=len(marker_genes_mct)
                        
                        
                        ## Make Axes
                        ## Number of needed rows and columns (based on the row with the most columns)
                        # nrow = len(marker_genes_mct)
                        nrow = nMCT
                        ncol = max([len(vs) for vs in marker_genes_mct.values()])
                        
                        fig, axs = plt.subplots(nrow, ncol, figsize=(3 * ncol, 3 * nrow))
                        
                        ## colour palette ##
                        pal_expr=sns.blend_palette(["lightgrey", sns.xkcd_rgb["scarlet"]], as_cmap=True)
                        
                        # ---
                        ## Plot expression for every marker gene on the corresponding Axes object ##
                        for row_idx, (cell_type, markers) in enumerate(marker_genes_mct.items()):
                            col_idx = 0
                            for marker in markers:
                                if (nMCT==1):
                                    row_idx=0
                                ax = axs[row_idx, col_idx]
                                
                                sc.pl.umap(adata, color=marker, cmap=pal_expr, ax=ax, show=False, frameon=False, s=15)
                                
                                
                                ## Add cell type as row label - here we simply add it as ylabel of the first Axes object in the row
                                if col_idx == 0:
                                    ## We disabled axis drawing in UMAP to have plots without background and border
                                    ## so we need to re-enable axis to plot the ylabel
                                    ax.axis("on")
                                    ax.tick_params(
                                        top="off",
                                        bottom="off",
                                        left="off",
                                        right="off",
                                        labelleft="on",
                                        labelbottom="off"
                                    )
                                    ax.set_ylabel(cell_type + "\n", rotation=90, fontsize=14)
                                    ax.set(frame_on=False)
                                col_idx += 1
                            ## Remove unused column Axes in the current row
                            while col_idx < ncol:
                                axs[row_idx, col_idx].remove()
                                col_idx += 1
                        ## Alignment within the Figure
                        fig.tight_layout()
                        
                        ## plot title ##
                        fig.suptitle("GCAR_PT01: ["+SMP+"] - {"+mct_name+"}", fontsize=16, y=1.02) #y=0.98(default)
                        # ---
                        
                        ## save as PNG
                        plt.savefig(FN_PLOT_UMAP, orientation='landscape', bbox_inches='tight', pad_inches=0.1)
                        if (os.path.exists(FN_PLOT_UMAP)):
                            print("[DONE] save as PNG: "+FN_PLOT_UMAP)
                    else:
                        print("[SKIP] PNG file exists: "+FN_PLOT_UMAP)
            
            ## for the "merged" data ##
            elif (SMP=="mrg" or SMP=="mrg7"):

                # ---
                if (MARKER_NATURECANCER==1 or MARKER_FLOWPANEL==1 or MARKER_TCELL==1):
                    for cell_type in marker_genes.keys():
                        print("--- key: ["+cell_type+"] ---")
                        list_mg=marker_genes[cell_type]
                        mct=str(cell_type)
                        nKey=len(marker_genes.keys())

                        for mg in list_mg:
                            if (nKey==1):
                                FN_PLOT_UMAP=RESDIR_PRJ+"1c_UMAP_cellType_raw_"+PRJ+"_"+SMP+"_"+marker_gene_list_name+"__"+mg+".png" #changed to PNG(default) due to large file size issues
                            elif (nKey>1):
                                FN_PLOT_UMAP=RESDIR_PRJ+"1c_UMAP_cellType_raw_"+PRJ+"_"+SMP+"_"+marker_gene_list_name+"_"+mct+"__"+mg+".png" #changed to PNG(default) due to large file size issues
                            
                            isExist=os.path.exists(FN_PLOT_UMAP)
                            if not isExist:

                                ## plot every marker gene for "subset" of the merged data ##
                                ## make 4 UMAP subplots - w/ different colour palette and/or cmap parameters
                                ncol=7
                                nrow=3
                                figsize=5
                                wspace=0.1
                                
                                # Adapt figure size based on number of rows and columns and added space between them
                                # (e.g. wspace between columns)
                                fig,axs=plt.subplots(
                                    nrows=nrow, ncols=ncol,
                                    figsize=(ncol * figsize + (ncol - 1) * wspace * figsize, nrow * figsize)
                                )                                
                                plt.subplots_adjust(wspace=wspace)
                                
                                ## note: "ValueError: Groups and mask arguments are incompatible."
                                
                                ## row #1 ##
                                umap_1=sc.pl.umap(adata, mask_obs=(adata.obs["orig.ident"]=="Enriched_Apheresis"), color=mg, cmap=pal_expr_num, ax=axs[0,0], show=False, frameon=False, s=15, title="Enriched", vmin=0)
                                umap_2=sc.pl.umap(adata, mask_obs=(adata.obs["orig.ident"]=="GCAR1_Product"), color=mg, cmap=pal_expr_num, ax=axs[0,1], show=False, frameon=False, s=15, title="Harvest", vmin=0)
                                umap_3=sc.pl.umap(adata, mask_obs=(adata.obs["orig.ident"]=="D08"), color=mg, cmap=pal_expr_num, ax=axs[0,2], show=False, frameon=False, s=15, title="D08", vmin=0)
                                umap_4=sc.pl.umap(adata, mask_obs=(adata.obs["orig.ident"]=="D15"), color=mg, cmap=pal_expr_num, ax=axs[0,3], show=False, frameon=False, s=15, title="D15", vmin=0)
                                umap_5=sc.pl.umap(adata, mask_obs=(adata.obs["orig.ident"]=="D22"), color=mg, cmap=pal_expr_num, ax=axs[0,4], show=False, frameon=False, s=15, title="D22", vmin=0)
                                umap_6=sc.pl.umap(adata, mask_obs=(adata.obs["orig.ident"]=="D28"), color=mg, cmap=pal_expr_num, ax=axs[0,5], show=False, frameon=False, s=15, title="D28", vmin=0)
                                umap_7=sc.pl.umap(adata, mask_obs=(adata.obs["orig.ident"]=="D42"), color=mg, cmap=pal_expr_num, ax=axs[0,6], show=False, frameon=False, s=15, title="D42", vmin=0)
                                
                                
                                ## row #2 ##
                                umap_p1=sc.pl.umap(adata[adata.obs["orig.ident"]=="Enriched_Apheresis"], mask_obs=(adata[adata.obs["orig.ident"]=="Enriched_Apheresis"].obs.expr_GCAR=="pos"), color=mg, cmap=pal_expr_pos, ax=axs[1,0], show=False, frameon=False, s=15, title="", vmin=0)
                                umap_p2=sc.pl.umap(adata[adata.obs["orig.ident"]=="GCAR1_Product"], mask_obs=(adata[adata.obs["orig.ident"]=="GCAR1_Product"].obs.expr_GCAR=="pos"), color=mg, cmap=pal_expr_pos, ax=axs[1,1], show=False, frameon=False, s=15, title="", vmin=0)
                                umap_p3=sc.pl.umap(adata[adata.obs["orig.ident"]=="D08"], mask_obs=(adata[adata.obs["orig.ident"]=="D08"].obs.expr_GCAR=="pos"), color=mg, cmap=pal_expr_pos, ax=axs[1,2], show=False, frameon=False, s=15, title="", vmin=0)
                                umap_p4=sc.pl.umap(adata[adata.obs["orig.ident"]=="D15"], mask_obs=(adata[adata.obs["orig.ident"]=="D15"].obs.expr_GCAR=="pos"), color=mg, cmap=pal_expr_pos, ax=axs[1,3], show=False, frameon=False, s=15, title="", vmin=0)
                                umap_p5=sc.pl.umap(adata[adata.obs["orig.ident"]=="D22"], mask_obs=(adata[adata.obs["orig.ident"]=="D22"].obs.expr_GCAR=="pos"), color=mg, cmap=pal_expr_pos, ax=axs[1,4], show=False, frameon=False, s=15, title="", vmin=0)
                                umap_p6=sc.pl.umap(adata[adata.obs["orig.ident"]=="D28"], mask_obs=(adata[adata.obs["orig.ident"]=="D28"].obs.expr_GCAR=="pos"), color=mg, cmap=pal_expr_pos, ax=axs[1,5], show=False, frameon=False, s=15, title="", vmin=0)
                                umap_p7=sc.pl.umap(adata[adata.obs["orig.ident"]=="D42"], mask_obs=(adata[adata.obs["orig.ident"]=="D42"].obs.expr_GCAR=="pos"), color=mg, cmap=pal_expr_pos, ax=axs[1,6], show=False, frameon=False, s=15, title="", vmin=0)
                                
                                ## row #3 ##
                                umap_n1=sc.pl.umap(adata[adata.obs["orig.ident"]=="Enriched_Apheresis"], mask_obs=(adata[adata.obs["orig.ident"]=="Enriched_Apheresis"].obs.expr_GCAR=="neg"), color=mg, cmap=pal_expr_neg, ax=axs[2,0], show=False, frameon=False, s=15, title="", vmin=0)
                                umap_n2=sc.pl.umap(adata[adata.obs["orig.ident"]=="GCAR1_Product"], mask_obs=(adata[adata.obs["orig.ident"]=="GCAR1_Product"].obs.expr_GCAR=="neg"), color=mg, cmap=pal_expr_neg, ax=axs[2,1], show=False, frameon=False, s=15, title="", vmin=0)
                                umap_n3=sc.pl.umap(adata[adata.obs["orig.ident"]=="D08"], mask_obs=(adata[adata.obs["orig.ident"]=="D08"].obs.expr_GCAR=="neg"), color=mg, cmap=pal_expr_neg, ax=axs[2,2], show=False, frameon=False, s=15, title="", vmin=0)
                                umap_n4=sc.pl.umap(adata[adata.obs["orig.ident"]=="D15"], mask_obs=(adata[adata.obs["orig.ident"]=="D15"].obs.expr_GCAR=="neg"), color=mg, cmap=pal_expr_neg, ax=axs[2,3], show=False, frameon=False, s=15, title="", vmin=0)
                                umap_n5=sc.pl.umap(adata[adata.obs["orig.ident"]=="D22"], mask_obs=(adata[adata.obs["orig.ident"]=="D22"].obs.expr_GCAR=="neg"), color=mg, cmap=pal_expr_neg, ax=axs[2,4], show=False, frameon=False, s=15, title="", vmin=0)
                                umap_n6=sc.pl.umap(adata[adata.obs["orig.ident"]=="D28"], mask_obs=(adata[adata.obs["orig.ident"]=="D28"].obs.expr_GCAR=="neg"), color=mg, cmap=pal_expr_neg, ax=axs[2,5], show=False, frameon=False, s=15, title="", vmin=0)
                                umap_n7=sc.pl.umap(adata[adata.obs["orig.ident"]=="D42"], mask_obs=(adata[adata.obs["orig.ident"]=="D42"].obs.expr_GCAR=="neg"), color=mg, cmap=pal_expr_neg, ax=axs[2,6], show=False, frameon=False, s=15, title="", vmin=0)
                                # ---
                                
                                ## Add cell type as row label - here we simply add it as ylabel of the first Axes object in the row
                                ## row #1 ##
                                ax = axs[0,0]
                                ## We disabled axis drawing in UMAP to have plots without background and border
                                ## so we need to re-enable axis to plot the ylabel
                                ax.axis("on")
                                ax.tick_params(
                                    top="off",
                                    bottom="off",
                                    left="off",
                                    right="off",
                                    labelleft="on",
                                    labelbottom="on" #"off"
                                )
                                ylab="[all]"+"\n"
                                ax.set_ylabel(ylab, rotation=90, fontsize=14)
                                ax.set_xlabel("")
                                ax.set(frame_on=False)
                                
                                ## row #2 ##
                                ax = axs[1,0]
                                ## We disabled axis drawing in UMAP to have plots without background and border
                                ## so we need to re-enable axis to plot the ylabel
                                ax.axis("on")
                                ax.tick_params(
                                    top="off",
                                    bottom="off",
                                    left="off",
                                    right="off",
                                    labelleft="on",
                                    labelbottom="on" #"off"
                                )
                                # ax.set_ylabel("GCAR-pos" + "\n", rotation=90, fontsize=14)
                                ylab="[GCAR-pos]"+"\n"+mg+"-pos"+"\n"
                                ax.set_ylabel(ylab, rotation=90, fontsize=14)
                                ax.set_xlabel("")
                                ax.set(frame_on=False)
                                
                                ## row #3 ##
                                ax = axs[2,0]
                                ## We disabled axis drawing in UMAP to have plots without background and border
                                ## so we need to re-enable axis to plot the ylabel
                                ax.axis("on")
                                ax.tick_params(
                                    top="off",
                                    bottom="off",
                                    left="off",
                                    right="off",
                                    labelleft="on",
                                    labelbottom="off"
                                )
                                # ax.set_ylabel("GCAR-neg" + "\n", rotation=90, fontsize=14)
                                ylab="[GCAR-neg]"+"\n"+mg+"-pos"+"\n"
                                ax.set_ylabel(ylab, rotation=90, fontsize=14)                                
                                # ax.set_xlabel("") 
                                ax.set(frame_on=False)
                                
                                ## plot title ##
                                fig.suptitle("GCAR_PT01/"+SMP+": ["+cell_type+"] - {"+mg+"}", fontsize=16, y=0.93) #y=0.98(default)

                                ## save as PNG
                                plt.savefig(FN_PLOT_UMAP, orientation='landscape', bbox_inches='tight', pad_inches=0.1)
                                print("[DONE] save as PNG: "+FN_PLOT_UMAP)
                                # ---

                            else:
                                print("[SKIP] PNG file exists: "+FN_PLOT_UMAP)

        if (MK_PLT_1d==1):
            ## using a single GOI ##
            ## set the GOI name
            GOI="GCAR"
            OBS_NAME_GOI="expr_"+GOI
            
            if (SMP!="mrg" and SMP!="mrg7"):
                FN_PLOT_UMAP=RESDIR_PRJ+"1d_UMAP_exprGCAR_"+PRJ+"_"+SMP+"_w"+GOI+".png"
                isExist=os.path.exists(FN_PLOT_UMAP)
                if not isExist:
                ## UMAP plot ##
                    if (SMP!="Enriched_Apheresis" and SMP!="GCAR1_Product"):
                        mytitle="["+SMP+"] GCAR-expressing cells"
                    elif (SMP=="Enriched_Apheresis"):
                        mytitle="["+"Enriched"+"] GCAR-expressing cells"
                    elif (SMP=="GCAR1_Product"):
                        mytitle="["+"Harvest"+"] GCAR-expressing cells"
                    umap=sc.pl.umap(adata, 
                                    color=OBS_NAME_GOI, groups=["pos"], 
                                    palette=pal_expr_cat, 
                                    legend_fontsize=5, 
                                    s=size_dot, 
                                    show=False,
                                    title=mytitle
                                    )
                    ## We can change the 'NA' in the legend that represents all cells outside of the specified groups
                    legend_texts = umap.get_legend().get_texts()
                    ## Find legend object whose text is "NA" and change it
                    for legend_text in legend_texts:
                        if legend_text.get_text() == "NA":
                            legend_text.set_text("neg")
                            
                plt.savefig(FN_PLOT_UMAP, orientation='landscape', bbox_inches='tight', pad_inches=0.1)
                if (os.path.exists(FN_PLOT_UMAP)):
                    print("[DONE] save as PDF: "+FN_PLOT_UMAP)
                
            elif (SMP=="mrg" or SMP=="mrg7"):
                ## merged version (single UMAP plot) ##
                FN_PLOT_UMAP=RESDIR_PRJ+"1d_UMAP_exprGCAR_"+PRJ+"_"+SMP+"_w"+GOI+".png"
                isExist=os.path.exists(FN_PLOT_UMAP)
                if not isExist: 
                    ## make 1 UMAP plot ##
                    figsize=5

                    mytitle="[mrg] GCAR expression (nSmp=7)"
                    umap0=sc.pl.umap(adata, color=OBS_NAME_GOI, groups=["pos"], palette=pal_expr_cat, show=False, legend_fontsize=5, s=size_dot, title=mytitle)
                    
                    ## edit UMAP plot legend labels
                    legend_texts = umap0.get_legend().get_texts()
                    for legend_text in legend_texts:
                        if legend_text.get_text() == "NA":
                            legend_text.set_text("neg")
                    
                    plt.savefig(FN_PLOT_UMAP, orientation='landscape', bbox_inches='tight', pad_inches=0.1)
                    if (os.path.exists(FN_PLOT_UMAP)):
                        print("[DONE] save as PDF: "+FN_PLOT_UMAP)

                ## horizontal version (4 x 2) ##
                FN_PLOT_UMAP=RESDIR_PRJ+"1d_UMAP_exprGCAR_"+PRJ+"_"+SMP+"_w"+GOI+"__hori"+".png"
                isExist=os.path.exists(FN_PLOT_UMAP)
                if not isExist: 
                    ## make 8 UMAP subplots - w/ different colour palette and/or cmap parameters ##
                    ncol=4
                    nrow=2
                    figsize=5
                    wspace=0.5
                    # Adapt figure size based on number of rows and columns and added space between them
                    # (e.g. wspace between columns)
                    fig,axs=plt.subplots(
                        nrows=nrow, ncols=ncol,
                        figsize=(ncol * figsize + (ncol - 1) * wspace * figsize, nrow * figsize)
                    )
                    plt.subplots_adjust(wspace=wspace)

                    # mytitle="["+SMP+"] nSmp=7"
                    # umap0=sc.pl.umap(adata, color="orig.ident", palette=pal_smp, show=False, legend_fontsize=5, s=size_dot, ax=axs[0,0], title=mytitle)                    
                    mytitle="[mrg] GCAR expression (nSmp=7)"
                    umap0=sc.pl.umap(adata, color=OBS_NAME_GOI, groups=["pos"], palette=pal_expr_cat, show=False, legend_fontsize=5, s=size_dot, ax=axs[0,0], title=mytitle)

                    umap1=sc.pl.umap(adata[adata.obs["orig.ident"]=="Enriched_Apheresis"], color=OBS_NAME_GOI, groups=["pos"], palette=pal_expr_cat, show=False, legend_fontsize=5, s=size_dot, ax=axs[0,1], title="Enriched")
                    umap2=sc.pl.umap(adata[adata.obs["orig.ident"]=="GCAR1_Product"], color=OBS_NAME_GOI, groups=["pos"], palette=pal_expr_cat, show=False, legend_fontsize=5, s=size_dot, ax=axs[0,2], title="Harvest")
                    umap3=sc.pl.umap(adata[adata.obs["orig.ident"]=="D08"], color=OBS_NAME_GOI, groups=["pos"], palette=pal_expr_cat, show=False, legend_fontsize=5, s=size_dot, ax=axs[0,3], title="D08")
                    
                    umap4=sc.pl.umap(adata[adata.obs["orig.ident"]=="D15"], color=OBS_NAME_GOI, groups=["pos"], palette=pal_expr_cat, show=False, legend_fontsize=5, s=size_dot, ax=axs[1,0], title="D15")
                    umap5=sc.pl.umap(adata[adata.obs["orig.ident"]=="D22"], color=OBS_NAME_GOI, groups=["pos"], palette=pal_expr_cat, show=False, legend_fontsize=5, s=size_dot, ax=axs[1,1], title="D22")
                    umap6=sc.pl.umap(adata[adata.obs["orig.ident"]=="D28"], color=OBS_NAME_GOI, groups=["pos"], palette=pal_expr_cat, show=False, legend_fontsize=5, s=size_dot, ax=axs[1,2], title="D28")
                    umap7=sc.pl.umap(adata[adata.obs["orig.ident"]=="D42"], color=OBS_NAME_GOI, groups=["pos"], palette=pal_expr_cat, show=False, legend_fontsize=5, s=size_dot, ax=axs[1,3], title="D42")
                    
                    ## edit UMAP plot legend labels
                    legend_texts = umap0.get_legend().get_texts()
                    for legend_text in legend_texts:
                        if legend_text.get_text() == "NA":
                            legend_text.set_text("neg")
                    
                    legend_texts = umap1.get_legend().get_texts()
                    for legend_text in legend_texts:
                        if legend_text.get_text() == "NA":
                            legend_text.set_text("neg")
                    
                    legend_texts = umap2.get_legend().get_texts()
                    for legend_text in legend_texts:
                        if legend_text.get_text() == "NA":
                            legend_text.set_text("neg")

                    legend_texts = umap3.get_legend().get_texts()
                    for legend_text in legend_texts:
                        if legend_text.get_text() == "NA":
                            legend_text.set_text("neg")

                    legend_texts = umap4.get_legend().get_texts()
                    for legend_text in legend_texts:
                        if legend_text.get_text() == "NA":
                            legend_text.set_text("neg")

                    legend_texts = umap5.get_legend().get_texts()
                    for legend_text in legend_texts:
                        if legend_text.get_text() == "NA":
                            legend_text.set_text("neg")

                    legend_texts = umap6.get_legend().get_texts()
                    for legend_text in legend_texts:
                        if legend_text.get_text() == "NA":
                            legend_text.set_text("neg")

                    legend_texts = umap7.get_legend().get_texts()
                    for legend_text in legend_texts:
                        if legend_text.get_text() == "NA":
                            legend_text.set_text("neg")
                
                    plt.savefig(FN_PLOT_UMAP, orientation='landscape', bbox_inches='tight', pad_inches=0.1)
                    if (os.path.exists(FN_PLOT_UMAP)):
                        print("[DONE] save as PDF: "+FN_PLOT_UMAP)
                
                ## vertical version (1 x 8) ##
                FN_PLOT_UMAP=RESDIR_PRJ+"1d_UMAP_exprGCAR_"+PRJ+"_"+SMP+"_w"+GOI+"__vert"+".png"
                isExist=os.path.exists(FN_PLOT_UMAP)
                if not isExist: 
                    ## make 8 UMAP subplots - w/ different colour palette and/or cmap parameters ##
                    ncol=1
                    nrow=8
                    figsize=5
                    wspace=0.5
                    # Adapt figure size based on number of rows and columns and added space between them
                    # (e.g. wspace between columns)
                    fig,axs=plt.subplots(
                        nrows=nrow, ncols=ncol,
                        figsize=(ncol * figsize + (ncol - 1) * wspace * figsize, nrow * figsize)
                    )
                    plt.subplots_adjust(wspace=wspace)

                    # mytitle="["+SMP+"] nSmp=7"
                    # umap0=sc.pl.umap(adata, color="orig.ident", palette=pal_smp, show=False, legend_fontsize=5, s=size_dot, ax=axs[0], title=mytitle)
                    mytitle="[mrg] GCAR expression (nSmp=7)"
                    umap0=sc.pl.umap(adata, color=OBS_NAME_GOI, groups=["pos"], palette=pal_expr_cat, show=False, legend_fontsize=5, s=size_dot, ax=axs[0], title=mytitle)
                    
                    umap1=sc.pl.umap(adata[adata.obs["orig.ident"]=="Enriched_Apheresis"], color=OBS_NAME_GOI, groups=["pos"], palette=pal_expr_cat, show=False, legend_fontsize=5, s=size_dot, ax=axs[1], title="Enriched")
                    umap2=sc.pl.umap(adata[adata.obs["orig.ident"]=="GCAR1_Product"], color=OBS_NAME_GOI, groups=["pos"], palette=pal_expr_cat, show=False, legend_fontsize=5, s=size_dot, ax=axs[2], title="Harvest")
                    umap3=sc.pl.umap(adata[adata.obs["orig.ident"]=="D08"], color=OBS_NAME_GOI, groups=["pos"], palette=pal_expr_cat, show=False, legend_fontsize=5, s=size_dot, ax=axs[3], title="D08")
                    
                    umap4=sc.pl.umap(adata[adata.obs["orig.ident"]=="D15"], color=OBS_NAME_GOI, groups=["pos"], palette=pal_expr_cat, show=False, legend_fontsize=5, s=size_dot, ax=axs[4], title="D15")
                    umap5=sc.pl.umap(adata[adata.obs["orig.ident"]=="D22"], color=OBS_NAME_GOI, groups=["pos"], palette=pal_expr_cat, show=False, legend_fontsize=5, s=size_dot, ax=axs[5], title="D22")
                    umap6=sc.pl.umap(adata[adata.obs["orig.ident"]=="D28"], color=OBS_NAME_GOI, groups=["pos"], palette=pal_expr_cat, show=False, legend_fontsize=5, s=size_dot, ax=axs[6], title="D28")
                    umap7=sc.pl.umap(adata[adata.obs["orig.ident"]=="D42"], color=OBS_NAME_GOI, groups=["pos"], palette=pal_expr_cat, show=False, legend_fontsize=5, s=size_dot, ax=axs[7], title="D42")
                    
                    ## edit UMAP plot legend labels
                    legend_texts = umap0.get_legend().get_texts()
                    for legend_text in legend_texts:
                        if legend_text.get_text() == "NA":
                            legend_text.set_text("neg")
                    
                    legend_texts = umap1.get_legend().get_texts()
                    for legend_text in legend_texts:
                        if legend_text.get_text() == "NA":
                            legend_text.set_text("neg")
                    
                    legend_texts = umap2.get_legend().get_texts()
                    for legend_text in legend_texts:
                        if legend_text.get_text() == "NA":
                            legend_text.set_text("neg")

                    legend_texts = umap3.get_legend().get_texts()
                    for legend_text in legend_texts:
                        if legend_text.get_text() == "NA":
                            legend_text.set_text("neg")

                    legend_texts = umap4.get_legend().get_texts()
                    for legend_text in legend_texts:
                        if legend_text.get_text() == "NA":
                            legend_text.set_text("neg")

                    legend_texts = umap5.get_legend().get_texts()
                    for legend_text in legend_texts:
                        if legend_text.get_text() == "NA":
                            legend_text.set_text("neg")

                    legend_texts = umap6.get_legend().get_texts()
                    for legend_text in legend_texts:
                        if legend_text.get_text() == "NA":
                            legend_text.set_text("neg")

                    legend_texts = umap7.get_legend().get_texts()
                    for legend_text in legend_texts:
                        if legend_text.get_text() == "NA":
                            legend_text.set_text("neg")
                
                    plt.savefig(FN_PLOT_UMAP, orientation='landscape', bbox_inches='tight', pad_inches=0.1)
                    if (os.path.exists(FN_PLOT_UMAP)):
                        print("[DONE] save as PDF: "+FN_PLOT_UMAP)

    # ---
    ## STEP 3 | Cell type classification ####
    ## - source: [https://genentech.github.io/scimilarity/notebooks/cell_annotation_tutorial.html#3.-Cell-type-classification]
    
    # Two methods within the CellAnnotation class: 
    # 1. annotate_dataset - automatically computes embeddings. 
    # 2. get_predictions - more detailed control of annotation.

    # Description of inputs:
    # - X_scimilarity: embeddings from the model, which can be used to generate UMAPs in lieu of PCA and is generalized across datasets.
    
    # Description of outputs:
    # - predictions: cell type annotation predictions. 
    # - nn_idxs: indicies of cells in the SCimilarity reference. 
    # - nn_dists: the minimum distance within k=50 nearest neighbors. 
    # - nn_stats: a dataframe containing useful metrics such as: 
    #   - hits: the distribution of celltypes in k=50 nearest neighbors.
    
    ## step 3a: cell type classification - using Unconstrained annotation ####
    if (ANNO_UNCONSTRAIN==1):
        print("--- STEP 3a | cell type classification - using Unconstrained annotation ---")
        ## Cells can be classified as any type that is in the SCimilarity reference
        predictions, nn_idxs, nn_dists, nn_stats = ca.get_predictions_kNN(
            adata.obsm["X_scimilarity"]
        )
        adata.obs["predictions_unconstrained"] = predictions.values

        ## Since each cell is classified independently, there is higher classification noise, filtering out low count cells can reduce the noise in visualization.
        celltype_counts = adata.obs.predictions_unconstrained.value_counts()
        well_represented_celltypes = celltype_counts[celltype_counts > 20].index
        ## make a list of unique variables for the [well_represented_celltypes]
        list_well_represented_celltypes=set(well_represented_celltypes)
        nCType_unc=len(list_well_represented_celltypes)
        nCType_unc_str=str(nCType_unc)
        
        if (EXT_CTYPE==1):
            ## save as TXT: [well_represented_celltypes]
            OUTDIR_CSV=RESDIR_PRJ+"adata-to-csv_constrain"+"_res"+res_clst_str+"/"
            FN_TXT=OUTDIR_CSV+"list_nCellType"+nCType_unc_str+".txt"        
            df_well_represented_celltypes=pd.DataFrame(list_well_represented_celltypes, columns=['well_represented_celltypes'])
            df_well_represented_celltypes.write_csvs(dirname=OUTDIR_CSV, skip_data=True, sep="\t")
            ## checkpoint ##
            if (os.path.exists(FN_TXT)):
                print("[DONE] save the constrained prediction-labelled anndata as TXT: "+FN_TXT)
            else:
                print("[CHECK] constrained prediction-labelled TXT file did not create successfully: "+FN_TXT)
        
        
        ## plot 2: UMAP plot - "unconstrained" prediction-based cell types ####
        if (MK_PLT_2==1):
            FN_PLOT_UMAP=RESDIR_PRJ+"2_UMAP_cellType_pred-unconstrained_"+PRJ+"_"+SMP+".pdf"
            isExist=os.path.exists(FN_PLOT_UMAP)
            if not isExist:
                ## count the number of cell types included in the "unconstrained" prediction result
                list_ctype=list(set(well_represented_celltypes))
                nCType=len(list_ctype)
                
                ## make a customized colour palette
                colors=glasbey.extend_palette("Dark2", palette_size=nCType, colorblind_safe=True, cvd_severity=100) #default:lightness_bounds=(10,90)
                pal=sns.set_palette(sns.color_palette(colors))

                sc.pl.umap(
                    adata[adata.obs.predictions_unconstrained.isin(well_represented_celltypes)],
                    color="predictions_unconstrained",
                    # cmap=sns.color_palette(palette=pal, as_cmap=True),
                    palette=pal, 
                    legend_fontsize=5,
                    s=size_dot
                )
                # plt.savefig(FN_PLOT_UMAP)
                plt.savefig(FN_PLOT_UMAP, orientation='landscape', bbox_inches='tight', pad_inches=0.1)
                if (os.path.exists(FN_PLOT_UMAP)):
                    print("[DONE] save as PDF: "+FN_PLOT_UMAP)

        ## plot 2a: UMAP plot - "unconstrained" prediction-based cell types (for each cell type) ####
        if (MK_PLT_2a==1):
            ## make a list of cell types predicted from the "unconstrained" mode
            list_ctype=list(set(well_represented_celltypes))
            ## sort the list
            list_ctype.sort(reverse=False)
            
            ## count the number of cell types included in the "unconstrained" prediction result
            nCType=len(list_ctype)
            
            for ctype in list_ctype:
                ## define a customized function: [order]
                def order(item, list, n=0):
                    if list:
                        if list[0] == item:
                            return n
                        elif len(list) >= 2: # for python 2, use "else:"
                            return order(item, list[1:], n+1)
                
                ## add "zero" before the order number of each cell type if the number is 1-digit
                ordCType=order(item=ctype, list=list_ctype)+1
                if (ordCType <10):
                    ordCType_str="0"+str(ordCType)
                else:
                    ordCType_str=str(ordCType)

                print("---- selected cell type: ["+ctype+"]"+" ("+ordCType_str+" /"+str(nCType)+") ---")
                
                ## make a customized colour palette
                colors=glasbey.extend_palette("Dark2", palette_size=nCType, colorblind_safe=True, cvd_severity=100) #default:lightness_bounds=(10,90)
                pal=sns.set_palette(sns.color_palette(colors))

                col_ct=colors[ordCType-1]
                print(col_ct)
                
                ax = sc.pl.umap(
                    adata[adata.obs.predictions_unconstrained.isin(well_represented_celltypes)],
                    color=["predictions_unconstrained"], groups=[ctype], 
                    palette=pal, 
                    show=False, #Show the plot, do not return axis. #default: None
                    legend_fontsize=5,
                    s=size_dot,
                    na_color='lightgrey'
                    )
                
                ## make another version of "cell type" label - i.e., each space(" ") is replaced by an underscore("_")
                ctype_spl=ctype.split(" ")
                ctype_rpl="_".join(ctype_spl)
                    
                FN_PLOT_UMAP=RESDIR_PRJ+"2a_UMAP_cellType_pred-unconstrained_"+PRJ+"_"+SMP+"_"+str(nCType)+"-"+ordCType_str+"_"+ctype_rpl+".pdf"
                isExist=os.path.exists(FN_PLOT_UMAP)
                if not isExist:
                    # We can change the 'NA' in the legend that represents all cells outside of the specified groups
                    legend_texts = ax.get_legend().get_texts()
                    # Find legend object whose text is "NA" and change it
                    for legend_text in legend_texts:
                        if legend_text.get_text() == "NA":
                            legend_text.set_text("other cell types")
                            
                            plt.savefig(FN_PLOT_UMAP, orientation='landscape', bbox_inches='tight', pad_inches=0.1)
                            if (os.path.exists(FN_PLOT_UMAP)):
                                print("[DONE] save as PDF: "+FN_PLOT_UMAP)
                    
            ## Coloring cell subset ##
            ## Here we show how we can plot all cells as a background and then plot on top indivdual cell groups in color.
            ## We can color-in only specific cell groups when using categorical colors with the groups parameter.
            # ax = sc.pl.umap(adata, color=["bulk_labels"], groups=["Dendritic"], show=False)
            
            ## We can also plot continous values of an individual cell group using the obs_mask key word argument:
            # sc.pl.umap(adata, color="IGJ", mask_obs=(adata.obs.bulk_labels == "CD19+ B"), size=20)


        ## plot 2b: "unconstrained" prediction UMAP plot - w/ GOI expression (aggregated version) ####
        if (MK_PLT_2b==1):    
            ## using a single GOI ##
            ## set the GOI name
            GOI="GCAR"
            OBS_NAME_GOI="expr_"+GOI
            
            FN_PLOT_UMAP=RESDIR_PRJ+"2b_UMAP_cellType_pred-unconstrained_"+PRJ+"_"+SMP+"_w"+GOI+".pdf"
            isExist=os.path.exists(FN_PLOT_UMAP)
            if not isExist:
                ## make a list of variables for cell-colour annotation
                var_color=[
                    OBS_NAME_GOI,
                    "predictions_unconstrained"
                ]
                
                ## make a list of cell types predicted from the "unconstrained" mode
                list_ctype=list(set(well_represented_celltypes))
                ## sort the list
                list_ctype.sort(reverse=False)
                
                ## count the number of cell types included in the "constrained" prediction result
                nCType=len(list_ctype)

                ## make customized colour palettes
                colors=glasbey.extend_palette("Dark2", palette_size=nCType, colorblind_safe=True, cvd_severity=100) #default:lightness_bounds=(10,90)
                pal=sns.set_palette(sns.color_palette(colors))
                pal_expr={"neg":"lightgrey", "pos":"#bd0026"}

                ## make two UMAP plots as subplots
                ## option 1: if UMAP plots share the colour palette and/or cmap parameters ##
                # sc.pl.umap(
                #     # adata[adata.obs.predictions_unconstrained.isin(well_represented_celltypes)],
                #     adata,
                #     color=var_color, 
                #     cmap=sns.blend_palette(["lightgrey", sns.xkcd_rgb["scarlet"]], as_cmap=True),
                #     palette=pal,
                #     legend_fontsize=5,
                #     s=size_dot,
                #     wspace=0.5 #space btwn two UMAP plots #default:None
                # )
                
                ## option 2: if UMAP plots does NOT share the colour palette and/or cmap parameters ##
                ncol=2
                nrow=1
                figsize=5
                wspace=0.5
                # Adapt figure size based on number of rows and columns and added space between them
                # (e.g. wspace between columns)
                fig,axs=plt.subplots(
                    nrows=nrow, ncols=ncol,
                    figsize=(ncol * figsize + (ncol - 1) * wspace * figsize, nrow * figsize)
                )
                plt.subplots_adjust(wspace=wspace)

                ## UMAP plot (left side) ##
                umap1=sc.pl.umap(adata, color=OBS_NAME_GOI, groups=["pos"], palette=pal_expr_cat, legend_fontsize=5, s=size_dot, ax=axs[0], show=False)
                ## We can change the 'NA' in the legend that represents all cells outside of the specified groups
                legend_texts = umap1.get_legend().get_texts()
                ## Find legend object whose text is "NA" and change it
                for legend_text in legend_texts:
                    if legend_text.get_text() == "NA":
                        legend_text.set_text("neg")

                        ## UMAP plot (right side) ##
                        sc.pl.umap(adata[adata.obs.predictions_unconstrained.isin(well_represented_celltypes)], color="predictions_unconstrained", palette=pal, legend_fontsize=5, s=size_dot, ax=axs[1])
                        
                        plt.savefig(FN_PLOT_UMAP, orientation='landscape', bbox_inches='tight', pad_inches=0.1)
                        if (os.path.exists(FN_PLOT_UMAP)):
                            print("[DONE] save as PDF: "+FN_PLOT_UMAP)
        

        ## plot 2c: "unconstrained" prediction UMAP plot - w/ GOI expression (for each cell type) ####
        if (MK_PLT_2c==1):
            ## using a single GOI ##
            ## set the GOI name
            GOI="GCAR"
            OBS_NAME_GOI="expr_"+GOI
            
            ## make a list of cell types predicted from the "unconstrained" mode
            list_ctype=list(set(well_represented_celltypes))
            ## sort the list
            list_ctype.sort(reverse=False)     

            ## count the number of cell types included in the "unconstrained" prediction result
            nCType=len(list_ctype)
            
            for ctype in list_ctype:
                ## define a customized function: [order]
                def order(item, list, n=0):
                    if list:
                        if list[0] == item:
                            return n
                        elif len(list) >= 2: # for python 2, use "else:"
                            return order(item, list[1:], n+1)
                
                ## add "zero" before the order number of each cell type if the number is 1-digit
                ordCType=order(item=ctype, list=list_ctype)+1
                if (ordCType <10):
                    ordCType_str="0"+str(ordCType)
                else:
                    ordCType_str=str(ordCType)

                print("---- [expr_"+GOI+"] + selected cell type: ["+ctype+"]"+" ("+ordCType_str+" /"+str(nCType)+") ---")         

                ## make customized colour palettes
                colors=glasbey.extend_palette("Dark2", palette_size=nCType, colorblind_safe=True, cvd_severity=100) #default:lightness_bounds=(10,90)
                pal=sns.set_palette(sns.color_palette(colors))
                pal_expr={"neg":"lightgrey", "pos":"#bd0026"}          

                ## make 2 UMAP subplots - w/ different colour palette and/or cmap parameters ##
                ncol=2
                nrow=1
                figsize=5
                wspace=0.5
                # Adapt figure size based on number of rows and columns and added space between them
                # (e.g. wspace between columns)
                fig,axs=plt.subplots(
                    nrows=nrow, ncols=ncol,
                    figsize=(ncol * figsize + (ncol - 1) * wspace * figsize, nrow * figsize)
                )
                plt.subplots_adjust(wspace=wspace)         
                
                ## make another version of "cell type" label - i.e., each space(" ") is replaced by an underscore("_")
                ctype_spl=ctype.split(" ")
                ctype_rpl="_".join(ctype_spl)

                FN_PLOT_UMAP=RESDIR_PRJ+"2c_UMAP_cellType_pred-unconstrained_"+PRJ+"_"+SMP+"_w"+GOI+"_"+str(nCType)+"-"+ordCType_str+"_"+ctype_rpl+".pdf"
                isExist=os.path.exists(FN_PLOT_UMAP)
                if not isExist:                
                    ## UMAP plot #1 (left side) ##
                    umap1=sc.pl.umap(adata, 
                                    color=OBS_NAME_GOI, groups=["pos"], 
                                    palette=pal_expr_cat, 
                                    legend_fontsize=5, 
                                    s=size_dot, 
                                    ax=axs[0], 
                                    show=False
                                    )
                    ## We can change the 'NA' in the legend that represents all cells outside of the specified groups
                    legend_texts = umap1.get_legend().get_texts()
                    ## Find legend object whose text is "NA" and change it
                    for legend_text in legend_texts:
                        if legend_text.get_text() == "NA":
                            legend_text.set_text("neg")

                            ## UMAP plot #2 (right side) ##
                            umap2=sc.pl.umap(adata[adata.obs.predictions_unconstrained.isin(well_represented_celltypes)],
                                            color=["predictions_unconstrained"], groups=[ctype], 
                                            palette=pal, 
                                            show=False, #Show the plot, do not return axis. #default: None
                                            legend_fontsize=5,
                                            s=size_dot,
                                            na_color='lightgrey',
                                            ax=axs[1]
                                            )
                            # We can change the 'NA' in the legend that represents all cells outside of the specified groups
                            legend_texts = umap2.get_legend().get_texts()
                            # Find legend object whose text is "NA" and change it
                            for legend_text in legend_texts:
                                if legend_text.get_text() == "NA":
                                    legend_text.set_text("other cell types")

                            plt.savefig(FN_PLOT_UMAP, orientation='landscape', bbox_inches='tight', pad_inches=0.1)
                            ## checkpoint ##
                            if (os.path.exists(FN_PLOT_UMAP)):
                                print("[DONE] save as PDF: "+FN_PLOT_UMAP)



        ## using a GOI list ## - skip for now
        # list_goi=["ASPSCR1","TFE3","GPNMB","GCAR"]
        # nGOI=len(list_goi)
        
        # FN_PLOT_UMAP=RESDIR_PRJ+"2b_UMAP_cellType_pred-unconstrained_"+PRJ+"_"+SMP+"_GOI"+str(nGOI)+".pdf"
        # isExist=os.path.exists(FN_PLOT_UMAP)
        # if not isExist:
        #     sc.pl.umap(
        #         # adata[adata.obs.predictions_unconstrained.isin(well_represented_celltypes)],
        #         adata,
        #         color=list_goi,
        #         legend_fontsize=5
        #     )
        #     # plt.savefig(FN_PLOT_UMAP)
        #     plt.savefig(FN_PLOT_UMAP, orientation='landscape', bbox_inches='tight', pad_inches=0.1)
        #     if (os.path.exists(FN_PLOT_UMAP)):
        #         print("[DONE] save 'plot 2c' as PDF: "+FN_PLOT_UMAP)



        ## save the results ## 
        ## save the unconstrained prediction-labelled anndata as h5ad file
        if (SAVE_H5AD_UNCONSTRAIN==1):
            # pip install hdf5plugin
            # conda activate hdf5
            import hdf5plugin
            
            # PATH_DATA=DATADIR+PRJ+"_"+SMP+"_afterUMAP"+".h5ad"
            if (TCELL_ONLY==0):
                FN_H5AD_NEW=DATADIR+PRJ+"_"+SMP+"_afterUMAP"+"_SCim_unconstrain"+"_res"+res_clst_str+".h5ad"
            elif (TCELL_ONLY==1):
                FN_H5AD_NEW=DATADIR+PRJ+"_"+SMP+"_afterUMAP"+"_SCim_unconstrain"+"_VDJonly"+"_res"+res_clst_str+".h5ad"
            
            if not (os.path.exists(FN_H5AD_NEW)):
                adata.write_h5ad(
                    # anndata=adata,
                    filename=FN_H5AD_NEW
                    # compression=hdf5plugin.FILTERS["zstd"]
                )
                ## checkpoint ##
                if (os.path.exists(FN_H5AD_NEW)):
                    print("[DONE] save the unconstrained prediction-labelled anndata as h5ad: "+FN_H5AD_NEW)
                else:
                    print("[CHECK] unconstrained prediction-labelled h5ad file did not create successfully: "+FN_H5AD_NEW)
        
        ## save the specific column of the unconstrained prediction-labelled anndata as CSV file
        if (SAVE_CSV_UNCONSTRAIN==1):
            ## usage ##
            # dirname = Name of the directory to which to export.
            # skip_data = Skip the data matrix X.
            
            if (PRED_DUB!=1):
                OUTDIR_CSV=RESDIR_PRJ+"adata-to-csv_unconstrain_"+sr_dtype+"_res"+res_clst_str+"/"
            else:
                OUTDIR_CSV=RESDIR_PRJ+"adata-to-csv_unconstrain_"+sr_dtype+"_res"+res_clst_str+"_expRate"+rate_exp_doub_str+"/"            
            
            isExist=os.path.exists(OUTDIR_CSV)
            if not isExist:
                    os.makedirs(OUTDIR_CSV)
                    print("[DONE] create a directory: ", OUTDIR_CSV, sep="")
            else:
                print("[SKIP] directory exists: ", OUTDIR_CSV, sep="")

            adata[adata.obs.predictions_unconstrained.isin(well_represented_celltypes)].write_csvs(dirname=OUTDIR_CSV, skip_data=True, sep=',')

            FN_CSV=OUTDIR_CSV+"obs.csv"
            ## checkpoint ##
            if (os.path.exists(FN_CSV)):
                print("[DONE] save the unconstrained prediction-labelled anndata as CSV: "+FN_CSV)
            else:
                print("[CHECK] unconstrained prediction-labelled CSV file did not create successfully: "+FN_CSV)




    ## step 3b: cell type classification - using Constrained classification ####
    if (ANNO_CONSTRAIN==1):
        print("--- STEP 3b | cell type classification - using Constrained annotation ---")
        # By classifying against the full reference, we can get redundant cell types, such as activated CD8-positive, alpha-beta T cell and CD8-positive, alpha-beta T cell.

        # Alternatively, we can subset the reference to just the cell types we want to classify to. This also reduces noise in cell type annotation.

        # Note, subsetting can slow classification speeds as the kNN is optimized for the full reference.
        if (TCELL_ONLY==0):
            if (SMP!="mrg7"):
                # nCType=59
                target_celltypes = [
                    "alpha-beta T cell",
                    "B cell",
                    "CD14-low, CD16-positive monocyte",
                    "CD14-positive monocyte",
                    "CD14-positive, CD16-positive monocyte",
                    "CD16-negative, CD56-bright natural killer cell, human",
                    "CD16-positive, CD56-dim natural killer cell, human",
                    "CD4-positive helper T cell",
                    "CD4-positive, alpha-beta cytotoxic T cell",
                    "CD4-positive, alpha-beta memory T cell",
                    "CD4-positive, alpha-beta T cell",
                    "CD8-positive, alpha-beta cytotoxic T cell",
                    "CD8-positive, alpha-beta memory T cell",
                    "CD8-positive, alpha-beta T cell",
                    "central memory CD4-positive, alpha-beta T cell",
                    "central memory CD8-positive, alpha-beta T cell",
                    "class switched memory B cell",
                    "classical monocyte",
                    "common lymphoid progenitor",
                    "conventional dendritic cell",
                    "dendritic cell",
                    "dendritic cell, human",
                    "effector CD8-positive, alpha-beta T cell",
                    "effector memory CD4-positive, alpha-beta T cell",
                    "effector memory CD8-positive, alpha-beta T cell",
                    "effector memory CD8-positive, alpha-beta T cell, terminally differentiated",
                    "erythrocyte",
                    "erythroid progenitor cell",
                    "gamma-delta T cell",
                    "hematopoietic precursor cell",
                    "hematopoietic stem cell",
                    "IgA plasma cell", 
                    "IgG plasma cell",
                    "innate lymphoid cell",
                    "leukocyte",
                    "lymphocyte",
                    "macrophage",
                    "mast cell",       
                    "mature NK T cell",
                    "memory B cell",
                    "monocyte",
                    "mucosal invariant T cell",
                    "myeloid cell",
                    "myeloid dendritic cell",
                    "naive B cell",
                    "naive thymus-derived CD4-positive, alpha-beta T cell",
                    "naive thymus-derived CD8-positive, alpha-beta T cell",
                    "native cell",
                    "natural killer cell",
                    "neutrophil",
                    "non-classical monocyte",
                    "plasma cell",
                    "plasmablast",
                    "plasmacytoid dendritic cell",
                    "platelet",
                    "progenitor cell",
                    "regulatory T cell",
                    "T cell",
                    "T follicular helper cell"
                    ]
        
        elif (TCELL_ONLY==1):
            if (SMP!="mrg7"):
                # nCType=46
                target_celltypes = [
                    "B cell",
                    "CD14-low, CD16-positive monocyte",
                    "CD14-positive monocyte",
                    "CD14-positive, CD16-positive monocyte",
                    "CD16-negative, CD56-bright natural killer cell, human",
                    "CD16-positive, CD56-dim natural killer cell, human",
                    "CD4-positive helper T cell",
                    "CD4-positive, alpha-beta cytotoxic T cell",
                    "CD4-positive, alpha-beta memory T cell",
                    "CD4-positive, alpha-beta T cell",
                    "CD8-positive, alpha-beta cytotoxic T cell",
                    "CD8-positive, alpha-beta memory T cell",
                    "CD8-positive, alpha-beta T cell",
                    "central memory CD4-positive, alpha-beta T cell",
                    "central memory CD8-positive, alpha-beta T cell",
                    "classical monocyte",
                    "conventional dendritic cell",
                    "dendritic cell",
                    "dendritic cell, human",
                    "effector CD8-positive, alpha-beta T cell",
                    "effector memory CD4-positive, alpha-beta T cell",
                    "effector memory CD8-positive, alpha-beta T cell",
                    "effector memory CD8-positive, alpha-beta T cell, terminally differentiated",
                    "erythrocyte",
                    "gamma-delta T cell",
                    "hematopoietic stem cell",
                    "IgA plasma cell",
                    "leukocyte",
                    "lymphocyte",
                    "macrophage",
                    "mature NK T cell",
                    "monocyte",
                    "mucosal invariant T cell",
                    "naive B cell",
                    "naive thymus-derived CD4-positive, alpha-beta T cell",
                    "naive thymus-derived CD8-positive, alpha-beta T cell",
                    "native cell",
                    "natural killer cell",
                    "non-classical monocyte",
                    "plasma cell",
                    "plasmablast",
                    "plasmacytoid dendritic cell",
                    "platelet",
                    "regulatory T cell",
                    "T cell",
                    "T follicular helper cell"
                ]
            elif (SMP=="mrg7"):
                # nCType=48
                target_celltypes = [
                    "B cell",
                    "CD14-low, CD16-positive monocyte",
                    "CD14-positive monocyte",
                    "CD14-positive, CD16-positive monocyte",
                    "CD16-negative, CD56-bright natural killer cell, human",
                    "CD16-positive, CD56-dim natural killer cell, human",
                    "CD4-positive helper T cell",
                    "CD4-positive, alpha-beta cytotoxic T cell",
                    "CD4-positive, alpha-beta memory T cell",
                    "CD4-positive, alpha-beta T cell",
                    "CD8-positive, alpha-beta cytotoxic T cell",
                    "CD8-positive, alpha-beta memory T cell",
                    "CD8-positive, alpha-beta T cell",
                    "central memory CD4-positive, alpha-beta T cell",
                    "central memory CD8-positive, alpha-beta T cell",
                    "classical monocyte",
                    "conventional dendritic cell",
                    "dendritic cell",
                    "dendritic cell, human",
                    "effector CD8-positive, alpha-beta T cell",
                    "effector memory CD4-positive, alpha-beta T cell",
                    "effector memory CD8-positive, alpha-beta T cell",
                    "effector memory CD8-positive, alpha-beta T cell, terminally differentiated",
                    "erythrocyte",
                    "gamma-delta T cell",
                    "hematopoietic stem cell",
                    "IgA plasma cell",
                    "innate lymphoid cell",
                    "leukocyte",
                    "lymphocyte",
                    "macrophage",
                    "mature NK T cell",
                    "memory B cell",
                    "monocyte",
                    "mucosal invariant T cell",
                    "naive B cell",
                    "naive thymus-derived CD4-positive, alpha-beta T cell",
                    "naive thymus-derived CD8-positive, alpha-beta T cell",
                    "native cell",
                    "natural killer cell",
                    "non-classical monocyte",
                    "plasma cell",
                    "plasmablast",
                    "plasmacytoid dendritic cell",
                    "platelet",
                    "regulatory T cell",
                    "T cell",
                    "T follicular helper cell"
                ]
        
        # target_celltypes=list(set(well_represented_celltypes)) #from the "unconstrained" approach (nCellType=59)
        
        ca.safelist_celltypes(target_celltypes)

        adata = ca.annotate_dataset(adata)

        ## plot 3: UMAP plot - "constrained" prediction-based cell types ####
        if (MK_PLT_3==1):
            FN_PLOT_UMAP=RESDIR_PRJ+"3_UMAP_cellType_pred-constrained_"+PRJ+"_"+SMP+".pdf"
            isExist=os.path.exists(FN_PLOT_UMAP)
            if not isExist:
                ## count the number of cell types included in the "constrained" prediction result
                list_ctype=list(set(target_celltypes))
                nCType=len(list_ctype)

                ## make a customized colour palette
                # pal=glasbey.create_palette(palette_size=nCType)
                # pal=glasbey.extend_palette("tab10", palette_size=nCType, colorblind_safe=True, cvd_severity=100)
                colors=glasbey.extend_palette("Dark2", palette_size=nCType, colorblind_safe=True, cvd_severity=100) #default:lightness_bounds=(10,90)
                pal=sns.set_palette(sns.color_palette(colors))

                sc.pl.umap(
                    adata, 
                    color="celltype_hint", 
                    # cmap=sns.color_palette(palette=pal, as_cmap=True),
                    palette=pal, 
                    legend_fontsize=5,
                    s=size_dot
                    )
                # plt.savefig(FN_PLOT_UMAP)
                plt.savefig(FN_PLOT_UMAP, orientation='landscape', bbox_inches='tight', pad_inches=0.1)
                if (os.path.exists(FN_PLOT_UMAP)):
                    print("[DONE] save as PDF: "+FN_PLOT_UMAP)
        

        ## plot 3a: UMAP plot - "constrained" prediction-based cell types ####    
        if (MK_PLT_3a==1):
            ## set a list of cell types used in the "constrained" mode
            list_ctype=target_celltypes
            ## sort the list
            list_ctype.sort(reverse=False)
        
            ## count the number of cell types included in the "constrained" prediction result
            nCType=len(list_ctype) #53

            for ctype in list_ctype:
                ## define a customized function: [order]
                def order(item, list, n=0):
                    if list:
                        if list[0] == item:
                            return n
                        elif len(list) >= 2: # for python 2, use "else:"
                            return order(item, list[1:], n+1)
                
                ## add "zero" before the order number of each cell type if the number is 1-digit
                ordCType=order(item=ctype, list=list_ctype)+1
                if (ordCType <10):
                    ordCType_str="0"+str(order(item=ctype, list=list_ctype)+1)
                else:
                    ordCType_str=str(order(item=ctype, list=list_ctype)+1)

                print("---- selected cell type: ["+ctype+"]"+" ("+ordCType_str+" /"+str(nCType)+") ---")
                
                ## make a customized colour palette
                colors=glasbey.extend_palette("Dark2", palette_size=nCType, colorblind_safe=True, cvd_severity=100) #default:lightness_bounds=(10,90)
                pal=sns.set_palette(sns.color_palette(colors))

                ax = sc.pl.umap(
                    adata,
                    color=["celltype_hint"], groups=[ctype], 
                    # cmap=sns.color_palette(palette=pal, as_cmap=True),
                    palette=pal, 
                    show=False, #Show the plot, do not return axis. #default: None
                    legend_fontsize=5,
                    s=size_dot,
                    na_color='lightgrey'
                    )
                
                ## make another version of "cell type" label - i.e., each space(" ") is replaced by an underscore("_")
                ctype_spl=ctype.split(" ")
                ctype_rpl="_".join(ctype_spl)
                    
                FN_PLOT_UMAP=RESDIR_PRJ+"3a_UMAP_cellType_pred-constrained_"+PRJ+"_"+SMP+"_"+str(nCType)+"-"+ordCType_str+"_"+ctype_rpl+".pdf"
                isExist=os.path.exists(FN_PLOT_UMAP)
                if not isExist:
                    # We can change the 'NA' in the legend that represents all cells outside of the specified groups
                    legend_texts = ax.get_legend().get_texts()
                    # Find legend object whose text is "NA" and change it
                    for legend_text in legend_texts:
                        if legend_text.get_text() == "NA":
                            legend_text.set_text("other cell types")
                            
                            plt.savefig(FN_PLOT_UMAP, orientation='landscape', bbox_inches='tight', pad_inches=0.1)
                            if (os.path.exists(FN_PLOT_UMAP)):
                                print("[DONE] save as PDF: "+FN_PLOT_UMAP)


        ## plot 3b: "constrained" prediction UMAP plot - w/ GOI expression (aggregated version) ####
        if (MK_PLT_3b==1):
            ## using a single GOI ##
            ## set the GOI name
            GOI="GCAR"
            OBS_NAME_GOI="expr_"+GOI
            
            FN_PLOT_UMAP=RESDIR_PRJ+"3b_UMAP_cellType_pred-constrained_"+PRJ+"_"+SMP+"_w"+GOI+".pdf"
            isExist=os.path.exists(FN_PLOT_UMAP)
            if not isExist:
                ## set a list of cell types used in the "constrained" mode
                list_ctype=target_celltypes
                ## sort the list
                list_ctype.sort(reverse=False)
                
                ## count the number of cell types included in the "constrained" prediction result
                nCType=len(list_ctype)

                ## make customized colour palettes
                colors=glasbey.extend_palette("Dark2", palette_size=nCType, colorblind_safe=True, cvd_severity=100) #default:lightness_bounds=(10,90)
                pal=sns.set_palette(sns.color_palette(colors))
                pal_expr={"neg":"lightgrey", "pos":"#bd0026"}
                
                ## make 2 UMAP subplots - w/ different colour palette and/or cmap parameters ##
                ncol=2
                nrow=1
                figsize=5
                wspace=0.5
                # Adapt figure size based on number of rows and columns and added space between them
                # (e.g. wspace between columns)
                fig,axs=plt.subplots(
                    nrows=nrow, ncols=ncol,
                    figsize=(ncol * figsize + (ncol - 1) * wspace * figsize, nrow * figsize)
                )
                plt.subplots_adjust(wspace=wspace)

                ## UMAP plot (left side) ##
                umap1=sc.pl.umap(adata, color=OBS_NAME_GOI, groups=["pos"], palette=pal_expr_cat, legend_fontsize=5, s=size_dot, ax=axs[0], show=False)
                ## We can change the 'NA' in the legend that represents all cells outside of the specified groups
                legend_texts = umap1.get_legend().get_texts()
                ## Find legend object whose text is "NA" and change it
                for legend_text in legend_texts:
                    if legend_text.get_text() == "NA":
                        legend_text.set_text("neg")

                        ## UMAP plot (right side) ##
                        sc.pl.umap(adata, color="celltype_hint", palette=pal, legend_fontsize=5, s=size_dot, ax=axs[1])

                        plt.savefig(FN_PLOT_UMAP, orientation='landscape', bbox_inches='tight', pad_inches=0.1)
                        if (os.path.exists(FN_PLOT_UMAP)):
                            print("[DONE] save as PDF: "+FN_PLOT_UMAP)


        ## plot 3c: "constrained" prediction UMAP plot - w/ GOI expression (for each cell type) ####
        if (MK_PLT_3c==1):
            ## using a single GOI ##
            ## set the GOI name
            GOI="GCAR"
            OBS_NAME_GOI="expr_"+GOI
            
            ## set a list of cell types used in the "constrained" mode
            list_ctype=target_celltypes
            ## sort the list
            list_ctype.sort(reverse=False)     

            ## count the number of cell types included in the "unconstrained" prediction result
            nCType=len(list_ctype)
            
            for ctype in list_ctype:
                ## define a customized function: [order]
                def order(item, list, n=0):
                    if list:
                        if list[0] == item:
                            return n
                        elif len(list) >= 2: # for python 2, use "else:"
                            return order(item, list[1:], n+1)
                
                ## add "zero" before the order number of each cell type if the number is 1-digit
                ordCType=order(item=ctype, list=list_ctype)+1
                if (ordCType <10):
                    ordCType_str="0"+str(order(item=ctype, list=list_ctype)+1)
                else:
                    ordCType_str=str(order(item=ctype, list=list_ctype)+1)

                print("---- [expr_"+GOI+"] + selected cell type: ["+ctype+"]"+" ("+ordCType_str+" /"+str(nCType)+") ---")                     

                ## make customized colour palettes
                colors=glasbey.extend_palette("Dark2", palette_size=nCType, colorblind_safe=True, cvd_severity=100) #default:lightness_bounds=(10,90)
                pal=sns.set_palette(sns.color_palette(colors))
                pal_expr={"neg":"lightgrey", "pos":"#bd0026"}          

                ## make 2 UMAP subplots - w/ different colour palette and/or cmap parameters ##
                ncol=2
                nrow=1
                figsize=5
                wspace=0.5
                # Adapt figure size based on number of rows and columns and added space between them
                # (e.g. wspace between columns)
                fig,axs=plt.subplots(
                    nrows=nrow, ncols=ncol,
                    figsize=(ncol * figsize + (ncol - 1) * wspace * figsize, nrow * figsize)
                )
                plt.subplots_adjust(wspace=wspace)            
                
                ## make another version of "cell type" label - i.e., each space(" ") is replaced by an underscore("_")
                ctype_spl=ctype.split(" ")
                ctype_rpl="_".join(ctype_spl)
                
                FN_PLOT_UMAP=RESDIR_PRJ+"3c_UMAP_cellType_pred-constrained_"+PRJ+"_"+SMP+"_w"+GOI+"_"+str(nCType)+"-"+ordCType_str+"_"+ctype_rpl+".pdf"
                isExist=os.path.exists(FN_PLOT_UMAP)
                if not isExist:         
                    ## UMAP plot #1 (left side) ##
                    umap1=sc.pl.umap(adata, 
                                    color=OBS_NAME_GOI, groups=["pos"], 
                                    palette=pal_expr_cat, 
                                    legend_fontsize=5, 
                                    s=size_dot, 
                                    ax=axs[0], 
                                    show=False
                                    )
                    ## We can change the 'NA' in the legend that represents all cells outside of the specified groups
                    legend_texts = umap1.get_legend().get_texts()
                    ## Find legend object whose text is "NA" and change it
                    for legend_text in legend_texts:
                        if legend_text.get_text() == "NA":
                            legend_text.set_text("neg")

                            ## UMAP plot #2 (right side) ##
                            umap2=sc.pl.umap(adata,
                                            color=["celltype_hint"], groups=[ctype], 
                                            palette=pal, 
                                            show=False, #Show the plot, do not return axis. #default: None
                                            legend_fontsize=5,
                                            s=size_dot,
                                            na_color='lightgrey',
                                            ax=axs[1]
                                            )                                
                            # We can change the 'NA' in the legend that represents all cells outside of the specified groups
                            legend_texts = umap2.get_legend().get_texts()
                            # Find legend object whose text is "NA" and change it
                            for legend_text in legend_texts:
                                if legend_text.get_text() == "NA":
                                    legend_text.set_text("other cell types")

                            plt.savefig(FN_PLOT_UMAP, orientation='landscape', bbox_inches='tight', pad_inches=0.1)
                            ## checkpoint ##
                            if (os.path.exists(FN_PLOT_UMAP)):
                                print("[DONE] save as PDF: "+FN_PLOT_UMAP)


    ## save the results ## 
        ## save the constrained prediction-labelled anndata as h5ad file
        if (SAVE_H5AD_CONSTRAIN==1):
            # pip install hdf5plugin
            # conda activate hdf5
            import hdf5plugin
            
            # PATH_DATA=DATADIR+PRJ+"_"+SMP+"_afterUMAP"+".h5ad"
            if (TCELL_ONLY==0):
                FN_H5AD_NEW=DATADIR+PRJ+"_"+SMP+"_afterUMAP"+"_SCim_constrain"+"_res"+res_clst_str+".h5ad"
            elif (TCELL_ONLY==1):
                FN_H5AD_NEW=DATADIR+PRJ+"_"+SMP+"_afterUMAP"+"_SCim_constrain"+"_VDJonly"+"_res"+res_clst_str+".h5ad"
            if not (os.path.exists(FN_H5AD_NEW)):
                adata.write_h5ad(
                    # anndata=adata,
                    filename=FN_H5AD_NEW
                    # compression=hdf5plugin.FILTERS["zstd"]
                )
                ## checkpoint ##
                if (os.path.exists(FN_H5AD_NEW)):
                    print("[DONE] save the constrained prediction-labelled anndata as h5ad: "+FN_H5AD_NEW)
                else:
                    print("[CHECK] constrained prediction-labelled h5ad file did not create successfully: "+FN_H5AD_NEW)
        
        ## save the specific column of the constrained prediction-labelled anndata as CSV file
        if (SAVE_CSV_CONSTRAIN==1):
            ## usage ##
            # dirname = Name of the directory to which to export.
            # skip_data = Skip the data matrix X.
            
            if (PRED_DUB!=1):
                OUTDIR_CSV=RESDIR_PRJ+"adata-to-csv_constrain_"+sr_dtype+"_res"+res_clst_str+"/"
            else:
                OUTDIR_CSV=RESDIR_PRJ+"adata-to-csv_constrain_"+sr_dtype+"_res"+res_clst_str+"_expRate"+rate_exp_doub_str+"/"
            
            isExist=os.path.exists(OUTDIR_CSV)
            if not isExist:
                    os.makedirs(OUTDIR_CSV)
                    print("[DONE] create a directory: ", OUTDIR_CSV, sep="")
            else:
                print("[SKIP] directory exists: ", OUTDIR_CSV, sep="")
            
            ## save as CSV
            adata.write_csvs(dirname=OUTDIR_CSV, skip_data=True, sep=',') #worked well
            # adata[adata.obs.celltype_hint].write_csvs(dirname=OUTDIR_CSV, skip_data=True, sep=',')

            FN_CSV=OUTDIR_CSV+"obs.csv"
            ## checkpoint ##
            if (os.path.exists(FN_CSV)):
                print("[DONE] save the constrained prediction-labelled anndata as CSV: "+FN_CSV)
            else:
                print("[CHECK] constrained prediction-labelled CSV file did not create successfully: "+FN_CSV)



    ## STEP 4 | Annotation QC ####
    ## - source: [https://genentech.github.io/scimilarity/notebooks/cell_annotation_tutorial.html#Annotation-QC]
    
    # cell annotation also computes QC metrics for our annotations. 
    
    # One of which, min_dist, represents the minimum distance between a cell in the query dataset and all cells in the training set. 
    # The greater min_dist, (i.e., the further away from what the model has seen before) the less confidence we have in the model’s prediction.

    # Note, for different applications and questions different min_dist ranges have different implications.

    if (ANNO_QC==1):
        print("--- STEP 4 | annotation QC ---")

        ## check the "min_dist" values stored in the [adata]
        ## - note: the "min_dist" values are only calculated in the "constrained" approach
        # print(adata.obs["min_dist"])
        
        ## plot 4: UMAP plot - annotation QC ####
        if (MK_PLT_4==1):
            if (ANNO_CONSTRAIN==1):
                FN_PLOT_UMAP=RESDIR_PRJ+"4_UMAP_annotationQC_"+PRJ+"_"+SMP+".pdf"
                isExist=os.path.exists(FN_PLOT_UMAP)
                if not isExist:
                    sc.pl.umap(
                        adata, 
                        color="min_dist",
                        cmap=sns.blend_palette(["lightgrey", sns.xkcd_rgb["cobalt"]], as_cmap=True),
                        vmax=0.1,
                        legend_fontsize=5,
                        s=size_dot
                        )
                    plt.savefig(FN_PLOT_UMAP)
                    if (os.path.exists(FN_PLOT_UMAP)):
                        print("[DONE] save as PDF: "+FN_PLOT_UMAP)
            else:
                print("[CHECK] 'ANNO_CONSTRAIN' running switch must be turned ON for this step")

        ## plot 4a: UMAP plot - subplots w/ GCAR-expressing cells & "unconstrained" cell type prediction
        if (MK_PLT_4a==1):
            ## using a single GOI ##
            ## set the GOI name
            GOI="GCAR"
            OBS_NAME_GOI="expr_"+GOI

            ## make a list of cell types predicted from the "unconstrained" mode
            list_ctype=list(set(well_represented_celltypes))
            ## sort the list
            list_ctype.sort(reverse=False)
            
            ## count the number of cell types included in the "constrained" prediction result
            nCType=len(list_ctype)

            for ctype in list_ctype:       
                ## define a customized function: [order]
                def order(item, list, n=0):
                    if list:
                        if list[0] == item:
                            return n
                        elif len(list) >= 2: # for python 2, use "else:"
                            return order(item, list[1:], n+1)
                
                ## add "zero" before the order number of each cell type if the number is 1-digit
                ordCType=order(item=ctype, list=list_ctype)+1
                if (ordCType <10):
                    ordCType_str="0"+str(order(item=ctype, list=list_ctype)+1)
                else:
                    ordCType_str=str(order(item=ctype, list=list_ctype)+1)

                print("---- plot 4a: [expr_"+GOI+"] + selected cell type: ["+ctype+"]"+" ("+ordCType_str+" /"+str(nCType)+") ---")           

                ## make customized colour palettes
                colors=glasbey.extend_palette("Dark2", palette_size=nCType, colorblind_safe=True, cvd_severity=100) #default:lightness_bounds=(10,90)
                pal=sns.set_palette(sns.color_palette(colors))
                pal_expr={"neg":"lightgrey", "pos":"#bd0026"}
                
                ## make a customized color map
                ## - As with the convention in matplotlib, every continuous colormap has a reversed version, which has the suffix "_r": ["rocket"],["rocket_r"]
                # pal_dist=sns.blend_palette(["lightgrey", sns.xkcd_rgb["cobalt"]], as_cmap=True) #previously: ["scarlet"]
                pal_dist=sns.color_palette("crest_r", as_cmap=True)


                ## make another version of "cell type" label - i.e., each space(" ") is replaced by an underscore("_")
                ctype_spl=ctype.split(" ")
                ctype_rpl="_".join(ctype_spl)

                FN_PLOT_UMAP=RESDIR_PRJ+"4a_UMAP_annoQC_"+PRJ+"_"+SMP+"_w"+GOI+"_ucPred"+"_"+str(nCType)+"-"+ordCType_str+"_"+ctype_rpl+".png" #changed to PNG(default) due to large file size issues
                isExist=os.path.exists(FN_PLOT_UMAP)
                if not isExist:
                    ## make 4 UMAP subplots - w/ different colour palette and/or cmap parameters ##
                    ncol=3
                    nrow=1
                    figsize=5
                    wspace=0.5

                    # Adapt figure size based on number of rows and columns and added space between them
                    # (e.g. wspace between columns)
                    fig,axs=plt.subplots(
                        nrows=nrow, ncols=ncol,
                        figsize=(ncol * figsize + (ncol - 1) * wspace * figsize, nrow * figsize)
                    )
                    plt.subplots_adjust(wspace=wspace)

                    ## UMAP plot #1 (left/top) ##
                    if (SMP!="mrg"):
                        if (SMP!="Enriched_Apheresis" and SMP!="GCAR1_Product"):
                            mytitle="["+SMP+"] Seurat clusters (nClst="+nClst_str+")"
                        elif (SMP=="Enriched_Apheresis"):
                            mytitle="["+"Enriched"+"] Seurat clusters (nClst="+nClst_str+")"
                        elif (SMP=="GCAR1_Product"):
                            mytitle="["+"Harvest"+"] Seurat clusters (nClst="+nClst_str+")"
                        umap1=sc.pl.umap(adata, 
                                    color=cname_clst, 
                                    palette=pal_clst,
                                    legend_fontsize=5,
                                    s=size_dot,
                                    ax=axs[0],
                                    title=mytitle
                                    )

                    elif (SMP=="mrg" or SMP=="mrg7"):
                        mytitle="GCAR_PT01 (nSmp="+nSmp_str+")"
                        umap1=sc.pl.umap(adata, 
                                        color="orig.ident", #orig.ident #celltype_raw
                                        palette=pal_smp,
                                        legend_fontsize=5,
                                        s=size_dot,
                                        ax=axs[0],
                                        title=mytitle
                                        )

                    ## UMAP plot #2 (right/top) ##
                    umap2=sc.pl.umap(adata, 
                                    color=OBS_NAME_GOI, groups=["pos"], 
                                    palette=pal_expr_cat, 
                                    legend_fontsize=5, 
                                    s=size_dot, 
                                    ax=axs[1], 
                                    show=False
                                    )
                    ## We can change the 'NA' in the legend that represents all cells outside of the specified groups
                    legend_texts = umap2.get_legend().get_texts()
                    ## Find legend object whose text is "NA" and change it
                    for legend_text in legend_texts:
                        if legend_text.get_text() == "NA":
                            legend_text.set_text("neg")

                    ## UMAP plot #3 (left/bottom) ##
                    umap3=sc.pl.umap(adata[adata.obs.predictions_unconstrained.isin(well_represented_celltypes)],
                                    color=["predictions_unconstrained"], groups=[ctype], 
                                    palette=pal, 
                                    show=False, #Show the plot, do not return axis. #default: None
                                    legend_fontsize=5,
                                    s=size_dot,
                                    na_color='lightgrey',
                                    ax=axs[2]
                                    )
                    # We can change the 'NA' in the legend that represents all cells outside of the specified groups
                    legend_texts = umap3.get_legend().get_texts()
                    # Find legend object whose text is "NA" and change it
                    for legend_text in legend_texts:
                        if legend_text.get_text() == "NA":
                            legend_text.set_text("other cell types")
                    
                    plt.savefig(FN_PLOT_UMAP, orientation='landscape', bbox_inches='tight', pad_inches=0.1)
                    if (os.path.exists(FN_PLOT_UMAP)):
                        print("[DONE] save as PDF: "+FN_PLOT_UMAP)
            
            ## merge PDF files ## - skip for now
            # # import PyPDF2
            # DIR=RESDIR_PRJ
            # FN_PDF_MERGE=RESDIR_PRJ+"4a_UMAP_annoQC_"+PRJ+"_"+SMP+"_w"+GOI+"_ucPred"+"_nCellType"+str(nCType)+".pdf"
            # START="4a_UMAP_annoQC_"+PRJ+"_"+SMP+"_w"+GOI+"_ucPred"+"_"+str(nCType)+"-"
            # # END="_vmax"+vm+".png" #".pdf"

            # ## change the current working directory
            # print("Current working directory: {0}".format(os.getcwd()))
            # os.chdir(DIR)
            # print("New working directory: {0}".format(os.getcwd()))

            # merger=PyPDF2.PdfMerger()

            # for file in os.listdir(DIR):
            #     if file.startswith(START) & file.endswith(END):
            #             merger.append(file)
            # merger.write(FN_PDF_MERGE)

            # isExist=os.path.exists(FN_PDF_MERGE)
            # if isExist:
            #     print("[DONE] merge PDF files: "+FN_PDF_MERGE+" ("+nClst_str+" pages)")

            #     for fname in os.listdir(DIR):
            #         PATH_PLOT_UMAP=os.path.join(DIR,fname)
            #         if fname.startswith(START) & fname.endswith(END):
            #             os.remove(PATH_PLOT_UMAP)
            #             print("[DONE] remove intermediate PDF files: "+PATH_PLOT_UMAP+")")               


    ## plot 4b: UMAP plot - subplots w/ GCAR-expressing cells & "unconstrained" cell type prediction
        if (MK_PLT_4b==1):
            if (ANNO_CONSTRAIN==1):
                ## using a single GOI ##
                ## set the GOI name
                GOI="GCAR"
                OBS_NAME_GOI="expr_"+GOI
                
                ## set a list of cell types used in the "constrained" mode
                list_ctype=target_celltypes
                ## sort the list
                list_ctype.sort(reverse=False)
                
                ## count the number of cell types included in the "constrained" prediction result
                nCType=len(list_ctype)

                for ctype in list_ctype:
                    ## define a customized function: [order]
                    def order(item, list, n=0):
                        if list:
                            if list[0] == item:
                                return n
                            elif len(list) >= 2: # for python 2, use "else:"
                                return order(item, list[1:], n+1)
                    
                    ## add "zero" before the order number of each cell type if the number is 1-digit
                    ordCType=order(item=ctype, list=list_ctype)+1
                    if (ordCType <10):
                        ordCType_str="0"+str(order(item=ctype, list=list_ctype)+1)
                    else:
                        ordCType_str=str(order(item=ctype, list=list_ctype)+1)

                    print("---- plot 4b: [expr_"+GOI+"] + selected cell type: ["+ctype+"]"+" ("+ordCType_str+" /"+str(nCType)+") ---")
                    
                    ## make customized colour palettes
                    colors=glasbey.extend_palette("Dark2", palette_size=nCType, colorblind_safe=True, cvd_severity=100) #default:lightness_bounds=(10,90)
                    pal=sns.set_palette(sns.color_palette(colors))
                    pal_expr={"neg":"lightgrey", "pos":"#bd0026"}
                                
                    ## make a customized color map
                    ## - As with the convention in matplotlib, every continuous colormap has a reversed version, which has the suffix "_r": ["rocket"],["rocket_r"]
                    # pal_dist=sns.blend_palette(["lightgrey", sns.xkcd_rgb["cobalt"]], as_cmap=True) #previously: ["scarlet"]
                    pal_dist=sns.color_palette("crest_r", as_cmap=True) 

                    ## make a list of [vmax] variables
                    list_vmax=["0.1","0.05","0.01"]

                    for vm in list_vmax:
                        ## make another version of "cell type" label - i.e., each space(" ") is replaced by an underscore("_")
                        ctype_spl=ctype.split(" ")
                        ctype_rpl="_".join(ctype_spl)

                        FN_PLOT_UMAP=RESDIR_PRJ+"4b_UMAP_annoQC_"+PRJ+"_"+SMP+"_w"+GOI+"_cPred"+"_"+str(nCType)+"-"+ordCType_str+"_"+ctype_rpl+"_vmax"+vm+".png"
                        isExist=os.path.exists(FN_PLOT_UMAP)
                        if not isExist:
                            ## make 4 UMAP subplots - w/ different colour palette and/or cmap parameters ##
                            ncol=2 #3
                            nrow=2 #1
                            figsize=5
                            wspace=0.5

                            # Adapt figure size based on number of rows and columns and added space between them
                            # (e.g. wspace between columns)
                            fig,axs=plt.subplots(
                                nrows=nrow, ncols=ncol,
                                figsize=(ncol * figsize + (ncol - 1) * wspace * figsize, nrow * figsize)
                            )
                            plt.subplots_adjust(wspace=wspace)

                            ## UMAP plot #1 (left/top) ##
                            if (SMP!="mrg"):
                                if (SMP!="Enriched_Apheresis" and SMP!="GCAR1_Product"):
                                    mytitle="["+SMP+"] Seurat clusters (nClst="+nClst_str+")"
                                elif (SMP=="Enriched_Apheresis"):
                                    mytitle="["+"Enriched"+"] Seurat clusters (nClst="+nClst_str+")"
                                elif (SMP=="GCAR1_Product"):
                                    mytitle="["+"Harvest"+"] Seurat clusters (nClst="+nClst_str+")"
                                umap1=sc.pl.umap(adata, 
                                            color=cname_clst, #"seurat_clusters"
                                            palette=pal_clst,
                                            legend_fontsize=5,
                                            s=size_dot,
                                            ax=axs[0,0],
                                            title=mytitle
                                            )

                            elif (SMP=="mrg" or SMP=="mrg7"):
                                mytitle="GCAR_PT01 (nSmp="+nSmp_str+")"
                                umap1=sc.pl.umap(adata, 
                                                color="orig.ident", #orig.ident #celltype_raw
                                                palette=pal_smp,
                                                legend_fontsize=5,
                                                s=size_dot,
                                                ax=axs[0,0],
                                                title=mytitle
                                                )

                            ## UMAP plot #2 (right/top) ##
                            umap2=sc.pl.umap(adata, 
                                            color=OBS_NAME_GOI, 
                                            groups=["pos"], 
                                            palette=pal_expr_cat, 
                                            legend_fontsize=5, 
                                            s=size_dot, 
                                            ax=axs[0,1], 
                                            show=False
                                            )
                            ## We can change the 'NA' in the legend that represents all cells outside of the specified groups
                            legend_texts = umap2.get_legend().get_texts()
                            ## Find legend object whose text is "NA" and change it
                            for legend_text in legend_texts:
                                if legend_text.get_text() == "NA":
                                    legend_text.set_text("neg")                            
                                    
                            ## UMAP plot #3 (left/bottom) ##
                            umap3=sc.pl.umap(adata, 
                                            color=["celltype_hint"], groups=[ctype], 
                                            palette=pal, 
                                            show=False, #Show the plot, do not return axis. #default: None
                                            legend_fontsize=5,
                                            s=size_dot,
                                            na_color='lightgrey',
                                            ax=axs[1,1] #changed as some cell type names are too long
                                            )                        
                                            
                            # We can change the 'NA' in the legend that represents all cells outside of the specified groups
                            legend_texts = umap3.get_legend().get_texts()
                            # Find legend object whose text is "NA" and change it
                            for legend_text in legend_texts:                            
                                if legend_text.get_text() == "NA":
                                    legend_text.set_text("other cell types")
                                

                            ## UMAP plot #4 (right/bottom) ##
                            mytitle="min_dist (vmax="+vm+")"
                            umap4=sc.pl.umap(adata[adata.obs.celltype_hint.isin(target_celltypes)], #success!
                                            color="min_dist",
                                            mask_obs=(adata[adata.obs.celltype_hint.isin(target_celltypes)].obs.celltype_hint==ctype), #success!
                                            cmap=pal_dist, 
                                            vmax=vm, 
                                            legend_fontsize=5, 
                                            s=size_dot, 
                                            ax=axs[1,0],
                                            title=mytitle
                                            )                                
                                        
                            plt.savefig(FN_PLOT_UMAP, orientation='landscape', bbox_inches='tight', pad_inches=0.1)
                            if (os.path.exists(FN_PLOT_UMAP)):
                                print("[DONE] save as PNG: "+FN_PLOT_UMAP)
                
                ## merge PDF files ##
                for vm in list_vmax:
                    # import PyPDF2
                    DIR=RESDIR_PRJ
                    FN_PDF_MERGE=RESDIR_PRJ+"4b_UMAP_annoQC_"+PRJ+"_"+SMP+"_w"+GOI+"_cPred"+"_nCellType"+str(nCType)+"_vmax"+vm+".pdf"
                    START="4b_UMAP_annoQC_"+PRJ+"_"+SMP+"_w"+GOI+"_cPred"+"_"+str(nCType)+"-"
                    END="_vmax"+vm+".png" #".pdf"

                    ## change the current working directory
                    print("Current working directory: {0}".format(os.getcwd()))
                    os.chdir(DIR)
                    print("New working directory: {0}".format(os.getcwd()))

                    merger=PyPDF2.PdfMerger()

                    for file in os.listdir(DIR):
                        if file.startswith(START) & file.endswith(END):
                                merger.append(file)
                    merger.write(FN_PDF_MERGE)

                    isExist=os.path.exists(FN_PDF_MERGE)
                    if isExist:
                        print("[DONE] merge PDF files: "+FN_PDF_MERGE+" ("+nClst_str+" pages)")

                        for fname in os.listdir(DIR):
                            PATH_PLOT_UMAP=os.path.join(DIR,fname)
                            if fname.startswith(START) & fname.endswith(END):
                                os.remove(PATH_PLOT_UMAP)
                                print("[DONE] remove intermediate PDF files: "+PATH_PLOT_UMAP+")")
            else:
                print("[CHECK] 'ANNO_CONSTRAIN' running switch must be turned ON for this step")
    # ===


#### marker gene expression visualization ####
if (VIS_EXPR==1):
    # import scanpy as sc #loaded at the beginning of this code
    from matplotlib.pyplot import rc_context
    
    ## load the h5ad file
    if (PRJ=="GCAR1_PT01"):
        if (TCELL_ONLY!=1):
            PATH_DATA=DATADIR+PRJ+"_"+SMP+"_afterUMAP"+"_res"+res_clst_str+".h5ad"
        elif (TCELL_ONLY==1):
            PATH_DATA=DATADIR+PRJ+"_"+SMP+"_afterUMAP"+"_raw_vdj"+"_res"+res_clst_str+".h5ad"
        if os.path.exists(PATH_DATA):
            print("[START] load the h5ad file: "+PATH_DATA)
        else:
            print("[CHECK] h5ad file does not exist: "+PATH_DATA)
    
    adata=sc.read(PATH_DATA)
    print(adata.shape) #checkpoint
    print("[END] load the h5ad file: "+PATH_DATA)
    
    # ---
    ## add a column to [obs]: [sample_new] ##
    ## create a list of conditions
    conditions=[
        (adata.obs["orig.ident"]=="Enriched_Apheresis"),
        (adata.obs["orig.ident"]=="GCAR1_Product"),
        (adata.obs["orig.ident"]=="D08"),
        (adata.obs["orig.ident"]=="D15"),
        (adata.obs["orig.ident"]=="D22"),
        (adata.obs["orig.ident"]=="D28"),
        (adata.obs["orig.ident"]=="D42")
    ]
    ## create a list of the "preferred" values to assign for each condition
    values=["Enriched_Apheresis","GCAR1_Product","D08","D15","D22","D28","D42"]
    
    ## create a new column to assign values
    adata.obs["sample_new"] = np.select(conditions,values)
    print(set(adata.obs["sample_new"]))
    print(pd.DataFrame(adata.obs).head()) #checkpoint
    print(pd.DataFrame(adata.obs).tail()) #checkpoint
    # ---
    
    ## subset the [adata] ####
    if (TCELL_ONLY!=1):
        if (EXPR_PLT=="all"):
            adata_input=adata
        elif (EXPR_PLT=="VDJ"):
            adata_input=adata[adata.obs.cell_VDJ=="VDJ"]
        elif (EXPR_PLT=="GCARpos"):
            adata_input=adata[adata.obs.expr_GCAR=="pos"]
        elif (EXPR_PLT=="GCARneg"):
            adata_input=adata[adata.obs.expr_GCAR=="neg"]
        elif (EXPR_PLT=="VDJ_GCARpos"):
            adata_input=adata[adata.obs.cell_VDJ=="VDJ"]
            adata_input=adata_input[adata_input.obs.expr_GCAR=="pos"]
        elif (EXPR_PLT=="VDJ_GCARneg"):
            adata_input=adata[adata.obs.cell_VDJ=="VDJ"]
            adata_input=adata_input[adata_input.obs.expr_GCAR=="neg"]
    
    elif (TCELL_ONLY==1):
        if (EXPR_PLT=="VDJ"):
            adata_input=adata
        elif (EXPR_PLT=="VDJ_GCARpos"):
            adata_input=adata[adata.obs.expr_GCAR=="pos"]
        elif (EXPR_PLT=="VDJ_GCARneg"):
            adata_input=adata[adata.obs.expr_GCAR=="neg"]


    ## count the number of the Seurat Clusters
    list_clst=set(adata_input.obs["seurat_clusters"])
    nClst=len(list_clst)
    nClst_str=str(nClst)

    # ---
    if (EXPR_PLT=="GCARpos" or EXPR_PLT=="GCARneg"):
        ## add a column to [obs]: [sample_new_GCAR] ##        
        ## create a list of conditions
        conditions=[
            (adata_input.obs["orig.ident"]=="Enriched_Apheresis"),
            (adata_input.obs["orig.ident"]=="GCAR1_Product"),
            (adata_input.obs["orig.ident"]=="D08"),
            (adata_input.obs["orig.ident"]=="D15"),
            (adata_input.obs["orig.ident"]=="D22"),
            (adata_input.obs["orig.ident"]=="D28"),
            (adata_input.obs["orig.ident"]=="D42")
        ]

        ## create a list of the "preferred" values to assign for each condition
        if (EXPR_PLT=="GCARpos"):
            values=[
                "Enriched_pos","Harvest_pos","D08_pos","D15_pos","D22_pos","D28_pos","D42_pos",
                ]
        elif (EXPR_PLT=="GCARneg"):
            values=[
                "Enriched_neg","Harvest_neg","D08_neg","D15_neg","D22_neg","D28_neg","D42_neg",
            ]
        
        ## create a new column to assign values
        adata_input.obs["sample_new_GCAR"] = np.select(conditions,values)
        print(set(adata_input.obs["sample_new_GCAR"]))
        print(pd.DataFrame(adata_input.obs).head()) #checkpoint
        print(pd.DataFrame(adata_input.obs).tail()) #checkpoint
        
    elif (EXPR_PLT=="VDJ_GCARpos"):
        ## add a column to [obs]: [sample_new_GCAR] ##        
        ## create a list of conditions
        conditions=[
            (adata_input.obs["orig.ident"]=="Enriched_Apheresis"),
            (adata_input.obs["orig.ident"]=="GCAR1_Product"),
            (adata_input.obs["orig.ident"]=="D08"),
            (adata_input.obs["orig.ident"]=="D15"),
            (adata_input.obs["orig.ident"]=="D22"),
            (adata_input.obs["orig.ident"]=="D28"),
            (adata_input.obs["orig.ident"]=="D42")
        ]

        ## create a list of the "preferred" values to assign for each condition
        values=[
            "Enriched_VDJ_pos","Harvest_VDJ_pos","D08_VDJ_pos","D15_VDJ_pos","D22_VDJ_pos","D28_VDJ_pos","D42_VDJ_pos",
            ]
        
        ## create a new column to assign values
        adata_input.obs["sample_new_VDJ_GCAR"] = np.select(conditions,values)
        print(set(adata_input.obs["sample_new_VDJ_GCAR"]))
        print(pd.DataFrame(adata_input.obs).head()) #checkpoint
        print(pd.DataFrame(adata_input.obs).tail()) #checkpoint
    
    elif (EXPR_PLT=="VDJ_GCARneg"):
        ## add a column to [obs]: [sample_new_GCAR] ##        
        ## create a list of conditions
        conditions=[
            (adata_input.obs["orig.ident"]=="Enriched_Apheresis"),
            (adata_input.obs["orig.ident"]=="GCAR1_Product"),
            (adata_input.obs["orig.ident"]=="D08"),
            (adata_input.obs["orig.ident"]=="D15"),
            (adata_input.obs["orig.ident"]=="D22"),
            (adata_input.obs["orig.ident"]=="D28"),
            (adata_input.obs["orig.ident"]=="D42")
        ]

        ## create a list of the "preferred" values to assign for each condition
        values=[
            "Enriched_VDJ_neg","Harvest_VDJ_neg","D08_VDJ_neg","D15_VDJ_neg","D22_VDJ_neg","D28_VDJ_neg","D42_VDJ_neg",
            ]
        
        ## create a new column to assign values
        adata_input.obs["sample_new_VDJ_GCAR"] = np.select(conditions,values)
        print(set(adata_input.obs["sample_new_VDJ_GCAR"]))
        print(pd.DataFrame(adata_input.obs).head()) #checkpoint
        print(pd.DataFrame(adata_input.obs).tail()) #checkpoint
    # ---

    if (EXPR_PLT=="VDJ" or EXPR_PLT=="VDJ_GCARpos" or EXPR_PLT=="VDJ_GCARneg"):
        if (TCELL_ONLY!=1):
            print(pd.DataFrame(adata_input[adata_input.obs.cell_VDJ=="VDJ"].obs).head()) #checkpoint
            print(pd.DataFrame(adata_input[adata_input.obs.cell_VDJ=="VDJ"].obs).tail()) #checkpoint
    
    
    ## set the order of the categorical variables ##
    ## [sample_new] ##
    if (EXPR_PLT=="all"):
        ord_smp=["Enriched","Harvest","D08","D15","D22","D28","D42"]
    elif (EXPR_PLT=="VDJ"):
        ord_smp=["Enriched_VDJ","Harvest_VDJ","D08_VDJ","D15_VDJ","D22_VDJ","D28_VDJ","D42_VDJ"]
    elif (EXPR_PLT=="GCARpos"):
        ord_smp=["Enriched_pos","Harvest_pos","D08_pos","D15_pos","D22_pos","D28_pos","D42_pos"]
    elif (EXPR_PLT=="GCARneg"):
        ord_smp=["Enriched_neg","Harvest_neg","D08_neg","D15_neg","D22_neg","D28_neg","D42_neg"]
    elif (EXPR_PLT=="VDJ_GCARpos"):
        ord_smp=["Enriched_VDJ_pos","Harvest_VDJ_pos","D08_VDJ_pos","D15_VDJ_pos","D22_VDJ_pos","D28_VDJ_pos","D42_VDJ_pos"]
    elif (EXPR_PLT=="VDJ_GCARneg"):
        ord_smp=["Enriched_VDJ_neg","Harvest_VDJ_neg","D08_VDJ_neg","D15_VDJ_neg","D22_VDJ_neg","D28_VDJ_neg","D42_VDJ_neg"]


    ## [seurat_clusters] ##
    ## make a list of "string"-format Seurat Cluster variables
    list_clst=sorted(set(adata_input.obs.seurat_clusters))
    print(str(list_clst)) #checkpoint
    nCluster=len(list_clst)
    ## sort numbers after converting them to "integers"
    list_clst=sorted([int(c) for c in list_clst]) #success!
    print(str(list_clst)) #checkpoint
    ord_clst=[str(clst) for clst in list_clst]
    nCluster=len(ord_clst)
    print(str(ord_clst)) #checkpoint
    if (EXPR_PLT=="all"):
        ord_clst=[str(nCluster) for nCluster in range(nCluster)]
    else:
        ord_clst=[str(clst) for clst in list_clst]
    # ---
    
    ## set the marker gene list
    if (EXPR_MK_DOT==1):
        if (EXPR_MARKER_SET=="CAR-T"):
            ## - source: Figure 4 - Anderson ND et al. Nature Medicine 2023 [PMID:37407840]
            dict_marker_genes={
                ## Stem-like memory ##
                "stem-like_memory": ["HNRNPLL","CCR7","SELL","LEF1","IL7R","TCF7"],
                ## Activation, cytotoxicity, effector function ##
                "activation_cytotoxicity_effector": ["CD28","CD27","GZMK","GZMA","PRF1","CCL4","GZMM","CCL5","NKG7","GZMB","GNLY","LYAR","CXCR3","CXCR4","TXNIP"],
                ## Pre-exhaustion, exhaustion ##
                "pre-exhaustion_exhaustion": ["TIGIT","ENTPD1","HAVCR2","CTLA4","LAG3","PDCD1","EOMES","TOX","PMCH","BATF","PRDM1","CXCR1"],
                ## Resident ##
                "resident": ["MCM5","TNFRSF9","ITGAE","CD69","NR4A2"]
            }
            ## w/o categories ##
            ## - the order of each marker gene will be used "as-is" in [plot 1c]
            dict_marker_genes_woCat={
                "CAR-T": [
                    "HNRNPLL","CCR7","SELL","LEF1","IL7R","TCF7",
                    "CD28","CD27","GZMK","GZMA","PRF1","CCL4","GZMM","CCL5","NKG7","GZMB","GNLY","LYAR","CXCR3","CXCR4","TXNIP",
                    "TIGIT","ENTPD1","HAVCR2","CTLA4","LAG3","PDCD1","EOMES","TOX","PMCH","BATF","PRDM1","CXCR1",
                    "MCM5","TNFRSF9","ITGAE","CD69","NR4A2"
                ]
            }
        
        ## count the number of "values" in the marker gene "dictionary"
        nMarkerCat=len(dict_marker_genes.items())
        nMarkerCat_str=str(nMarkerCat)
        
        ## set default "mytitle" variable
        mytitle=""


    ## set default "groupby" variable
    if (EXPR_PLT=="all"):
        var_grp="sample_new"
    elif (EXPR_PLT=="VDJ"):
        var_grp="sample_new_VDJ"
    elif (EXPR_PLT=="GCARpos" or EXPR_PLT=="GCARneg"):
        var_grp="sample_new_GCAR"
    elif (EXPR_PLT=="VDJ_GCARpos" or EXPR_PLT=="VDJ_GCARneg"):
        var_grp="sample_new_VDJ_GCAR"
    
    #### dot plots ####
    if (EXPR_MK_DOT==1):
        ## option 1: dot plot ##
        ## - source: [https://scanpy-tutorials.readthedocs.io/en/latest/plotting/core.html#dotplot]
        if (EXPR_MK_PLT_1==1):        
            FN_PLOT_DOT=RESDIR_PRJ_EXPR+"1_dot_exprMarker_"+PRJ+"_"+SMP+"_"+EXPR_MARKER_SET+"_nCat"+nMarkerCat_str+"_"+EXPR_PLT+".pdf"
            
            isExist=os.path.exists(FN_PLOT_DOT)
            if not isExist:
                ## set plot dimension
                plt.rcParams["figure.figsize"] = [8, 5] #figure size in inches
                
                ## dot plot ##
                fig=sc.pl.dotplot(
                    adata_input, dict_marker_genes, 
                    groupby=var_grp,
                    categories_order=ord_smp,
                    cmap=pal_expr_dot,
                    dendrogram=False
                    )
                
                ## save as PNG
                plt.savefig(FN_PLOT_DOT, orientation='landscape', bbox_inches='tight', pad_inches=0.1)
                print("[DONE] save as PNG: "+FN_PLOT_DOT)
        # ---

        # ---
        ## option 1a: dot plot (w/ Dendrograms) ##
        if (EXPR_MK_PLT_1a==1):        
            FN_PLOT_DOT=RESDIR_PRJ_EXPR+"1a_dot_exprMarker_"+PRJ+"_"+SMP+"_"+EXPR_MARKER_SET+"_nCat"+nMarkerCat_str+"_wDendro"+"_"+EXPR_PLT+".pdf"
            
            isExist=os.path.exists(FN_PLOT_DOT)
            if not isExist:
                ## set plot dimension
                plt.rcParams["figure.figsize"] = [8.5, 5] #figure size in inches
                
                ## dot plot ##
                fig=sc.pl.dotplot(
                    adata_input, dict_marker_genes, 
                    groupby=var_grp, 
                    cmap=pal_expr_dot,
                    # title=mytitle,
                    dendrogram=True
                    )
                
                ## save as PNG
                plt.savefig(FN_PLOT_DOT, orientation='landscape', bbox_inches='tight', pad_inches=0.1)
                print("[DONE] save as PNG: "+FN_PLOT_DOT)
        # ---
        
        # ---
        ## option 1b: dot plot (w/ Dendrograms) ##
        if (EXPR_MK_PLT_1b==1):        
            FN_PLOT_DOT=RESDIR_PRJ_EXPR+"1b_dot_exprMarker_"+PRJ+"_"+SMP+"_"+EXPR_MARKER_SET+"_nCat"+nMarkerCat_str+"_tx"+"_"+EXPR_PLT+".pdf"
            
            isExist=os.path.exists(FN_PLOT_DOT)
            if not isExist:
                ## set plot dimension
                plt.rcParams["figure.figsize"] = [5.5, 8] #figure size in inches
                
                ## dot plot ##
                mytitle="GCAR_PT01/"+SMP+": ["+EXPR_MARKER_SET+"] - {nCategory: "+nMarkerCat_str+"}"
                fig=sc.pl.dotplot(
                    adata_input, dict_marker_genes, 
                    groupby=var_grp, 
                    categories_order=ord_smp,
                    cmap=pal_expr_dot,
                    title=mytitle,
                    swap_axes=True, #False(default)
                    dendrogram=False
                    )
                
                ## save as PNG
                plt.savefig(FN_PLOT_DOT, orientation='landscape', bbox_inches='tight', pad_inches=0.1)
                print("[DONE] save as PNG: "+FN_PLOT_DOT)
        # ---
        
        # ---
        ## option 1c: dot plot (w/ Dendrograms) ##
        if (EXPR_MK_PLT_1c==1):
            if (EXPR_MARKER_SET=="FlowPanel"):
                FN_PLOT_DOT=RESDIR_PRJ_EXPR+"1c_dot_exprMarker_"+PRJ+"_"+SMP+"_"+EXPR_MARKER_SET+"_nCat"+nMarkerCat_str+"_woCat"+"_wDendro"+"_"+EXPR_PLT+".pdf"
                
                isExist=os.path.exists(FN_PLOT_DOT)
                if not isExist:
                    ## set plot dimension
                    plt.rcParams["figure.figsize"] = [5.5, 8] #figure size in inches
                    
                    ## dot plot ##
                    mytitle="GCAR_PT01/"+SMP+": ["+EXPR_MARKER_SET+"] - {nCategory: "+nMarkerCat_str+"}"
                    fig=sc.pl.dotplot(
                        adata_input, dict_marker_genes_woCat, 
                        groupby=var_grp,
                        categories_order=ord_smp,
                        cmap=pal_expr_dot,
                        title=mytitle,
                        swap_axes=False,
                        dendrogram=True
                        )
                    
                    ## save as PNG
                    plt.savefig(FN_PLOT_DOT, orientation='landscape', bbox_inches='tight', pad_inches=0.1)
                    print("[DONE] save as PNG: "+FN_PLOT_DOT)
        # ---
        
        # ---
        ## option 1d: dot plot (w/ Dendrograms) ##
        if (EXPR_MK_PLT_1d==1):
            if (EXPR_MARKER_SET=="FlowPanel"):
                FN_PLOT_DOT=RESDIR_PRJ_EXPR+"1d_dot_exprMarker_"+PRJ+"_"+SMP+"_"+EXPR_MARKER_SET+"_nCat"+nMarkerCat_str+"_wNumCell"+"_"+EXPR_PLT+".pdf"
                
                isExist=os.path.exists(FN_PLOT_DOT)
                if not isExist:
                    ## set plot dimension
                    plt.rcParams["figure.figsize"] = [9, 5] #figure size in inches
                    
                    ## dot plot ##
                    fig=sc.pl.dotplot(
                        adata_input, dict_marker_genes_woCat, 
                        groupby=var_grp, 
                        categories_order=ord_smp,
                        cmap=pal_expr_dot,
                        title=mytitle,
                        swap_axes=False,
                        dendrogram=False,
                        return_fig=True #important
                        )
                    ## add bar plots showing the "number of cells" in [groupby] category.
                    fig.add_totals()

                    ## save as PNG
                    plt.savefig(FN_PLOT_DOT, orientation='landscape', bbox_inches='tight', pad_inches=0.1)
                    print("[DONE] save as PNG: "+FN_PLOT_DOT)
        # ---
        
        # ---
        ## option 1e: dot plot by Cluster (w/o Dendrograms) ##
        if (EXPR_MK_PLT_1e==1):
            FN_PLOT_DOT=RESDIR_PRJ_EXPR+"1e_dot_exprMarker_"+PRJ+"_"+SMP+"_"+EXPR_MARKER_SET+"_nCat"+nMarkerCat_str+"_byClst"+"_"+EXPR_PLT+".pdf"
            
            isExist=os.path.exists(FN_PLOT_DOT)
            if not isExist:
                ## set plot dimension
                plt.rcParams["figure.figsize"] = [30, 8] #figure size in inches
                
                ## dot plot ##
                fig=sc.pl.dotplot(
                    adata_input, dict_marker_genes, 
                    groupby="seurat_clusters",
                    categories_order=ord_clst,
                    cmap=pal_expr_dot,
                    # title=mytitle,
                    swap_axes=False,
                    dendrogram=False
                    )
                
                ## save as PNG
                plt.savefig(FN_PLOT_DOT, orientation='landscape', bbox_inches='tight', pad_inches=0.1)
                print("[DONE] save as PNG: "+FN_PLOT_DOT)
        # ---
        
        # ---
        ## option 1f: dot plot by Cluster (w/ Dendrograms) ##
        if (EXPR_MK_PLT_1f==1):
            FN_PLOT_DOT=RESDIR_PRJ_EXPR+"1f_dot_exprMarker_"+PRJ+"_"+SMP+"_"+EXPR_MARKER_SET+"_nCat"+nMarkerCat_str+"_byClst"+"_wDendro"+"_"+EXPR_PLT+".pdf"
            
            isExist=os.path.exists(FN_PLOT_DOT)
            if not isExist:
                ## set plot dimension
                plt.rcParams["figure.figsize"] = [30, 8.5] #figure size in inches
                
                ## dot plot ##
                fig=sc.pl.dotplot(
                    adata_input, dict_marker_genes, 
                    groupby="seurat_clusters",
                    cmap=pal_expr_dot,
                    # title=mytitle,
                    swap_axes=False,
                    dendrogram=True
                    )
                
                ## save as PNG
                plt.savefig(FN_PLOT_DOT, orientation='landscape', bbox_inches='tight', pad_inches=0.1)
                print("[DONE] save as PNG: "+FN_PLOT_DOT)
        # ---

    # ---
    ## marker-based score ##
    ## - source: [https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.score_genes.html]
    
    # [scanpy.tl.score_genes]:
    # " Score a set of genes [Satija et al., 2015].

    # The score is the average expression of a set of genes subtracted with the average expression of a reference set of genes. 
    # The reference set is randomly sampled from the gene_pool for each binned expression value.

    # This reproduces the approach in Seurat [Satija et al., 2015] and has been implemented for Scanpy by Davide Cittaro. "

    if (EXPR_MK_SCORE==1):
        ## Make symmetric palette with vmin and vmax ##
        ## [CAR-T] ##
        ## - source: Figure 4 - Anderson ND et al. Nature Medicine 2023 [PMID:37407840]
        if (EXPR_MARKER_SET=="CAR-T"):
            list_mgset=[
                "CART_mem",
                "CART_eff",
                "CART_exh",
                "CART_res"                
            ]
        
        for mgset in list_mgset:
            ## [CAR-T] ##
            if (EXPR_MARKER_SET=="CAR-T"):
                if (mgset=="CART_mem"):
                    list_mg=["HNRNPLL","CCR7","SELL","LEF1","IL7R","TCF7"]
                    name_score_full="CAR-T: Stem-like memory"
                elif (mgset=="CART_eff"):
                    list_mg=["CD28","CD27","GZMK","GZMA","PRF1","CCL4","GZMM","CCL5","NKG7","GZMB","GNLY","LYAR","CXCR3","CXCR4","TXNIP"]
                    name_score_full="CAR-T: Activation, cytotoxicity, effector function"
                elif (mgset=="CART_exh"):
                    list_mg=["TIGIT","ENTPD1","HAVCR2","CTLA4","LAG3","PDCD1","EOMES","TOX","PMCH","BATF","PRDM1","CXCR1"]
                    name_score_full="CAR-T: Pre-exhaustion, exhaustion"
                elif (mgset=="CART_res"):
                    list_mg=["MCM5","TNFRSF9","ITGAE","CD69","NR4A2"]
                    name_score_full="CAR-T: Resident"
            
            ## count the number of genes in the selected marker gene set
            nMG=len(list_mg)
            nMG_str=str(nMG)
            ## define the name of the score
            name_score=mgset+nMG_str+"_score"
            print(name_score)
            
            ## make a customized color map
            ## - As with the convention in matplotlib, every continuous colormap has a reversed version, which has the suffix "_r": ["rocket"],["rocket_r"]
            pal_scr="RdBu_r" #high(dark_red) #low(dark_blue) [RdBu_r]
            
            ## plot 1: scoreUMAP (all clusters) ####
            if (SCR_PLT_1==1):
                ## re-set plot dimension
                plt.rcParams["figure.figsize"] = [5, 5] #figure size in inches
                
                FN_PLOT_UMAP=RESDIR_PRJ_SCR+"1_scoreUMAP_"+PRJ+"_"+SMP+"_"+EXPR_PLT+"_"+EXPR_MARKER_SET+"_"+mgset+"_nGene"+nMG_str+".png"
                
                isExist=os.path.exists(FN_PLOT_UMAP)
                if not isExist:
                    mytitle="[GCAR_PT01/"+SMP+"/"+EXPR_PLT+"] "+name_score+"\n(all clusters; nClst="+nClst_str+"; res="+res_clst_str+")"

                    ## Make mock column for plotting
                    sc.tl.score_genes(adata_input, list_mg, score_name=name_score)              
                    
                    ## To make a symmetric palette centerd around 0 we set vmax to maximal absolut value and vmin to the negative value of maxabs
                    maxabs = max(abs(adata_input.obs[name_score]))
                    size_dot=10
                    sc.pl.umap(
                        adata_input, color=name_score, cmap=pal_scr, s=size_dot, vmin=-maxabs, vmax=maxabs, title=mytitle                            
                    )
                    adata_input.obs.drop(name_score, axis=1, inplace=True)
                    
                    ## save as PNG
                    plt.savefig(FN_PLOT_UMAP, orientation='landscape', bbox_inches='tight', pad_inches=0.1)
                    print("[DONE] save as PNG: "+FN_PLOT_UMAP)
            
            
            ## plot 1a: scoreUMAP (by each cell cluster) ####                
            ## make a list of "string"-format Seurat Cluster variables
            if (EXPR_PLT=="VDJ"):

                if (SCR_PLT_1a==1):
                    list_clst=sorted(set(adata_input.obs.seurat_clusters))
                    print(str(list_clst)) #checkpoint
                    nCluster=len(list_clst)
                    
                    for c in list_clst:
                        ## add "zero" before the Seurat cluster number if the number is 1-digit
                        if (int(c) <10):
                            c_str="0"+str(c)
                        else:
                            c_str=str(c)
                        
                        FN_PLOT_UMAP=RESDIR_PRJ_SCR+"1a_scoreUMAP_"+PRJ+"_"+SMP+"_"+"_"+EXPR_PLT+"_"+EXPR_MARKER_SET+"_"+mgset+"_nGene"+nMG_str+"_nClst"+nClst_str+"-"+c_str+".png"
                        
                        isExist=os.path.exists(FN_PLOT_UMAP)
                        if not isExist:
                            ## Make mock column for plotting
                            sc.tl.score_genes(adata_input, list_mg, score_name=name_score)              
                            
                            ## To make a symmetric palette centerd around 0 we set vmax to maximal absolut value and vmin to the negative value of maxabs
                            maxabs = max(abs(adata_input.obs[name_score]))
                            
                            ## make 24 UMAP subplots (3x8) ##
                            ncol=8 #8
                            nrow=3 #1
                            figsize=5
                            wspace=0.25
                            
                            # Adapt figure size based on number of rows and columns and added space between them
                            # (e.g. wspace between columns)
                            fig,axs=plt.subplots(
                                nrows=nrow, ncols=ncol,
                                figsize=(ncol * figsize + (ncol - 1) * wspace * figsize, nrow * figsize)
                            )
                            plt.subplots_adjust(wspace=wspace)  

                            # ---
                            ## 3x8 subplots:
                            size_dot=15
                            
                            ## [all VDJ] ##
                            mytitle="["+EXPR_PLT+"] all"
                            umap0=sc.pl.umap(adata_input, color=name_score, cmap=pal_scr, s=size_dot, vmin=-maxabs, vmax=maxabs, mask_obs=(adata_input.obs.seurat_clusters==c), ax=axs[0,0], title=mytitle)
                            umap0_1=sc.pl.umap(adata_input[adata_input.obs.seurat_clusters==c], color=name_score, cmap=pal_scr, s=size_dot, vmin=-maxabs, vmax=maxabs, mask_obs=(adata_input[adata_input.obs.seurat_clusters==c].obs["sample_new"]=="Enriched"), ax=axs[0,1], title="Enriched")
                            umap0_2=sc.pl.umap(adata_input[adata_input.obs.seurat_clusters==c], color=name_score, cmap=pal_scr, s=size_dot, vmin=-maxabs, vmax=maxabs, mask_obs=(adata_input[adata_input.obs.seurat_clusters==c].obs["sample_new"]=="Harvest"), ax=axs[0,2], title="Harvest")
                            umap0_3=sc.pl.umap(adata_input[adata_input.obs.seurat_clusters==c], color=name_score, cmap=pal_scr, s=size_dot, vmin=-maxabs, vmax=maxabs, mask_obs=(adata_input[adata_input.obs.seurat_clusters==c].obs["sample_new"]=="D08"), ax=axs[0,3], title="D08")
                            umap0_4=sc.pl.umap(adata_input[adata_input.obs.seurat_clusters==c], color=name_score, cmap=pal_scr, s=size_dot, vmin=-maxabs, vmax=maxabs, mask_obs=(adata_input[adata_input.obs.seurat_clusters==c].obs["sample_new"]=="D15"), ax=axs[0,4], title="D15")
                            umap0_5=sc.pl.umap(adata_input[adata_input.obs.seurat_clusters==c], color=name_score, cmap=pal_scr, s=size_dot, vmin=-maxabs, vmax=maxabs, mask_obs=(adata_input[adata_input.obs.seurat_clusters==c].obs["sample_new"]=="D22"), ax=axs[0,5], title="D22")
                            umap0_6=sc.pl.umap(adata_input[adata_input.obs.seurat_clusters==c], color=name_score, cmap=pal_scr, s=size_dot, vmin=-maxabs, vmax=maxabs, mask_obs=(adata_input[adata_input.obs.seurat_clusters==c].obs["sample_new"]=="D28"), ax=axs[0,6], title="D28")
                            umap0_7=sc.pl.umap(adata_input[adata_input.obs.seurat_clusters==c], color=name_score, cmap=pal_scr, s=size_dot, vmin=-maxabs, vmax=maxabs, mask_obs=(adata_input[adata_input.obs.seurat_clusters==c].obs["sample_new"]=="D42"), ax=axs[0,7], title="D42")
                            
                            ## [VDJ: GCAR-pos] ##                                
                            mytitle="["+EXPR_PLT+"_"+"GCARpos"+"] all"
                            size_dot=20
                            umap1=sc.pl.umap(adata_input[adata_input.obs.expr_GCAR=="pos"], color=name_score, cmap=pal_scr, s=size_dot, vmin=-maxabs, vmax=maxabs ,mask_obs=(adata_input[adata_input.obs.expr_GCAR=="pos"].obs["seurat_clusters"]==c), ax=axs[1,0], title=mytitle)
                            umap1_1=sc.pl.umap(adata_input[(adata_input.obs.expr_GCAR=="pos") & (adata_input.obs.seurat_clusters==c)], color=name_score, cmap=pal_scr, s=size_dot, vmin=-maxabs, vmax=maxabs, mask_obs=(adata_input[(adata_input.obs.expr_GCAR=="pos") & (adata_input.obs.seurat_clusters==c)].obs["sample_new"]=="Enriched"), ax=axs[1,1], title="Enriched")
                            umap1_2=sc.pl.umap(adata_input[(adata_input.obs.expr_GCAR=="pos") & (adata_input.obs.seurat_clusters==c)], color=name_score, cmap=pal_scr, s=size_dot, vmin=-maxabs, vmax=maxabs, mask_obs=(adata_input[(adata_input.obs.expr_GCAR=="pos") & (adata_input.obs.seurat_clusters==c)].obs["sample_new"]=="Harvest"), ax=axs[1,2], title="Harvest")
                            umap1_3=sc.pl.umap(adata_input[(adata_input.obs.expr_GCAR=="pos") & (adata_input.obs.seurat_clusters==c)], color=name_score, cmap=pal_scr, s=size_dot, vmin=-maxabs, vmax=maxabs, mask_obs=(adata_input[(adata_input.obs.expr_GCAR=="pos") & (adata_input.obs.seurat_clusters==c)].obs["sample_new"]=="D08"), ax=axs[1,3], title="D08")
                            umap1_4=sc.pl.umap(adata_input[(adata_input.obs.expr_GCAR=="pos") & (adata_input.obs.seurat_clusters==c)], color=name_score, cmap=pal_scr, s=size_dot, vmin=-maxabs, vmax=maxabs, mask_obs=(adata_input[(adata_input.obs.expr_GCAR=="pos") & (adata_input.obs.seurat_clusters==c)].obs["sample_new"]=="D15"), ax=axs[1,4], title="D15")
                            umap1_5=sc.pl.umap(adata_input[(adata_input.obs.expr_GCAR=="pos") & (adata_input.obs.seurat_clusters==c)], color=name_score, cmap=pal_scr, s=size_dot, vmin=-maxabs, vmax=maxabs, mask_obs=(adata_input[(adata_input.obs.expr_GCAR=="pos") & (adata_input.obs.seurat_clusters==c)].obs["sample_new"]=="D22"), ax=axs[1,5], title="D22")
                            umap1_6=sc.pl.umap(adata_input[(adata_input.obs.expr_GCAR=="pos") & (adata_input.obs.seurat_clusters==c)], color=name_score, cmap=pal_scr, s=size_dot, vmin=-maxabs, vmax=maxabs, mask_obs=(adata_input[(adata_input.obs.expr_GCAR=="pos") & (adata_input.obs.seurat_clusters==c)].obs["sample_new"]=="D28"), ax=axs[1,6], title="D28")
                            umap1_7=sc.pl.umap(adata_input[(adata_input.obs.expr_GCAR=="pos") & (adata_input.obs.seurat_clusters==c)], color=name_score, cmap=pal_scr, s=size_dot, vmin=-maxabs, vmax=maxabs, mask_obs=(adata_input[(adata_input.obs.expr_GCAR=="pos") & (adata_input.obs.seurat_clusters==c)].obs["sample_new"]=="D42"), ax=axs[1,7], title="D42")

                            ## [VDJ: GCAR-neg] ##
                            mytitle="["+EXPR_PLT+"_"+"GCARneg"+"] all"
                            size_dot=20
                            umap2=sc.pl.umap(adata_input[adata_input.obs.expr_GCAR=="neg"], color=name_score, cmap=pal_scr, s=size_dot, vmin=-maxabs, vmax=maxabs ,mask_obs=(adata_input[adata_input.obs.expr_GCAR=="neg"].obs["seurat_clusters"]==c), ax=axs[2,0], title=mytitle)
                            umap2_1=sc.pl.umap(adata_input[(adata_input.obs.expr_GCAR=="neg") & (adata_input.obs.seurat_clusters==c)], color=name_score, cmap=pal_scr, s=size_dot, vmin=-maxabs, vmax=maxabs, mask_obs=(adata_input[(adata_input.obs.expr_GCAR=="neg") & (adata_input.obs.seurat_clusters==c)].obs["sample_new"]=="Enriched"), ax=axs[2,1], title="Enriched")
                            umap2_2=sc.pl.umap(adata_input[(adata_input.obs.expr_GCAR=="neg") & (adata_input.obs.seurat_clusters==c)], color=name_score, cmap=pal_scr, s=size_dot, vmin=-maxabs, vmax=maxabs, mask_obs=(adata_input[(adata_input.obs.expr_GCAR=="neg") & (adata_input.obs.seurat_clusters==c)].obs["sample_new"]=="Harvest"), ax=axs[2,2], title="Harvest")
                            umap2_3=sc.pl.umap(adata_input[(adata_input.obs.expr_GCAR=="neg") & (adata_input.obs.seurat_clusters==c)], color=name_score, cmap=pal_scr, s=size_dot, vmin=-maxabs, vmax=maxabs, mask_obs=(adata_input[(adata_input.obs.expr_GCAR=="neg") & (adata_input.obs.seurat_clusters==c)].obs["sample_new"]=="D08"), ax=axs[2,3], title="D08")
                            umap2_4=sc.pl.umap(adata_input[(adata_input.obs.expr_GCAR=="neg") & (adata_input.obs.seurat_clusters==c)], color=name_score, cmap=pal_scr, s=size_dot, vmin=-maxabs, vmax=maxabs, mask_obs=(adata_input[(adata_input.obs.expr_GCAR=="neg") & (adata_input.obs.seurat_clusters==c)].obs["sample_new"]=="D15"), ax=axs[2,4], title="D15")
                            umap2_5=sc.pl.umap(adata_input[(adata_input.obs.expr_GCAR=="neg") & (adata_input.obs.seurat_clusters==c)], color=name_score, cmap=pal_scr, s=size_dot, vmin=-maxabs, vmax=maxabs, mask_obs=(adata_input[(adata_input.obs.expr_GCAR=="neg") & (adata_input.obs.seurat_clusters==c)].obs["sample_new"]=="D22"), ax=axs[2,5], title="D22")
                            umap2_6=sc.pl.umap(adata_input[(adata_input.obs.expr_GCAR=="neg") & (adata_input.obs.seurat_clusters==c)], color=name_score, cmap=pal_scr, s=size_dot, vmin=-maxabs, vmax=maxabs, mask_obs=(adata_input[(adata_input.obs.expr_GCAR=="neg") & (adata_input.obs.seurat_clusters==c)].obs["sample_new"]=="D28"), ax=axs[2,6], title="D28")
                            umap2_7=sc.pl.umap(adata_input[(adata_input.obs.expr_GCAR=="neg") & (adata_input.obs.seurat_clusters==c)], color=name_score, cmap=pal_scr, s=size_dot, vmin=-maxabs, vmax=maxabs, mask_obs=(adata_input[(adata_input.obs.expr_GCAR=="neg") & (adata_input.obs.seurat_clusters==c)].obs["sample_new"]=="D42"), ax=axs[2,7], title="D42")
                            
                            adata_input.obs.drop(name_score, axis=1, inplace=True)

                            ## plot title ##
                            mysuptitle="[GCAR_PT01/"+SMP+"/"+EXPR_PLT+"] - {"+name_score_full+"}"+" (nGene:"+nMG_str+") - Cluster #"+c_str+" (nClst="+nClst_str+"; res="+res_clst_str+")"
                            fig.suptitle(mysuptitle, fontsize=16, y=0.93) #y=0.98(default) #3x8 subplots

                            
                            ## save as PNG
                            plt.savefig(FN_PLOT_UMAP, orientation='landscape', bbox_inches='tight', pad_inches=0.1)
                            print("[DONE] save as PNG: "+FN_PLOT_UMAP)
                            # ---
# ===

# ===
## gene expression heatmap - using Scanpy's built-in function ##
## - source: [https://scanpy.readthedocs.io/en/stable/generated/scanpy.pl.heatmap.html]
if (EXPR_MK_HMP_SC==1):
    ## import relevant modules
    # pip install marsilea #need-to-be-run-on-CLI
    # import marsilea as ma
    # import marsilea.plotter as mp
    
    if (TCELL_ONLY!=1):
        PATH_DATA=DATADIR+PRJ+"_"+SMP+"_afterUMAP"+"_raw"+"_res"+res_clst_str+".h5ad"
    elif (TCELL_ONLY==1):
        PATH_DATA=DATADIR+PRJ+"_"+SMP+"_afterUMAP"+"_raw_vdj"+"_res"+res_clst_str+".h5ad"
    
    ## step 0 | load the [adata] object ##
    PATH_DATA=DATADIR+"GCAR1_PT01_mrg7_afterUMAP_raw_vdj_res2.4.h5ad"
    adata=sc.read(PATH_DATA)
    print(adata.shape) #checkpoint
    # (77651, 32510)
    print(pd.DataFrame(adata.obs).head()) #checkpoint
    print(pd.DataFrame(adata.obs).tail()) #checkpoint

    ## set the Seurat clustering resolution
    res_clst=2.4
    res_clst_str=str(res_clst)
    
    ## count the number of the [Seurat Cell Clusters]
    list_clst=set(adata.obs["seurat_clusters"])
    print(list_clst)
    nClst=len(list_clst)
    nClst_str=str(nClst)

    
    # ---
    ## set the input [adata]
    ## step 1 | select cells that meet TRB-based filtering criteria ##
    adata_input=adata[adata.obs["select_combLen_chain1_chain2_gene"]=="Y_select"]
    print(adata_input.shape) #checkpoint
    # (68862, 32510)
    print(pd.DataFrame(adata_input.obs).head()) #checkpoint
    print(pd.DataFrame(adata_input.obs).tail()) #checkpoint
    # ---
    
    
    ## step 1a | add a column based on mulitple conditions ##
        ## create a list of conditions
    conditions=[
        ## pos ##
        (adata_input.obs["orig.ident"]=="Enriched_Apheresis") & (adata_input.obs["expr_GCAR"]=="pos"),
        (adata_input.obs["orig.ident"]=="GCAR1_Product") & (adata_input.obs["expr_GCAR"]=="pos"),
        (adata_input.obs["orig.ident"]=="D08") & (adata_input.obs["expr_GCAR"]=="pos"),
        (adata_input.obs["orig.ident"]=="D15") & (adata_input.obs["expr_GCAR"]=="pos"),
        (adata_input.obs["orig.ident"]=="D22") & (adata_input.obs["expr_GCAR"]=="pos"),
        (adata_input.obs["orig.ident"]=="D28") & (adata_input.obs["expr_GCAR"]=="pos"),
        (adata_input.obs["orig.ident"]=="D42") & (adata_input.obs["expr_GCAR"]=="pos"),
        ## neg ##
        (adata_input.obs["orig.ident"]=="Enriched_Apheresis") & (adata_input.obs["expr_GCAR"]=="neg"),
        (adata_input.obs["orig.ident"]=="GCAR1_Product") & (adata_input.obs["expr_GCAR"]=="neg"),
        (adata_input.obs["orig.ident"]=="D08") & (adata_input.obs["expr_GCAR"]=="neg"),
        (adata_input.obs["orig.ident"]=="D15") & (adata_input.obs["expr_GCAR"]=="neg"),
        (adata_input.obs["orig.ident"]=="D22") & (adata_input.obs["expr_GCAR"]=="neg"),
        (adata_input.obs["orig.ident"]=="D28") & (adata_input.obs["expr_GCAR"]=="neg"),
        (adata_input.obs["orig.ident"]=="D42") & (adata_input.obs["expr_GCAR"]=="neg")
    ]
    ## create a list of the "preferred" values to assign for each condition
    values=[
        "Enriched_pos","Harvest_pos","D08_pos","D15_pos","D22_pos","D28_pos","D42_pos",
        "Enriched_neg","Harvest_neg","D08_neg","D15_neg","D22_neg","D28_neg","D42_neg"
        ]
    
    ## create a new column to assign values
    adata_input.obs["sample_new__expr_GCAR"] = np.select(conditions,values)
    print(set(adata_input.obs["sample_new__expr_GCAR"]))
    print(pd.DataFrame(adata_input.obs).head()) #checkpoint
    print(pd.DataFrame(adata_input.obs).tail()) #checkpoint
    # ---
    
    
    ## step 2 | make the heatmaps for each cell cluster ##
    ## example cell clusters:
    ##   - [cluster #3] CD15-high, CD8em_td = [ CD8em_td__D15_pos ]
    ##   - [cluster #10] CD15-high, CD8em = [ CD8em__D15_pos ]
    for coi in list_clst:
        ## set a "cluster-of-interest" (i.e., [coi])
        coi=int(coi) #3 #10
        coi_str=str(coi)
        print("--- cluster-of-interest: [ "+coi_str+" ] ---")
        
        if (coi <10):
            coi_str_lab="0"+str(coi)
        else:
            coi_str_lab=str(coi)
    
        adata_input_coi=adata_input[adata_input.obs["RNA_snn_res.2.4"]==coi_str]
        print('[adata_input_coi]:')
        print(adata_input_coi.shape) #checkpoint
        print(pd.DataFrame(adata_input_coi.obs).head()) #checkpoint
        print(pd.DataFrame(adata_input_coi.obs).tail()) #checkpoint
        # ---

        ## set the "markers-of-interest"
        ## - source: [CAR-T immunophenotyping marker list]
        ## Activation, cytotoxicity, effector function ##
        markers = ["CD28","CD27","GZMK","GZMA","PRF1","CCL4","GZMM","CCL5","NKG7","GZMB","GNLY","LYAR","CXCR3","CXCR4","TXNIP"]

        # ---
        ## plot 1 | heatmap using GEX data - for each cell cluster (w/o dendrogram) ##
        FN_PLOT_HMAP=RESDIR_PRJ+"1_HMAP_expr_raw_"+PRJ+"_"+SMP+"_byGCARexpr"+"_res"+res_clst_str+"_nClst"+nClst_str+"_clst"+coi_str_lab+"_woDendro"+".png"
        isExist=os.path.exists(FN_PLOT_HMAP)
        if not isExist:  
            sc.pl.heatmap(adata_input_coi, markers, groupby='sample_new__expr_GCAR', swap_axes=True)
            plt.savefig(FN_PLOT_HMAP, orientation='landscape', bbox_inches='tight', pad_inches=0.1)
            if (os.path.exists(FN_PLOT_HMAP)):
                print("[DONE] save as PNG: "+FN_PLOT_HMAP)

        ## plot 1a | heatmap using GEX data - for each cell cluster (w/ dendrograms) ##
        FN_PLOT_HMAP=RESDIR_PRJ+"1a_HMAP_expr_raw_"+PRJ+"_"+SMP+"_byGCARexpr"+"_res"+res_clst_str+"_nClst"+nClst_str+"_clst"+coi_str_lab+"_wDendro"+".png"
        isExist=os.path.exists(FN_PLOT_HMAP)
        if not isExist: 
            sc.pl.heatmap(adata_input_coi, markers, groupby='sample_new__expr_GCAR', swap_axes=True, dendrogram=True)
            plt.savefig(FN_PLOT_HMAP, orientation='landscape', bbox_inches='tight', pad_inches=0.1)
            if (os.path.exists(FN_PLOT_HMAP)):
                print("[DONE] save as PNG: "+FN_PLOT_HMAP)
        # ---
# ===


# ===
## cluster-level annotation ##
if (VIS_CLST==1):
    ## load the h5ad file
    if (PRJ=="GCAR1_PT01"):
        if (TCELL_ONLY!=1):
            PATH_DATA=DATADIR+PRJ+"_"+SMP+"_afterUMAP"+"_res"+res_clst_str+".h5ad"
        else:
            PATH_DATA=DATADIR+PRJ+"_"+SMP+"_afterUMAP"+"_SCim_constrain_VDJonly_res"+res_clst_str+".h5ad"
        if os.path.exists(PATH_DATA):
            print("[START] load the h5ad file: "+PATH_DATA)
        else:
            print("[CHECK] h5ad file does not exist: "+PATH_DATA)
    
    adata=sc.read(PATH_DATA)
    print(adata.shape) #checkpoint
    print("[END] load the h5ad file: "+PATH_DATA)

    # ---
    ## add a column to [obs]: [sample_new] ##
    ## create a list of conditions
    conditions=[
        (adata.obs["orig.ident"]=="Enriched_Apheresis"),
        (adata.obs["orig.ident"]=="GCAR1_Product"),
        (adata.obs["orig.ident"]=="D08"),
        (adata.obs["orig.ident"]=="D15"),
        (adata.obs["orig.ident"]=="D22"),
        (adata.obs["orig.ident"]=="D28"),
        (adata.obs["orig.ident"]=="D42")
    ]
    ## create a list of the "preferred" values to assign for each condition
    values=["Enriched","Harvest","D08","D15","D22","D28","D42"]
    
    ## create a new column to assign values
    adata.obs["sample_new"] = np.select(conditions,values)
    print(set(adata.obs["sample_new"]))
    print(pd.DataFrame(adata.obs).head()) #checkpoint
    print(pd.DataFrame(adata.obs).tail()) #checkpoint
    # ---
    
    # ---
    ## add a column to [obs]: [cluster_celltype] ##        
    ## create a list of conditions
    conditions=[
        (adata.obs["seurat_clusters"]=="0"),
        (adata.obs["seurat_clusters"]=="1"),
        (adata.obs["seurat_clusters"]=="2"),
        (adata.obs["seurat_clusters"]=="3"),
        (adata.obs["seurat_clusters"]=="4"),
        (adata.obs["seurat_clusters"]=="5"),
        (adata.obs["seurat_clusters"]=="6"),
        (adata.obs["seurat_clusters"]=="7"),
        (adata.obs["seurat_clusters"]=="8"),
        (adata.obs["seurat_clusters"]=="9"),
        
        (adata.obs["seurat_clusters"]=="10"),
        (adata.obs["seurat_clusters"]=="11"),
        (adata.obs["seurat_clusters"]=="12"),
        (adata.obs["seurat_clusters"]=="13"),
        (adata.obs["seurat_clusters"]=="14"),
        (adata.obs["seurat_clusters"]=="15"),
        (adata.obs["seurat_clusters"]=="16"),
        (adata.obs["seurat_clusters"]=="17"),
        (adata.obs["seurat_clusters"]=="18"),
        (adata.obs["seurat_clusters"]=="19"),
        
        (adata.obs["seurat_clusters"]=="20"),
        (adata.obs["seurat_clusters"]=="21"),
        (adata.obs["seurat_clusters"]=="22"),
        (adata.obs["seurat_clusters"]=="23"),
        (adata.obs["seurat_clusters"]=="24"),
        (adata.obs["seurat_clusters"]=="25"),
        (adata.obs["seurat_clusters"]=="26"),
        (adata.obs["seurat_clusters"]=="27"),
        (adata.obs["seurat_clusters"]=="28"),
        (adata.obs["seurat_clusters"]=="29"),

        (adata.obs["seurat_clusters"]=="30"),
        (adata.obs["seurat_clusters"]=="31"),
        (adata.obs["seurat_clusters"]=="32"),
        (adata.obs["seurat_clusters"]=="33"),
        (adata.obs["seurat_clusters"]=="34"),
        (adata.obs["seurat_clusters"]=="35"),
        (adata.obs["seurat_clusters"]=="36"),
        (adata.obs["seurat_clusters"]=="37"),
        (adata.obs["seurat_clusters"]=="38"),
        (adata.obs["seurat_clusters"]=="39"),

        (adata.obs["seurat_clusters"]=="40"),
        (adata.obs["seurat_clusters"]=="41"),
        (adata.obs["seurat_clusters"]=="42"),
        (adata.obs["seurat_clusters"]=="43"),
        (adata.obs["seurat_clusters"]=="44"),
        (adata.obs["seurat_clusters"]=="45"),
        (adata.obs["seurat_clusters"]=="46")
        ]

    ## create a list of the "preferred" values to assign for each condition
    if (VIS_CLST_CCLASS==1):
        values=[
            "CD8", #0       
            "CD8", #1       
            "CD4", #2       
            "CD8", #3       
            "CD4", #4       
            "CD4", #5       
            "CD4", #6       
            "CD4", #7       
            "Monocyte", #8  
            "CD4", #9       
            
            "CD8", #10      
            "MAIT", #11     
            "Monocyte", #12 
            "Monocyte", #13 
            "CD8", #14      
            "CD4", #15      
            "CD4", #16      
            "CD8", #17      
            "NKT", #18      
            "CD8", #19      
            
            "Monocyte", #20 
            "CD8", #21      
            "Monocyte", #22 
            "CD8", #23      
            "Monocyte", #24 
            "NK", #25
            "Monocyte", #26 
            "Monocyte", #27 
            "CD4", #28      
            "mix_CD14monocyte_CD4naive", #29
            
            "CD8", #30      
            "Lymphocyte", #31
            "CD4", #32      
            "CD4", #33      
            "Lymphocyte", #34
            "Lymphocyte", #35
            "CD8", #36      
            "Monocyte", #37 
            "CD4", #38      
            "DCs", #39      
            
            "B-cells", #40  
            "Monocyte", #41 
            "CD4", #42      
            "CD4", #43      
            "CD8", #44      
            "HSC", #45      
            "CD8" #46
            ]
    
    elif (VIS_CLST_TPOINT==1):
        values=[
            "circulating", #0       
            "circulating", #1       
            "circulating", #2       
            "D15", #3       
            "Harvest", #4   
            "circulating", #5       
            "circulating", #6       
            "circulating", #7       
            "circulating", #8       
            "Enriched", #9  
            
            "D15", #10      
            "circulating", #11      
            "circulating", #12      
            "circulating", #13      
            "circulating", #14      
            "circulating", #15      
            "Enriched", #16 
            "D22", #17      
            "circulating", #18      
            "D15-28", #19   
            
            "circulating", #20      
            "Enriched", #21 
            "D22", #22      
            "circulating", #23      
            "D22", #24      
            "circulating", #25      
            "circulating", #26      
            "circulating", #27      
            "Harvest", #28  
            "D22_D42", #29  
            
            "D15", #30      
            "circulating", #31      
            "circulating", #32      
            "circulating", #33      
            "D15", #34      
            "D15", #35      
            "Enriched", #36 
            "D15", #37      
            "Harvest", #38  
            "circulating", #39      
            
            "circulating", #40      
            "D42", #41      
            "Harvest", #42  
            "Enriched", #43 
            "D15", #44      
            "circulating", #45      
            "D15" #46
            ]
        
    elif (VIS_CLST_CTYPE==1):
        values=[
            "CD8__circulating__neg", #0
            "CD8__circulating__neg", #1
            "CD4__circulating__neg", #2
            "CD8__D15__pos", #3
            "CD4__Harvest__pos", #4
            "CD4__circulating__neg", #5
            "CD4__circulating__neg", #6
            "CD4__circulating__neg", #7
            "Monocyte__circulating__neg", #8
            "CD4__Enriched__neg", #9
            
            "CD8__D15__pos", #10
            "MAIT__circulating__neg", #11
            "Monocyte__circulating__neg", #12
            "Monocyte__circulating__neg", #13
            "CD8__circulating__neg", #14
            "CD4__circulating__neg", #15
            "CD4__Enriched__neg", #16
            "CD8__D22__pos", #17
            "NKT__circulating__neg", #18
            "CD8__D15-28__pos", #19
            
            "Monocyte__circulating__neg", #20
            "CD8__Enriched__neg", #21
            "Monocyte__D22__neg", #22
            "CD8__circulating__neg", #23
            "Monocyte__D22__neg", #24
            "NK__circulating__neg", #25
            "Monocyte__circulating__neg", #26
            "Monocyte__circulating__neg", #27
            "CD4__Harvest__pos", #28
            "mix_CD14monocyte_CD4naive__D22_D42__neg", #29 #excluded
            
            "CD8__D15__pos", #30
            "Lymphocyte__circulating__neg", #31
            "CD4__circulating__neg", #32
            "CD4__circulating__neg", #33
            "Lymphocyte__D15__pos", #34
            "Lymphocyte__D15__pos", #35
            "CD8__Enriched__neg", #36
            "Monocyte__D15__pos", #37
            "CD4__Harvest__pos", #38
            "DCs__circulating__neg", #39
            
            "B-cells__circulating__neg", #40
            "Monocyte__D42__neg", #41
            "CD4__Harvest__pos", #42
            "CD4__Enriched__neg", #43 #excluded
            "CD8__D15__pos", #44
            "HSC__circulating__neg", #45
            "CD8__D15__pos" #46
            ]
    
    ## create a new column to assign values
    print("- length of [conditions]: "+str(len(conditions)))
    print("- length of [values]: "+str(len(values)))
    adata.obs["cluster_celltype"] = np.select(conditions,values)
    print(set(adata.obs["cluster_celltype"]))
    print(pd.DataFrame(adata.obs).head()) #checkpoint
    print(pd.DataFrame(adata.obs).tail()) #checkpoint
    # ---
    
    ## count the number of the Seurat Clusters
    list_clst=set(adata.obs["seurat_clusters"])
    nClst=len(list_clst)
    nClst_str=str(nClst)
    ## count the number of the categorical variables
    nCtgr=len(set(adata.obs["cluster_celltype"]))
    nCtgr_str=str(nCtgr)
    
    ## colour palettes ##
    pal_name="Paired" #"tab10" "Paired"
    colors_clst=glasbey.extend_palette(pal_name, palette_size=nClst, lightness_bounds=(30, 60), colorblind_safe=True, cvd_severity=100) #default:lightness_bounds=(10,90) #lightness_bounds=(30, 60) -- from the manual page
    pal_clst=sns.set_palette(sns.color_palette(colors_clst))

    # ---
    ## plot 1: UMAP at cluster-level ##
    if (CLST_MK_PLT_1==1):
        if (SMP=="mrg" or SMP=="mrg7"):
            if (VIS_CLST_CCLASS==1):
                FN_PLOT_UMAP=RESDIR_PRJ_CLST+"1_UMAP_clusterLevel_cellClass_"+PRJ+"_"+SMP+"_res"+res_clst_str+"_nClst"+nClst_str+"_nCClass"+nCtgr_str+".pdf" #changed to "pdf" due to low-resolution of the legends in the "png" format                
                mytitle="Cluster-level annotation (nCluster="+nClst_str+")\n(nCClass="+nCtgr_str+")"
            elif (VIS_CLST_TPOINT==1):
                FN_PLOT_UMAP=RESDIR_PRJ_CLST+"1_UMAP_clusterLevel_timePoint_"+PRJ+"_"+SMP+"_res"+res_clst_str+"_nClst"+nClst_str+"_nTPoint"+nCtgr_str+".pdf" #changed to "pdf" due to low-resolution of the legends in the "png" format                
                mytitle="Cluster-level annotation (nCluster="+nClst_str+")\n(nTimepoint="+nCtgr_str+")"
            elif (VIS_CLST_CTYPE==1):
                FN_PLOT_UMAP=RESDIR_PRJ_CLST+"1_UMAP_clusterLevel_cellType_"+PRJ+"_"+SMP+"_res"+res_clst_str+"_nClst"+nClst_str+"_nCType"+nCtgr_str+".pdf" #changed to "pdf" due to low-resolution of the legends in the "png" format                
                mytitle="Cluster-level annotation (nCluster="+nClst_str+")\n(nCellType="+nCtgr_str+")"

            size_dot=2.5

            isExist=os.path.exists(FN_PLOT_UMAP)
            if not isExist:
                # mytitle="Cluster-level annotation (nCluster="+nClst_str+")\n nCtgr="+nCtgr_str+")"
                size_dot=5
                umap=sc.pl.umap(
                    adata, 
                    color="cluster_celltype", 
                    palette=pal_clst,
                    legend_fontsize=5,
                    s=size_dot,
                    title=mytitle
                    )

                plt.savefig(FN_PLOT_UMAP, orientation='landscape', bbox_inches='tight', pad_inches=0.1)
                if (os.path.exists(FN_PLOT_UMAP)):
                    print("[DONE] save as PDF: "+FN_PLOT_UMAP)
    # ---
    
    # ---
    ## plot 1a: UMAP at cluster-level ##
    if (CLST_MK_PLT_1a==1):
        
        
        if (SMP=="mrg" or SMP=="mrg7"):
            if (VIS_CLST_CCLASS==1):
                FN_PLOT_UMAP=RESDIR_PRJ_CLST+"1a_UMAP_clusterLevel_cellClass_"+PRJ+"_"+SMP+"_res"+res_clst_str+"_nClst"+nClst_str+"_nCClass"+nCtgr_str+".pdf" #changed to "pdf" due to low-resolution of the legends in the "png" format                
                mytitle="Cluster-level annotation (nCluster="+nClst_str+")\n(nCClass="+nCtgr_str+")"
            elif (VIS_CLST_TPOINT==1):
                FN_PLOT_UMAP=RESDIR_PRJ_CLST+"1a_UMAP_clusterLevel_timePoint_"+PRJ+"_"+SMP+"_res"+res_clst_str+"_nClst"+nClst_str+"_nTPoint"+nCtgr_str+".pdf" #changed to "pdf" due to low-resolution of the legends in the "png" format                
                mytitle="Cluster-level annotation (nCluster="+nClst_str+")\n(nTimepoint="+nCtgr_str+")"
            elif (VIS_CLST_CTYPE==1):
                FN_PLOT_UMAP=RESDIR_PRJ_CLST+"1a_UMAP_clusterLevel_cellType_"+PRJ+"_"+SMP+"_res"+res_clst_str+"_nClst"+nClst_str+"_nCType"+nCtgr_str+".pdf" #changed to "pdf" due to low-resolution of the legends in the "png" format                
                mytitle="Cluster-level annotation (nCluster="+nClst_str+")\n(nCellType="+nCtgr_str+")"

            size_dot=2.5

            isExist=os.path.exists(FN_PLOT_UMAP)
            if not isExist:
                size_dot=5
                umap=sc.pl.umap(
                    adata, 
                    color="cluster_celltype", 
                    palette=pal_clst,
                    legend_fontsize=5,
                    s=size_dot,
                    title=mytitle
                    )

                plt.savefig(FN_PLOT_UMAP, orientation='landscape', bbox_inches='tight', pad_inches=0.1)
                if (os.path.exists(FN_PLOT_UMAP)):
                    print("[DONE] save as PDF: "+FN_PLOT_UMAP)
    # ---
    
    
    # ---
    ## plot 2: UMAP at cluster-level -- by sample/GCAR-pos/GCAR-neg ##
    if (CLST_MK_PLT_2==1):
        if (SMP=="mrg" or SMP=="mrg7"):
            
            ## horizontal version (1 x 4) ##
            FN_PLOT_UMAP=RESDIR_PRJ_CLST+"2_UMAP_clusterLevel_cellType_"+PRJ+"_"+SMP+"_res"+res_clst_str+"_bySmp_exprGCAR"+".png"
            isExist=os.path.exists(FN_PLOT_UMAP)
            if not isExist:

                ## make 4 UMAP subplots - w/ different colour palette and/or cmap parameters ##
                ncol=4
                nrow=1
                figsize=5
                wspace=0.3
                
                # Adapt figure size based on number of rows and columns and added space between them
                # (e.g. wspace between columns)
                fig,axs=plt.subplots(
                    nrows=nrow, ncols=ncol,
                    figsize=(ncol * figsize + (ncol - 1) * wspace * figsize, nrow * figsize)
                )
                plt.subplots_adjust(wspace=wspace)

                # note: ValueError: Groups and mask arguments are incompatible.
                mytitle="all (nCellType="+nCtgr_str+")"
                size_dot=5
                umap0=sc.pl.umap(adata, color="cluster_celltype", palette=pal_clst, show=False, legend_fontsize=5, s=size_dot, ax=axs[3], title=mytitle) #position:right(due to the length of the cell types)
                umap1=sc.pl.umap(adata, mask_obs=(adata.obs["cell_VDJ"]=="VDJ"), color="sample_new", palette=pal_smp, show=False, legend_fontsize=5, s=size_dot, ax=axs[0], title="by sample (T cells)")
                umap2=sc.pl.umap(adata[adata.obs.cell_VDJ=="VDJ"], color="expr_GCAR", groups=["pos"], show=False, palette=pal_expr_cat_pos, legend_fontsize=5, s=size_dot, ax=axs[1], title="T cells (GCAR-pos)")
                umap3=sc.pl.umap(adata[adata.obs.cell_VDJ=="VDJ"], color="expr_GCAR", groups=["neg"], show=False, palette=pal_expr_cat_neg, legend_fontsize=5, s=size_dot, ax=axs[2], title="T cells (GCAR-neg)")
                
                ## change "NA" values in the UMAP plot to a specific value
                legend_texts = umap1.get_legend().get_texts()
                for legend_text in legend_texts:
                    if legend_text.get_text() == "NA":
                        legend_text.set_text("non-T cells")
                
                ## change "NA" values in the UMAP plot to a specific value
                legend_texts = umap2.get_legend().get_texts()
                for legend_text in legend_texts:
                    if legend_text.get_text() == "NA":
                        legend_text.set_text("neg")
                        
                ## change "NA" values in the UMAP plot to a specific value
                legend_texts = umap3.get_legend().get_texts()
                for legend_text in legend_texts:
                    if legend_text.get_text() == "NA":
                        legend_text.set_text("pos")
                
                plt.savefig(FN_PLOT_UMAP, orientation='landscape', bbox_inches='tight', pad_inches=0.1)
                if (os.path.exists(FN_PLOT_UMAP)):
                    print("[DONE] save as PDF: "+FN_PLOT_UMAP)
    # ---
    
    # ---
    ## plot 3: UMAP at cluster-level -- by sample/GCAR-pos/GCAR-neg ##
    if (CLST_MK_PLT_3==1):
        if (SMP=="mrg" or SMP=="mrg7"):
            
            ## horizontal version (8 x 4) ##
            FN_PLOT_UMAP=RESDIR_PRJ_CLST+"3_UMAP_clusterLevel_cellType_"+PRJ+"_"+SMP+"_res"+res_clst_str+"_bySmp_exprGCAR"+"_4x8"+".png"
            isExist=os.path.exists(FN_PLOT_UMAP)
            if not isExist:

                ## make 32 UMAP subplots - w/ different colour palette and/or cmap parameters ##
                ncol=8
                nrow=4
                figsize=5
                wspace=0.25
                
                # Adapt figure size based on number of rows and columns and added space between them
                # (e.g. wspace between columns)
                fig,axs=plt.subplots(
                    nrows=nrow, ncols=ncol,
                    figsize=(ncol * figsize + (ncol - 1) * wspace * figsize, nrow * figsize)
                )
                plt.subplots_adjust(wspace=wspace)

                # note: ValueError: Groups and mask arguments are incompatible.
                ## cluster-level cell type annotation ##
                mytitle="all (nCellType="+nCtgr_str+")"
                size_dot=5
                umap0=sc.pl.umap(adata, color="cluster_celltype", palette=pal_clst, show=False, legend_fontsize=5, s=size_dot, ax=axs[0,0], title=mytitle, legend_loc=None) #position:right(due to the length of the cell types)
                umap0_1=sc.pl.umap(adata[adata.obs.seurat_clusters.isin(list_clst)], mask_obs=(adata[adata.obs.seurat_clusters.isin(list_clst)].obs["sample_new"]=="Enriched"), color="cluster_celltype", palette=pal_clst, show=False, legend_fontsize=5, s=size_dot, ax=axs[0,1], title="Enriched", legend_loc=None)
                umap0_2=sc.pl.umap(adata[adata.obs.seurat_clusters.isin(list_clst)], mask_obs=(adata[adata.obs.seurat_clusters.isin(list_clst)].obs["sample_new"]=="Harvest"), color="cluster_celltype", palette=pal_clst, show=False, legend_fontsize=5, s=size_dot, ax=axs[0,2], title="Harvest", legend_loc=None)
                umap0_3=sc.pl.umap(adata[adata.obs.seurat_clusters.isin(list_clst)], mask_obs=(adata[adata.obs.seurat_clusters.isin(list_clst)].obs["sample_new"]=="D08"), color="cluster_celltype", palette=pal_clst, show=False, legend_fontsize=5, s=size_dot, ax=axs[0,3], title="D08", legend_loc=None)
                umap0_4=sc.pl.umap(adata[adata.obs.seurat_clusters.isin(list_clst)], mask_obs=(adata[adata.obs.seurat_clusters.isin(list_clst)].obs["sample_new"]=="D15"), color="cluster_celltype", palette=pal_clst, show=False, legend_fontsize=5, s=size_dot, ax=axs[0,4], title="D15", legend_loc=None)
                umap0_5=sc.pl.umap(adata[adata.obs.seurat_clusters.isin(list_clst)], mask_obs=(adata[adata.obs.seurat_clusters.isin(list_clst)].obs["sample_new"]=="D22"), color="cluster_celltype", palette=pal_clst, show=False, legend_fontsize=5, s=size_dot, ax=axs[0,5], title="D22", legend_loc=None)
                umap0_6=sc.pl.umap(adata[adata.obs.seurat_clusters.isin(list_clst)], mask_obs=(adata[adata.obs.seurat_clusters.isin(list_clst)].obs["sample_new"]=="D28"), color="cluster_celltype", palette=pal_clst, show=False, legend_fontsize=5, s=size_dot, ax=axs[0,6], title="D28", legend_loc=None)
                umap0_7=sc.pl.umap(adata[adata.obs.seurat_clusters.isin(list_clst)], mask_obs=(adata[adata.obs.seurat_clusters.isin(list_clst)].obs["sample_new"]=="D42"), color="cluster_celltype", palette=pal_clst, show=False, legend_fontsize=5, s=size_dot, ax=axs[0,7], title="D42")
                                
                umap1=sc.pl.umap(adata, mask_obs=(adata.obs["cell_VDJ"]=="VDJ"), color="sample_new", palette=pal_smp, show=False, legend_fontsize=5, s=size_dot, ax=axs[1,0], title="by sample (T cells)", legend_loc=None)
                umap1_1=sc.pl.umap(adata[adata.obs.cell_VDJ=="VDJ"], mask_obs=(adata[adata.obs.cell_VDJ=="VDJ"].obs["sample_new_VDJ"]=="Enriched_VDJ"), color="sample_new", palette=pal_smp, show=False, legend_fontsize=5, s=size_dot, ax=axs[1,1], title="Enriched (T cells)", legend_loc=None)
                umap1_2=sc.pl.umap(adata[adata.obs.cell_VDJ=="VDJ"], mask_obs=(adata[adata.obs.cell_VDJ=="VDJ"].obs["sample_new_VDJ"]=="Harvest_VDJ"), color="sample_new", palette=pal_smp, show=False, legend_fontsize=5, s=size_dot, ax=axs[1,2], title="Harvest (T cells)", legend_loc=None)
                umap1_3=sc.pl.umap(adata[adata.obs.cell_VDJ=="VDJ"], mask_obs=(adata[adata.obs.cell_VDJ=="VDJ"].obs["sample_new_VDJ"]=="D08_VDJ"), color="sample_new", palette=pal_smp, show=False, legend_fontsize=5, s=size_dot, ax=axs[1,3], title="D08 (T cells)", legend_loc=None)
                umap1_4=sc.pl.umap(adata[adata.obs.cell_VDJ=="VDJ"], mask_obs=(adata[adata.obs.cell_VDJ=="VDJ"].obs["sample_new_VDJ"]=="D15_VDJ"), color="sample_new", palette=pal_smp, show=False, legend_fontsize=5, s=size_dot, ax=axs[1,4], title="D15 (T cells)", legend_loc=None)
                umap1_5=sc.pl.umap(adata[adata.obs.cell_VDJ=="VDJ"], mask_obs=(adata[adata.obs.cell_VDJ=="VDJ"].obs["sample_new_VDJ"]=="D22_VDJ"), color="sample_new", palette=pal_smp, show=False, legend_fontsize=5, s=size_dot, ax=axs[1,5], title="D22 (T cells)", legend_loc=None)
                umap1_6=sc.pl.umap(adata[adata.obs.cell_VDJ=="VDJ"], mask_obs=(adata[adata.obs.cell_VDJ=="VDJ"].obs["sample_new_VDJ"]=="D28_VDJ"), color="sample_new", palette=pal_smp, show=False, legend_fontsize=5, s=size_dot, ax=axs[1,6], title="D28 (T cells)", legend_loc=None)
                umap1_7=sc.pl.umap(adata[adata.obs.cell_VDJ=="VDJ"], mask_obs=(adata[adata.obs.cell_VDJ=="VDJ"].obs["sample_new_VDJ"]=="D42_VDJ"), color="sample_new", palette=pal_smp, show=False, legend_fontsize=5, s=size_dot, ax=axs[1,7], title="D42 (T cells)")
                
                umap2=sc.pl.umap(adata[adata.obs.cell_VDJ=="VDJ"], color="expr_GCAR", groups=["pos"], show=False, palette=pal_expr_cat_pos, legend_fontsize=5, s=size_dot, ax=axs[2,0], title="T cells (GCAR-pos)", legend_loc=None)
                umap2_1=sc.pl.umap(adata[adata.obs["sample_new_VDJ"]=="Enriched_VDJ"], color="expr_GCAR", groups=["pos"], show=False, palette=pal_expr_cat_pos, legend_fontsize=5, s=size_dot, ax=axs[2,1], title="T cells (GCAR-pos) - Enriched", legend_loc=None)
                umap2_2=sc.pl.umap(adata[adata.obs["sample_new_VDJ"]=="Harvest_VDJ"], color="expr_GCAR", groups=["pos"], show=False, palette=pal_expr_cat_pos, legend_fontsize=5, s=size_dot, ax=axs[2,2], title="T cells (GCAR-pos) - Harvest", legend_loc=None)
                umap2_3=sc.pl.umap(adata[adata.obs["sample_new_VDJ"]=="D08_VDJ"], color="expr_GCAR", groups=["pos"], show=False, palette=pal_expr_cat_pos, legend_fontsize=5, s=size_dot, ax=axs[2,3], title="T cells (GCAR-pos) - D08", legend_loc=None)
                umap2_4=sc.pl.umap(adata[adata.obs["sample_new_VDJ"]=="D15_VDJ"], color="expr_GCAR", groups=["pos"], show=False, palette=pal_expr_cat_pos, legend_fontsize=5, s=size_dot, ax=axs[2,4], title="T cells (GCAR-pos) - D15", legend_loc=None)
                umap2_5=sc.pl.umap(adata[adata.obs["sample_new_VDJ"]=="D22_VDJ"], color="expr_GCAR", groups=["pos"], show=False, palette=pal_expr_cat_pos, legend_fontsize=5, s=size_dot, ax=axs[2,5], title="T cells (GCAR-pos) - D22", legend_loc=None)
                umap2_6=sc.pl.umap(adata[adata.obs["sample_new_VDJ"]=="D28_VDJ"], color="expr_GCAR", groups=["pos"], show=False, palette=pal_expr_cat_pos, legend_fontsize=5, s=size_dot, ax=axs[2,6], title="T cells (GCAR-pos) - D28", legend_loc=None)
                umap2_7=sc.pl.umap(adata[adata.obs["sample_new_VDJ"]=="D42_VDJ"], color="expr_GCAR", groups=["pos"], show=False, palette=pal_expr_cat_pos, legend_fontsize=5, s=size_dot, ax=axs[2,7], title="T cells (GCAR-pos) - D42")
                                
                umap3=sc.pl.umap(adata[adata.obs.cell_VDJ=="VDJ"], color="expr_GCAR", groups=["neg"], show=False, palette=pal_expr_cat_neg, legend_fontsize=5, s=size_dot, ax=axs[3,0], title="T cells (GCAR-neg)", legend_loc=None)
                umap3_1=sc.pl.umap(adata[adata.obs["sample_new_VDJ"]=="Enriched_VDJ"], color="expr_GCAR", groups=["neg"], show=False, palette=pal_expr_cat_neg, legend_fontsize=5, s=size_dot, ax=axs[3,1], title="T cells (GCAR-neg) - Enriched", legend_loc=None)
                umap3_2=sc.pl.umap(adata[adata.obs["sample_new_VDJ"]=="Harvest_VDJ"], color="expr_GCAR", groups=["neg"], show=False, palette=pal_expr_cat_neg, legend_fontsize=5, s=size_dot, ax=axs[3,2], title="T cells (GCAR-neg) - Harvest", legend_loc=None)
                umap3_3=sc.pl.umap(adata[adata.obs["sample_new_VDJ"]=="D08_VDJ"], color="expr_GCAR", groups=["neg"], show=False, palette=pal_expr_cat_neg, legend_fontsize=5, s=size_dot, ax=axs[3,3], title="T cells (GCAR-neg) - D08", legend_loc=None)
                umap3_4=sc.pl.umap(adata[adata.obs["sample_new_VDJ"]=="D15_VDJ"], color="expr_GCAR", groups=["neg"], show=False, palette=pal_expr_cat_neg, legend_fontsize=5, s=size_dot, ax=axs[3,4], title="T cells (GCAR-neg) - D15", legend_loc=None)
                umap3_5=sc.pl.umap(adata[adata.obs["sample_new_VDJ"]=="D22_VDJ"], color="expr_GCAR", groups=["neg"], show=False, palette=pal_expr_cat_neg, legend_fontsize=5, s=size_dot, ax=axs[3,5], title="T cells (GCAR-neg) - D22", legend_loc=None)
                umap3_6=sc.pl.umap(adata[adata.obs["sample_new_VDJ"]=="D28_VDJ"], color="expr_GCAR", groups=["neg"], show=False, palette=pal_expr_cat_neg, legend_fontsize=5, s=size_dot, ax=axs[3,6], title="T cells (GCAR-neg) - D28", legend_loc=None)
                umap3_7=sc.pl.umap(adata[adata.obs["sample_new_VDJ"]=="D42_VDJ"], color="expr_GCAR", groups=["neg"], show=False, palette=pal_expr_cat_neg, legend_fontsize=5, s=size_dot, ax=axs[3,7], title="T cells (GCAR-neg) - D42")

                ## for [by Sample] plots ##
                legend_texts = umap0_7.get_legend().get_texts()
                for legend_text in legend_texts:
                    if legend_text.get_text() == "NA":
                        legend_text.set_text("other samples")

                legend_texts = umap1_7.get_legend().get_texts()
                for legend_text in legend_texts:
                    if legend_text.get_text() == "NA":
                        legend_text.set_text("non-T cells")
                
                ## for [GCAR-pos] plots ##
                ## change "NA" values in the UMAP plot to a specific value
                legend_texts = umap2_7.get_legend().get_texts()
                for legend_text in legend_texts:
                    if legend_text.get_text() == "NA":
                        legend_text.set_text("neg")
                
                ## for [GCAR-neg] plots ##
                ## change "NA" values in the UMAP plot to a specific value
                legend_texts = umap3_7.get_legend().get_texts()
                for legend_text in legend_texts:
                    if legend_text.get_text() == "NA":
                        legend_text.set_text("pos")                        
                
                plt.savefig(FN_PLOT_UMAP, orientation='landscape', bbox_inches='tight', pad_inches=0.1)
                if (os.path.exists(FN_PLOT_UMAP)):
                    print("[DONE] save as PDF: "+FN_PLOT_UMAP)
    # ---
# ===