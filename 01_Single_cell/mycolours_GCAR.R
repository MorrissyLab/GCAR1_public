##############################
#### my colours           ####
#### - project: [ GCAR1 ] ####
##############################
# version: 1.12 (2024-12-13)
# by hyojinsong
today <- Sys.Date()
# ---

#### package ####
mypackages <- c(
  "RColorBrewer", "pals" #for colour palette
  ,"stringr" #for [str_pad] function
)

## [RColorBrewer] ##
# install.packages("RColorBrewer")

## [pals] ##
# install.packages("pals")

lapply(mypackages, library, character.only = TRUE)
# ---


#### GCAR ####
## #1 | longitudinal sample ##
## - note: created based on [Dark2] (n=7)
col_gcar_e <- "#666666" #dark_grey
col_gcar_h <- "#bf5b17" #brown
col_gcar_d08 <- "#a6d854" #light_green
col_gcar_d15 <- "#e7298a" #pink
col_gcar_d22 <- "#e6ab02" #yellow
col_gcar_d28 <- "#1b9e77" #dark_green
col_gcar_d42 <- "#386cb0" #blue
smpColor <- setNames(c(col_gcar_p,col_gcar_t,col_gcar_d08,col_gcar_d15,col_gcar_d22,col_gcar_d28,col_gcar_d42), c("Enriched_Apheresis","GCAR1_Product","D08","D15","D22","D28","D42")) #sample_new
# ---


## #2 | expression ##
## - source: [https://xkcd.com/color/rgb/]
# pal_expr=sns.blend_palette(["lightgrey", sns.xkcd_rgb["scarlet"]], as_cmap=True)
col_pos <- "#be0119" #scarlet
col_neg <- "#d8dcd6" #lightgrey
exprColor <- setNames(c(col_pos,col_neg), c("pos","neg"))
# ---


## #3 | Seurat cell cluster ##
# nClst <- length(unique(sobj_flt_umap$seurat_clusters))
# if (nClst<=25) {
#     clstColor <- setNames(cols25(n = 25)[1:nClst], sort(unique(as.factor(sobj_flt_umap$seurat_clusters))))
# } else if (nClst>25) {
#     clstColor <- setNames(c(cols25(n = 25),kelly(n=22)[3:22])[1:nClst], sort(unique(as.factor(sobj_flt_umap$seurat_clusters))))
# } else if (nClst>45) {
#     clstColor <- setNames(c(cols25(n = 25),kelly(n=22))[1:nClst], sort(unique(as.factor(sobj_flt_umap$seurat_clusters))))
# }

# col_C_00   <- "#1F78C8"
# col_C_01   <- "#ff0000"
# col_C_02   <- "#33a02c"
# col_C_03   <- "#6A33C2"
# col_C_04   <- "#ff7f00"
# col_C_05   <- "#565656"
# col_C_06   <- "#FFD700"
# col_C_07   <- "#a6cee3"
# col_C_08   <- "#FB6496"
# col_C_09   <- "#b2df8a"
# col_C_10   <- "#CAB2D6"
# col_C_11   <- "#FDBF6F"
# col_C_12   <- "#999999"
# col_C_13   <- "#EEE685"
# col_C_14   <- "#C8308C"
# col_C_15   <- "#FF83FA"
# col_C_16   <- "#C814FA"
# col_C_17   <- "#0000FF"
# col_C_18   <- "#36648B"
# col_C_19   <- "#00E2E5"
# col_C_20   <- "#00FF00"
# col_C_21   <- "#778B00"
# col_C_22   <- "#BEBE00"
# col_C_23   <- "#8B3B00"
# col_C_24   <- "#A52A3C"
# col_C_25   <- "#F2F3F4"
# col_C_26   <- "#222222"
# col_C_27   <- "#F3C300"
# col_C_28   <- "#875692"
# col_C_29   <- "#F38400" #excluded from the analyses
# col_C_30   <- "#A1CAF1"
# col_C_31   <- "#BE0032"
# col_C_32   <- "#C2B280"
# col_C_33   <- "#848482"
# col_C_34   <- "#008856"
# col_C_35   <- "#E68FAC"
# col_C_36   <- "#0067A5"
# col_C_37   <- "#F99379"
# col_C_38   <- "#604E97"
# col_C_39   <- "#F6A600"
# col_C_40   <- "#B3446C"
# col_C_41   <- "#DCD300"
# col_C_42   <- "#882D17"
# col_C_43   <- "#8DB600" #excluded from the analyses
# col_C_44   <- "#654522"
# col_C_45   <- "#E25822"
# col_C_46   <- "#2B3D26"

nClst <- 47 #res=2.4
list_clst <- paste0("C_",str_pad(seq(0,(nClst-1),by=1),2,pad="0"))
clstColor <- setNames(c(cols25(n = 25),kelly(n=22))[1:nClst], list_clst) #up to 47 #from [C_00] to [C_46]
# ---


## #4 | cell class ##
list_CellClass_specific <- sort(unique(obs_c__exCC43_CC29$Specific_cell_class))
nCellClass_specific <- length(list_CellClass_specific) #23
ccColor21 <- setNames( cols25(n = 25)[1:nCellClass_specific], list_CellClass_specific ) #21
#                 ccColor21
# B-cells           #1F78C8
# CD4               #ff0000
# CD4_Treg          #33a02c
# CD4em             #6A33C2
# CD4naive          #ff7f00
# CD4naive_helper   #565656
# CD8e              #FFD700
# CD8e_td           #a6cee3
# CD8em             #FB6496
# CD8em_MAIT        #b2df8a
# CD8em_NK          #CAB2D6
# CD8em_td          #FDBF6F
# CD8em_td_MAIT     #999999
# DCs               #EEE685
# HSC               #C8308C
# Lymphocyte        #FF83FA
# MAIT              #C814FA
# Monocyte          #0000FF
# ncMonocyte        #36648B
# NK                #00E2E5
# NKT               #00FF00
# ---