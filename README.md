# GCAR1 (GPNMB-targeting CAR T therapy)


- Study Title: "Development and first-in-human CAR T therapy against the pathognomonic MiT-fusion driven protein GPNMB" (_manuscript submitted_)

  > _"We developed a GPNMB(glycoprotein non-metastatic B)-targeting CAR T therapy called GCAR1 that shows activity against patient-matched cells, organoids and xenograft models."_

---
## Data deposition
1. single-cell transcriptomic data (scRNA & TCR-seq)
   - `[GSE283385]`: Raw and processed single-cell 5' and V(D)J Enrichment data -- Raw and processed data include fastq files and CellRanger output, respectively, for 14 samples (7 scRNA-seq and 7 TCR-seq).

3. spatial transcriptomic data (Visium HD)
   - `[GSE282057]`: Raw and processed VisiumHD spatial gene expression data -- Raw and processed data include fastq files and SpaceRanger output respectively, for 4 samples.

5. whole genome & transcriptomic data (WGS & WTS)
   - `[EGAD00000000000]`: Whole genome and transcriptome sequencing data and deep sequencing panel data (accession number generation _in progress_)

- All other data are included in the main text, the Supplementary Information or are available from the authors, as are unique reagents used in this article. Source data (including the raw numbers for charts and graphs) are available in the Source Data file whenever possible. Source data are provided in this paper.


---
## Sub-directory / Code description
### 01_Single_cell/
   - `<00_cellranger.sh>`: code to generate _Cell Ranger_ output (from FASTQ files)
   - `01_GEX_scRNA/`: includes analysis & visualization codes ---
      - `<01_Seurat.R>`: scRNA-seq data analysis using _Seurat_'s built-in functions
      - `<02_SCimilarity.py>`: cell-level cell type annotation using _SCimilarity_'s built-in model & relevant visualization codes
      - `<02a_calc_minDistScore.R>`: function to calculate the "cluster label score" (See _Methods_ in the manuscript)
      - `<02b_sum_propCell_exprGCAR.R>`: visualization code for [Supplementary Figure #10]
   - `02_VDJ_scTCR/`: includes analysis & visualization codes ---
      - `<01_scRepertoire.R>`: scTCR-seq data analysis using _scRepertoire_
   - `<mycolours_GCAR.R>`: main colour palette used in the manuscript figures
   - `<ref/GCAR.fa>`: FASTA file used to generate our 'customized' reference for the _Cell Ranger_ run

### 02_Visium_HD/
   - `<02_analysis/01_cnmf_prep.R>`: Produce the count matrix, metadata table, and list of over-dispersed genes necessary for cNMF. 
   - `<02_analysis/02_cnmf_init.sh>`: Slurm script that, (1) Prepares the data for MosaicMPI factorization, (2) Submits another Slurm script to facilitate parallelized data factorization, (3) Factorization post-processing, and (4) Submits another Slurm script to create spatial plots.
   - `<02_analysis/03_cnmf_factorize.sh>`: Slurm script that parallelizes the factorization of the prepared data across N jobs, where N is 500 in the current script.
   - `<02_analysis/04_cnmf_usage_plots.sh>`: Slurm script that submits jobs to produce spatial plots (usage, gene expression, HE) for each rank ran.
   - `<02_analysis/05_plots.py>`: Python file that enables the creation of the spatial plots. Please note that the top of the file contains some arguments that can be enabled/disabled or modified to produce the desired plots.
   - `<02_analysis/06_pdfcombine.py>`: Python file that combines all spatial usage plots generated into one master PDF file. We recommend that reviewers look to the low-resolution versions of the spatial plots due to the file size of the PDF-based images, which may crash common PDF viewers due to the immense memory needed to open the file.
   - `</02_analysis/07_cnmf_hd_analysis.R>`: Includes code blocks to reproduce each figure supporting the spatial data analysis results. Please note, that an initial function has been added to facilitate the loading of the Seurat objects and cNMF usage results of the primary, biopsy cores, and xenograft samples for easy reproducibility. This data can be found in the GitHub as “.rds” files under the `<./02_Visium_HD/04_datasets>`. However, all the data and code necessary to reproduce all the figures from scratch have also been provided, with data uploaded to GEO and the code in this GitHub.
