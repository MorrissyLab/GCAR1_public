#!/usr/bin/env bash

#SBATCH --job-name=SpaceRanger
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=245G
#SBATCH --signal=2
#SBATCH --no-requeue
#SBATCH --partition=cpu2023,cpu2022-bf24,cpu2022,cpu2021-bf24,cpu2021,cpu2019-bf05,cpu2019,cpu2017-bf05,parallel,single,sherlock,bigmem
#SBATCH --time=05:00:00
#SBATCH -o %j.stdout.txt
#SBATCH -e %j.stderr.txt
#SBATCH --export=ALL

####### Set environment variables ###############
BASEDIR="../datasets/VisiumCytAssist"
TRANSCRIPTOME="../datasets/GCAR_reference/GCAR1_reference_V2"
PROBE_SET="../datasets/GCAR_reference/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A_V2.csv"

ID=(SPSQ_biopsy_S2_DS_RA_M)
SAMPLE=(KHSA2_M)
FASTQ_PATH=(${BASEDIR}/"merged_biopsy/")
SLIDE=(H1-Q9DYC4M)
AREA=(D1)

####### Run your script #########################
spaceranger count --id=${ID} \
                  --sample=${SAMPLE} \
                  --description="${ID} sample." \
                  --transcriptome=${TRANSCRIPTOME} \
                  --probe-set=${PROBE_SET} \
                  --fastqs=${FASTQ_PATH} \
                  --cytaimage=${BASEDIR}/"assay_CAVG10599_2024-08-14_12-30-05_H1-Q9DYC4M_1723662702_CytAssist/CAVG10599_2024-08-14_13-11-42_2024-08-14_12-30-05_H1-Q9DYC4M_D1_FHS23-020962-A2.tif" \
                  --image=${BASEDIR}/"scanned images/FHR_A2_40X.tif" \
                  --slide=${SLIDE} \
                  --area=${AREA} \
                  --localcores=${SLURM_CPUS_PER_TASK} \
                  --localmem=245 \
                  --create-bam=true \
                  --custom-bin-size=24 \
                  --loupe-alignment=${BASEDIR}/"H1-Q9DYC4M-D1-fiducials-image-registration_V2.json"
