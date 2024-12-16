#!/bin/bash

#SBATCH --job-name=cNMF_plots
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mincpus=1
#SBATCH --mem=36G
#SBATCH --signal=2
#SBATCH --no-requeue
#SBATCH --partition=cpu2023,cpu2022-bf24,cpu2022,cpu2021-bf24,cpu2021,cpu2019-bf05,cpu2019,cpu2017-bf05,parallel,single
#SBATCH --time=05:00:00
#SBATCH -o %A_%a.stdout.txt
#SBATCH -e %A_%a.stderr.txt
#SBATCH --export=ALL

####### Set environment variables ###############
source ~/software/init-conda
conda activate mosaic

RUN_NAME=$1
PLOT_DIRECTORY=$2
METADATA_FILE_CSV=$3
SAMPLES_FILE=$4
SR_OUTPUT=$5
RUN_RANGE=(${@:6})

RANK=${RUN_RANGE[${SLURM_ARRAY_TASK_ID}]}
USAGE_FILE=${RUN_NAME}"/"${RUN_NAME}".usages.k_"${RANK}".dt_2_0.consensus.txt"
SPECTRA=${RUN_NAME}"/"${RUN_NAME}".gene_spectra_score.k_"${RANK}".dt_2_0.txt"

ROWS=1
COLS=3
MAX_CUTOFF=0
SEP="_"

####### Run your script #########################
python plots.py ${RANK} ${USAGE_FILE} ${METADATA_FILE_CSV} ${SAMPLES_FILE} ${SR_OUTPUT} ${PLOT_DIRECTORY} ${ROWS} ${COLS} ${MAX_CUTOFF} ${SEP}
wait

# The combined plots for Visium HD data are multi-GB files, which can be difficult to open on any computer.
# python pdfcombine.py ${RANK} ${PLOT_DIRECTORY}
