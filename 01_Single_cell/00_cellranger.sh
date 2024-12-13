#!/bin/bash
###################################
#### Cell Ranger (scRNA/scTCR) ####
#### project: [GCAR1]          ####
###################################
# version 0.5.5 (2024-12-13)
# by hyojinsong

# ---
## running switch ##
MULTI=1
# ---

## variable ##
SP=$1
PRJ=$2
SMP=$3

## project ##
PRJID="GCAR1_PT01"

## library info ##
# [1000263] "Chromium Next GEM Single Cell 5' Kit v2, 16 rxns"
# [1000286] "Chromium Next GEM Chip K Single Cell Kit, 48 rxns"
# [1000252] "Chromium Single Cell Human TCR Amplification Kit, 16 rxns"
# [Single or Dual Indexing] = "" #note: need to check the detail
# [ssDNA or dsDNA] = "" #note: need to check the detail

## sequencer info ##
# [Sequencing Service (instrument and run)] = "NovaSeq 200 cycle S4 v1.5" -- #example from the earlier code
# [Sequencing Service (instrument and run)] = "" -- #note: need to check the detail

## custom reference required ##
## GCAR construct info:
## among the entire <*.dna> file:
## - from: 3428 (CD8a Leader)
## - to: 4930 (CD3z)
## sequence length: 1,503 bp

#### species info ####
SP="hs38_GCAR"
PRJSET="GCAR1_CLIC-YYC"


#### directory structure ####
## default directories - need to modify these directory variables regarding to your running environment
HOMEDIR="/Users/your_username"
DATADIR="$HOMEDIR/data"
REFDIR="$HOMEDIR/reference"
RESDIR="$HOMEDIR/result"
TOOLDIR="$HOMEDIR/tool/cellranger-7.0.1"

## set the working directory
DIR="$DATADIR/$PRJSET"
DIR_REF="$REFDIR/cellranger"
DIR_GCAR="$RESDIR/GCAR1"
DIR_DF="$DATADIR/GCAR1_PT01"    

## set an input (i.e., FASTQ) directory
DIR_FQ_ROOT="${DATADIR}/${PRJID}/GEX_VDJ"

## set an output directory
if [ "$MULTI" -eq 1 ]; then
    OUTDIR_RUN="${RESDIR}/cellranger/${PRJID}/runs_${SP}_multi"
fi

## make the output directory
if [ ! -d $OUTDIR_RUN ]; then
    mkdir -p $OUTDIR_RUN
    echo "[DONE] create a directory (for output): ${OUTDIR_RUN}"
else
    echo "[SKIP] directory exists: ${OUTDIR_RUN}"
fi
# ---


## in/output ##
## set genome name ##
if [ "$SP" == "hs38_GCAR" ]; then
    GNAME="customref-gex-GRCh38-2024-A_GCAR"
    DIR_REFDATA="${REFDIR}/cellranger/customref-gex-GRCh38-2024-A_GCAR"
    GNAME_VDJ="customref-vdj-GRCh38-2024-A_GCAR"
    DIR_REFDATA_VDJ="${REFDIR}/cellranger/customref-vdj-GRCh38-2024-A_GCAR"
fi


## CONFIG file ##
if [ "$MULTI" -eq 1 ]; then
    CSV_CONFIG="${DIR_DF}/config_multi_${PRJID}_${SMP}.csv"
fi
# ---


## STEP 0 | make a custom reference - using bash/command-line ####
# FA_GCAR_ORIG="$REFDIR/GCAR_construct_1503bp_from-CD8aLeader_to-CD3z.fa"
# FA_GCAR="$REFDIR/GCAR.fa"
# cat $FA_GCAR_ORIG | sed 's/G\-CLIC\ \CTA\ \(CD8a\ \Leader\ \-\ \CD3z\)\ \ \(1503\ \bp\)/GCAR' >GCAR_construct_1503bp.fa

# To find the number of bases in this sequence,
# we will use the grep -v "^>" command to search all lines that don't start with the > character, which removes line returns with tr -d "\n" so they aren't counted, and then counts the number of characters with the command wc -c.
# Each command is sent to the next step with the pipe | command.
# cat GCAR.fa | grep -v "^>" | tr -d "\n" | wc -c #checkpoint #1503

# echo -e 'GCAR\tunknown\texon\t1\t1503\t.\t+\t.\tgene_id "GCAR"; transcript_id "GCAR"; gene_name "GCAR"; gene_biotype "protein_coding";' >GCAR.gtf

## [FA] file - fasta ##
# cd $REFDIR/cellranger/customref-gex-GRCh38-2024-A_GCAR/fasta
# cp genome.fa genome_GCAR.fa
# cat $REFDIR/GCAR.fa >>genome_GCAR.fa
# grep ">" genome_GCAR.fa #checkpoint
# cat $REFDIR/GCAR.fa >>genome.fa #can overwrite to this file since it was already copied from the original ref file
# grep ">" genome.fa

## [GTF] file - genes ##
# cd $REFDIR/cellranger/customref-gex-GRCh38-2024-A_GCAR/genes
# cat $DIR_GCAR/ref/GCAR.gtf >>genes.gtf
# tail genes.gtf #checkpoint

## make an index file: <*.fa.fai> - after adding the GCAR construct sequence
# samtools faidx ref.fasta [region1 [...]]
# samtools --fai-idx FILE
# FA="$REFDIR/cellranger/customref-gex-GRCh38-2024-A_GCAR/fasta/genome.fa"
# samtools faidx $FA


###############
#### multi ####
###############

## STEP 1 | run <cellranger multi> ####
## - usage:
# mkdir /home/jdoe/runs
# cd /home/jdoe/runs
# cellranger multi --id=sample345 --csv=/home/jdoe/sample345.csv 

if [ "$MULTI" -eq 1 ]; then
    cd $OUTDIR_RUN
    
    cmd_multi="cellranger multi \
        --id=$SMP \
        --csv=$CSV_CONFIG
    "
    echo $cmd_multi
    eval $cmd_multi
    echo "--- [DONE] Step 1: cellranger multi"
fi
