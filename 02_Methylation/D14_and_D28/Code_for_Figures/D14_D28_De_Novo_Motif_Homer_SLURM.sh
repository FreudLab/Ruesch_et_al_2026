#!/bin/bash
#SBATCH --account=PAS2527
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --partition=nextgen
#SBATCH --mail-type=ALL
#SBATCH --output=HOMER_DeNovoPeaks.slurm-%j.out

# Load conda & activate homer env
module load miniconda3/24.1.2-py310
source activate homer

# Add HOMER to PATH
export PATH=$PATH:/fs/ess/PAS2527/HOMER/bin/

# Define working directory
WORKDIR=$HOME/DNAmethylation/NKDI_NK_IV_Analysis/De_Novo_Motifs_Homer/NK_15Cv15D
TMP_WORKDIR=$TMPDIR/homer_loop
mkdir -p $TMP_WORKDIR
cd $TMP_WORKDIR

# List of foreground/background pairs
pairs=(
  "CD16Neg_hyper"
  "CD16Neg_hypo"
  "CD16Pos_hyper"
  "CD16Pos_hypo"
  "D14_hyper"
  "D14_hypo"
)

# Copy inputs to TMPDIR
for p in "${pairs[@]}"; do
  cp $WORKDIR/fg_${p}.txt $TMP_WORKDIR/
  cp $WORKDIR/bg_${p}.txt $TMP_WORKDIR/
done

# Loop through pairs
for p in "${pairs[@]}"; do
  echo ">>> Running HOMER for ${p}"
  
  /usr/bin/time -v findMotifsGenome.pl fg_${p}.txt hg38 out_${p} \
    -bg bg_${p}.txt \
    -size given \
    -cpg \
    -len 8,10,12 \
    -p 40 \
    -mask \
    2>&1 | tee homer_${p}.log
done

# Copy results back to WORKDIR
cp -r out_* $WORKDIR/
cp homer_*.log $WORKDIR/

# Clean up
rm -rf $TMP_WORKDIR

echo "All HOMER runs completed."

