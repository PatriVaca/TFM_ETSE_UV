#!/bin/bash

#SBATCH --job-name=synhc_MAST

#SBATCH --partition=bigmem

#SBATCH --cpus-per-task 10

#SBATCH --mem 230G

#SBATCH --output=../processed_data/synhc_MAST_FM.txt

#SBATCH --error=../processed_data/synhc_MAST_error_FM.txt

start=$(date +%s.%N)

# --------------------------------

module load R/4.2.2-foss-2021b

echo "Vamos a lanzar MAST"
# Rscript MAST_1.R
# Rscript MAST_2_F.R
# Rscript MAST_3_M.R 
Rscript MAST_4_FM.R

# --------------------------------

end=$(date +%s.%N)
runtime=$(echo "$end - $start" | bc)

hours=$(echo "scale=0; $runtime/3600" | bc)
remaining=$(echo "$runtime%3600" | bc)
minutes=$(echo "scale=0; $remaining/60" | bc)
seconds=$(echo "$remaining%60" | bc)

printf "Tiempo de ejecución: %02d:%02d:%09.6f\n" $hours $minutes $seconds
