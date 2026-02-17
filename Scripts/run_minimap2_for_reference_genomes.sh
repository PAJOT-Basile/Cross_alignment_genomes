#! /bin/bash

#SBATCH --mem=250G
#SBATCH --partition=fast
#SBATCH --job-name=Minimap

# Import variables from command line
ref_genome_1="${1}"
ref_genome_2="${2}"
out_name="${3}"


# Run minimap2
/shared/software/miniconda/envs/minimap2-2.28/bin/minimap2 "${ref_genome_1}" "${ref_genome_2}" > "${out_name}"
echo "Done!"