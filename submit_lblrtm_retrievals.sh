#!/bin/bash
#SBATCH --account=rrg-yihuang-ad
#SBATCH --job-name=submit_lblrtm_jobs
#SBATCH --output=logs/launcher_%j.out
#SBATCH --error=logs/launcher_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00

# Load modules
module load matlab

# Create log directory
mkdir -p logs

# Base directory where profiles live
PROFILE_DIR="/home/binmenja/direct/field_campaigns/whafffers"

# Get list of profile files
FILES=(${PROFILE_DIR}/profile_2025*.mat)

echo "Submitting ${#FILES[@]} LBLRTM retrieval jobs..."

for FILE in "${FILES[@]}"; do
    BASENAME=$(basename "$FILE" .mat)
    DATESTRING=$(echo "$BASENAME" | sed 's/profile_//' | sed 's/UTC//')
    JOBNAME="retrv_${DATESTRING}"

    sbatch <<EOF
#!/bin/bash
#SBATCH --account=rrg-yihuang-ad
#SBATCH --job-name=$JOBNAME
#SBATCH --output=logs/${JOBNAME}_%j.out
#SBATCH --error=logs/${JOBNAME}_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=12:00:00

module load matlab
matlab -nodesktop -nosplash -r "try; run_multiple_lblrtm_retrievals('$DATESTRING'); catch ME; disp(getReport(ME)); end; exit"
EOF

done

echo "All LBLRTM retrieval jobs submitted."
