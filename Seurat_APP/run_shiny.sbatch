#!/bin/bash
#SBATCH --job-name=scshiny
#SBATCH --output=scshiny_%j.out
#SBATCH --error=scshiny_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=interactive  # Change to your appropriate partition

# Load the singularity/apptainer module if needed (comment out if not used on your HPC)
module load singularity

# The port that Shiny server will use INSIDE the container
PORT=3838

# Print info for user
echo "Launching Shiny app at port $PORT in singularity container"
echo "Bind-mounting current directory as /data in the container."

# Run the container with the app, binding the present working directory for file uploads/downloads
singularity exec --bind $(pwd):/data \
    --env SHINY_PORT=$PORT \
    sc_shiny.sif shiny-server

echo "If running on a login/interactive node, or after job starts, you can forward the port like:"
echo "ssh -N -L 3838:localhost:3838 your_hpc_login"

echo "Then, on your laptop, open http://localhost:3838 in your browser."

