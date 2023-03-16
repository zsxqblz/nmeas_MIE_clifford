#!/bin/bash
#SBATCH --job-name=nd4
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G         # memory per cpu-core (4G is default)
#SBATCH --time=23:59:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --array=0-499
#SBATCH --mail-user=yz4281@princeton.edu

julia scan_nmeas_depth.jl 470 500 1 0 250 5 10 500 10 20 data/230313/230313_nd4_$SLURM_ARRAY_TASK_ID