#!/bin/bash
#SBATCH --job-name=scan_dy
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G         # memory per cpu-core (4G is default)
#SBATCH --time=23:59:59          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --array=0-99
#SBATCH --mail-user=yz4281@princeton.edu

let i=4

let id=400+i
let date=250322

julia scan_dy.jl $i data/${date}/${date}_n${id}_$SLURM_ARRAY_TASK_ID