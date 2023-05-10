#!/bin/bash
#SBATCH --job-name=2D
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G         # memory per cpu-core (4G is default)
#SBATCH --time=21:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --array=0-99
#SBATCH --mail-user=yz4281@princeton.edu

let i=4
let dx=10*i
let dy=10*i
let n_meas_start=0
let n_meas_end=dx*dy
let n_meas_step=i*i
let depth_start=0
let depth_end=2*dx
let depth_step=2*i
let nsim=100

let id=i
let date=230507

julia scan_nmeas_depth_2D.jl $n_meas_start $n_meas_end $n_meas_step $depth_start $depth_end $depth_step $dx $dy $nsim data/${date}/${date}_n${id}_$SLURM_ARRAY_TASK_ID