#!/bin/bash
#SBATCH --job-name=1D
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G         # memory per cpu-core (4G is default)
#SBATCH --time=23:59:59          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --array=0-9
#SBATCH --mail-user=yz4281@princeton.edu

let i=15
let n_Asites=20+4*i
let n_Bsites=100+20*i
let n_Csites=20+4*i
let n_meas_start=45+10*i
let n_meas_end=55+10*i
let n_meas_step=1
let depth_start=20
let depth_end=120
# let depth_end=1
let depth_step=2
let nsim=10

let id=i
let date=240203

julia scan_neras_depth_1D.jl $n_meas_start $n_meas_end $n_meas_step $depth_start $depth_end $depth_step $n_Asites $n_Bsites $n_Csites $nsim data/${date}/${date}_n${id}_$SLURM_ARRAY_TASK_ID