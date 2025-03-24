include("dependencies.jl")
include("sim.jl")
include("exp.jl")
include("exp_2D_2norm.jl")

const i = parse(Int64,ARGS[1])
const file_name = ARGS[2]

dx = 10
nsim = 10000
dy_start = 2
dy_end = 20
dy_step = 2
noise = 0.01+0.05*(i)
depth = 10

dy_l = floor.(Int,collect(range(dy_start,stop=dy_end,step=dy_step)))
dist_ave_arr, dist_std_arr = scandY2D(sim2DNoisyBW,dx,dy_start,dy_end,dy_step,noise,depth,nsim,true)
save1DData(dy_l,dist_ave_arr, dist_std_arr,file_name)