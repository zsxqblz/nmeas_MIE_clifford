include("dependencies.jl")
include("sim.jl")
include("exp.jl")
include("exp_2D_2norm.jl")

for i in [4]
    dx = 10
    nsim = 10000
    dy_start = 2
    dy_end = 20
    dy_step = 2
    noise = 0.02+0.005*(i-1)
    depth = 10

    save_idx_start = 2600
    dy_l = floor.(Int,collect(range(dy_start,stop=dy_end,step=dy_step)))
    dist_ave_arr, dist_std_arr = scandY2D(sim2DNoisyBW,dx,dy_start,dy_end,dy_step,noise,depth,nsim,true)
    save1DData(dy_l,dist_ave_arr, dist_std_arr,string("data/250319/250319_",save_idx_start+i))
end