include("dependencies.jl")
include("sim.jl")
include("exp.jl")
include("exp_1D.jl")

for i = 1:5
    n_Asites = 10
    n_Bsites = 20*i
    n_Csites = 10
    nsim = 1000
    n_meas_start = 0
    n_meas_end = n_Bsites
    n_meas_step = 2*i

    save_idx_start = 10
    n_meas_l = floor.(Int,collect(range(n_meas_start,stop=n_meas_end,step=n_meas_step)))
    cmi_ave,cmi_std = scanNmeasHPP(n_Asites,n_Bsites,n_Csites,n_meas_start,n_meas_end,n_meas_step,nsim,true)
    save1DData(n_meas_l,cmi_ave,cmi_std,string("data/230317/230317_",save_idx_start+i))
end

let
    n_Asites = 10
    n_Bsites = 100
    n_Csites = 10
    nsim = 100
    n_meas_start = 0
    n_meas_end = 100
    n_meas_step = 2
    depth_start = 0
    depth_end = 5
    depth_step = 1

    save_idx = 3
    n_meas_l = floor.(Int,collect(range(n_meas_start,stop=n_meas_end,step=n_meas_step)))
    n_meas_length = length(n_meas_l)
    depth_l = floor.(Int,collect(range(depth_start,stop=depth_end,step=depth_step)))
    depth_length = length(n_meas_l)
    cmi_ave,cmi_std = scanNmeasDepthHPBDiff(n_Asites,n_Bsites,n_Csites,n_meas_start,n_meas_end,n_meas_step,depth_start,depth_end,depth_step,nsim,true)
    save2DData(n_meas_l,depth_l,cmi_ave,cmi_std,string("data/230322/230322_",save_idx))
end

let
    run(`clear`)
    n_Asites = 1
    n_Bsites = 10
    n_Csites = 1
    nsim = 1000
    n_meas_start = 0
    n_meas_end = 10
    n_meas_step = 1
    depth_start = 0
    depth_end = 80
    depth_step = 4

    save_idx = 7
    n_meas_l = floor.(Int,collect(range(n_meas_start,stop=n_meas_end,step=n_meas_step)))
    n_meas_length = length(n_meas_l)
    depth_l = floor.(Int,collect(range(depth_start,stop=depth_end,step=depth_step)))
    depth_length = length(n_meas_l)
    cmi_ave,cmi_std = scanNmeasDepth(simProdBWMeas,n_Asites,n_Bsites,n_Csites,n_meas_start,n_meas_end,n_meas_step,depth_start,depth_end,depth_step,nsim,true)
    # @time scanNmeasDepthHPBDiff(n_Asites,n_Bsites,n_Csites,n_meas_start,n_meas_end,n_meas_step,depth_start,depth_end,depth_step,nsim,true)
    save2DData(n_meas_l,depth_l,cmi_ave,cmi_std,string("data/230414/230414_",save_idx))
end

# exp_2D
let
    run(`clear`)
    dx = 10
    dy = 10
    nsim = 10
    n_meas_start = 0
    n_meas_end = dx*dy
    n_meas_step = 1
    depth_start = 0
    depth_end = 2*dx
    depth_step = 2

    save_idx = 1
    n_meas_l = floor.(Int,collect(range(n_meas_start,stop=n_meas_end,step=n_meas_step)))
    n_meas_length = length(n_meas_l)
    depth_l = floor.(Int,collect(range(depth_start,stop=depth_end,step=depth_step)))
    depth_length = length(n_meas_l)
    cmi_ave,cmi_std = scanNmeasDepth2D(sim2DBW,dx,dy,n_meas_start,n_meas_end,n_meas_step,depth_start,depth_end,depth_step,nsim,true)
    # @time scanNmeasDepthHPBDiff(n_Asites,n_Bsites,n_Csites,n_meas_start,n_meas_end,n_meas_step,depth_start,depth_end,depth_step,nsim,true)
    save2DData(n_meas_l,depth_l,cmi_ave,cmi_std,string("data/230507/230507_",save_idx))
end

# exp_1D
let
    run(`clear`)
    n_Asites = 10
    n_Bsites = 100
    n_Csites = 10
    nsim = 100
    n_meas_start = 70
    n_meas_end = 100
    n_meas_step = 3
    depth_start = 0
    depth_end = 60
    depth_step = 4

    save_idx = 2
    n_meas_l = floor.(Int,collect(range(n_meas_start,stop=n_meas_end,step=n_meas_step)))
    n_meas_length = length(n_meas_l)
    depth_l = floor.(Int,collect(range(depth_start,stop=depth_end,step=depth_step)))
    depth_length = length(n_meas_l)
    SA_ave_arr,SA_std_arr,SC_ave_arr,SC_std_arr,SAC_ave_arr,SAC_std_arr = scanNmeasDepth1D(simBWMeas1D,n_Asites,n_Bsites,n_Csites,n_meas_start,n_meas_end,n_meas_step,depth_start,depth_end,depth_step,nsim,true)
    save2DData(n_meas_l,depth_l,SA_ave_arr,SA_std_arr,string("data/230510/230510_",save_idx,"_SA"))
    save2DData(n_meas_l,depth_l,SC_ave_arr,SC_std_arr,string("data/230510/230510_",save_idx,"_SC"))
    save2DData(n_meas_l,depth_l,SAC_ave_arr,SAC_std_arr,string("data/230510/230510_",save_idx,"_SAC"))
end

let 
    reg = genInitABC(1,4,1)
    apply!(reg.stab,random_clifford(2),[2,3])
    reg.stab
end