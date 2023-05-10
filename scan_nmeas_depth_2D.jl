include("dependencies.jl")
include("sim.jl")
include("exp.jl")
include("exp_2D.jl")

const n_meas_start = parse(Int64,ARGS[1])
const n_meas_end = parse(Int64,ARGS[2])
const n_meas_step = parse(Int64,ARGS[3])
const depth_start = parse(Int64,ARGS[4])
const depth_end = parse(Int64,ARGS[5])
const depth_step = parse(Int64,ARGS[6])
const dx = parse(Int64,ARGS[7])
const dy = parse(Int64,ARGS[8])
const nsim = parse(Int64,ARGS[9])
const file_name = ARGS[10]

const n_meas_l = floor.(Int,collect(range(n_meas_start,stop=n_meas_end,step=n_meas_step)))
const n_meas_length = length(n_meas_l)
const depth_l = floor.(Int,collect(range(depth_start,stop=depth_end,step=depth_step)))
const depth_length = length(depth_l)
cmi_ave,cmi_std = scanNmeasDepth2D(sim2DBW,dx,dy,n_meas_start,n_meas_end,n_meas_step,depth_start,depth_end,depth_step,nsim,true)
save2DData(n_meas_l,depth_l,cmi_ave,cmi_std,file_name)
