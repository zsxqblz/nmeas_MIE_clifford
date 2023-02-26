include("dependencies.jl")
include("sim.jl")
include("exp.jl")

const n_meas_start = parse(Int64,ARGS[1])
const n_meas_end = parse(Int64,ARGS[2])
const n_meas_step = parse(Int64,ARGS[3])
const n_Asites = 10
const n_Bsites = 200
const n_Csites = 10
const nsim = 100
const depth_start = parse(Int64,ARGS[4])
const depth_end = parse(Int64,ARGS[5])
const depth_step = parse(Int64,ARGS[6])
const file_name = ARGS[7]

const n_meas_l = floor.(Int,collect(range(n_meas_start,stop=n_meas_end,step=n_meas_step)))
const n_meas_length = length(n_meas_l)
const depth_l = floor.(Int,collect(range(depth_start,stop=depth_end,step=depth_step)))
const depth_length = length(depth_l)
cmi_ave,cmi_std = scanNmeasDepthHPBW(n_Asites,n_Bsites,n_Csites,n_meas_start,n_meas_end,n_meas_step,depth_start,depth_end,depth_step,nsim,true)
save2DData(n_meas_l,depth_l,cmi_ave,cmi_std,file_name)
