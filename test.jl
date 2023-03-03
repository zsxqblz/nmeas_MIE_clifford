include("dependencies.jl")
include("sim.jl")
include("exp.jl")

for i = 2:6
    n_Asites = 10
    n_Bsites = 20*i
    n_Csites = 10
    nsim = 1000
    n_meas_start = n_Bsites-38
    n_meas_end = n_Bsites
    n_meas_step = 2*i

    save_idx_start = 10
    n_meas_l = floor.(Int,collect(range(n_meas_start,stop=n_meas_end,step=n_meas_step)))
    cmi_ave,cmi_std = scanNmeasTMB(n_Asites,n_Bsites,n_Csites,n_meas_start,n_meas_end,n_meas_step,nsim,true)
    save1DData(n_meas_l,cmi_ave,cmi_std,string("data/230226/230226_",save_idx_start+i))
end

let
    n_Asites = 10
    n_Bsites = 100
    n_Csites = 10
    nsim = 100
    n_meas_start = 80
    n_meas_end = 100
    n_meas_step = 2
    depth_start = 60
    depth_end = 100
    depth_step = 4

    save_idx = 2
    n_meas_l = floor.(Int,collect(range(n_meas_start,stop=n_meas_end,step=n_meas_step)))
    n_meas_length = length(n_meas_l)
    depth_l = floor.(Int,collect(range(n_meas_start,stop=n_meas_end,step=n_meas_step)))
    depth_length = length(n_meas_l)
    cmi_ave,cmi_std = scanNmeasDepthBW(n_Asites,n_Bsites,n_Csites,n_meas_start,n_meas_end,n_meas_step,depth_start,depth_end,depth_step,nsim,true)
    save2DData(n_meas_l,depth_l,cmi_ave,cmi_std,string("data/230225/230225_2d_",save_idx))
end

plot(n_meas_l,abs.(cmi_ave.+0.001),yaxis=:log)
# p = genIStr(2)
# p = S"X"
# p = Destabilizer(p)
# p = CliffordOperator(p)
# u = random_clifford(3)
# p ⊗ u
# reg = genInitABC(1,3,1)

# apply!(reg,u)
# p = zero(CliffordOperator,2)
# p[1,1] = (fa)

# p = genIOp(2)
# p ⊗ u


# reg = genInitHP(1,4,1)
# randBellMeasB(reg,1,4,1,2)
# applyRandClif(reg,1,3,1)
# reg
# floor.(Int,[1.1,1.2])
# cmi(reg,1,4,1)

string("data/230213_",1)

let 
    typeof(random_clifford(2))
end