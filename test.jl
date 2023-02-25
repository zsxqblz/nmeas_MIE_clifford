include("dependencies.jl")
include("sim.jl")
include("exp.jl")

for i = 11:20
    n_Asites = i
    n_Bsites = 40
    n_Csites = i
    nsim = 1000
    n_meas_start = 2
    n_meas_end = 40
    n_meas_step = 2

    save_idx_start = 171
    n_meas_l = floor.(Int,collect(range(n_meas_start,stop=n_meas_end,step=n_meas_step)))
    cmi_ave,cmi_std = scanNmeasHP(n_Asites,n_Bsites,n_Csites,n_meas_start,n_meas_end,n_meas_step,nsim,true)
    save1DData(n_meas_l,cmi_ave,cmi_std,string("data/230213_",save_idx_start+i))
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