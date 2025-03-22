include("dependencies.jl")
include("sim.jl")
include("exp.jl")
include("exp_2D_2norm.jl")


include("exp_1D.jl")
include("exp_2D.jl")

let
    run(`clear`)
    n_Asites = 10
    n_Bsites = 20
    n_Csites = 10
    nsites = n_Asites+n_Bsites+n_Csites
    s = genProdZState(nsites)
    boolStr = Vector{Bool}(undef,nsites)
    reg = Register(MixedDestabilizer(s), boolStr)
    U = random_clifford(nsites)
    apply!(reg,U)
    SAB = entanglement_entropy(reg.stab,collect(1:n_Asites+n_Bsites),Val(:rref))
    SBC = entanglement_entropy(reg.stab,collect(n_Asites+1:n_Asites+n_Bsites+n_Csites),Val(:rref))
    SB = entanglement_entropy(reg.stab,collect(n_Asites+1:n_Asites+n_Bsites),Val(:rref))
    SABC = entanglement_entropy(reg.stab,collect(1:n_Asites+n_Bsites+n_Csites),Val(:rref))
    @show SAB,SBC,SB,SABC
    randErasB(reg,n_Asites,n_Bsites,n_Csites,10)
    SAB = entanglement_entropy(reg.stab,collect(1:n_Asites+n_Bsites),Val(:rref))
    SBC = entanglement_entropy(reg.stab,collect(n_Asites+1:n_Asites+n_Bsites+n_Csites),Val(:rref))
    SB = entanglement_entropy(reg.stab,collect(n_Asites+1:n_Asites+n_Bsites),Val(:rref))
    SABC = entanglement_entropy(reg.stab,collect(1:n_Asites+n_Bsites+n_Csites),Val(:rref))
    @show SAB,SBC,SB,SABC
end

for i = 1:5
    n_Asites = 10
    n_Bsites = 20*i
    n_Csites = 10
    nsim = 100
    n_meas_start = 0
    n_meas_end = n_Bsites
    n_meas_step = 2

    save_idx_start = 0
    n_meas_l = floor.(Int,collect(range(n_meas_start,stop=n_meas_end,step=n_meas_step)))
    # cmi_ave,cmi_std = scanNmeasHPP(n_Asites,n_Bsites,n_Csites,n_meas_start,n_meas_end,n_meas_step,nsim,true)
    cmi_ave,cmi_std = scanNmeas(simProdHPClifMeas,n_Asites,n_Bsites,n_Csites,n_meas_start,n_meas_end,n_meas_step,nsim,true)
    save1DData(n_meas_l,cmi_ave,cmi_std,string("data/240206/240206_",save_idx_start+i))
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

# exp_2D entropy
let
    run(`clear`)
    dx = 5
    dy = 15
    dA = 5
    dC = 5
    dB = dy-dA-dC
    nsim = 100
    n_meas_start = 0
    n_meas_end = dx*dB
    n_meas_step = 3
    depth_start = 0
    depth_end = 2*dx
    depth_step = 2

    save_idx = 4
    n_meas_l = floor.(Int,collect(range(n_meas_start,stop=n_meas_end,step=n_meas_step)))
    n_meas_length = length(n_meas_l)
    depth_l = floor.(Int,collect(range(depth_start,stop=depth_end,step=depth_step)))
    depth_length = length(n_meas_l)
    SA_ave_arr,SA_std_arr,SC_ave_arr,SC_std_arr,SAC_ave_arr,SAC_std_arr = scanNmeasDepth2D(sim2DBW,dx,dy,dA,dC,n_meas_start,n_meas_end,n_meas_step,depth_start,depth_end,depth_step,nsim,true)
    save2DData(n_meas_l,depth_l,SA_ave_arr,SA_std_arr,string("data/230513/230513_",save_idx,"_SA"))
    save2DData(n_meas_l,depth_l,SC_ave_arr,SC_std_arr,string("data/230513/230513_",save_idx,"_SC"))
    save2DData(n_meas_l,depth_l,SAC_ave_arr,SAC_std_arr,string("data/230513/230513_",save_idx,"_SAC"))
end

# temp
let
    dx = 20
    dy = 20
    dA = 4
    dC = 4
    dB = dy-dA-dC
    nsim = 100
    n_meas_start = 0
    n_meas_end = dx*dB
    n_meas_step = 16
    depth_start = 0
    depth_end = 2*dx
    depth_step = 2

    save_idx = 4
    n_meas_l = floor.(Int,collect(range(n_meas_start,stop=n_meas_end,step=n_meas_step)))
    n_meas_length = length(n_meas_l)
    depth_l = floor.(Int,collect(range(depth_start,stop=depth_end,step=depth_step)))
    depth_length = length(n_meas_l)
    SA_ave_arr,SA_std_arr,SC_ave_arr,SC_std_arr,SAC_ave_arr,SAC_std_arr = scanNmeasDepth2D(sim2DBW,dx,dy,dA,dC,n_meas_start,n_meas_end,n_meas_step,depth_start,depth_end,depth_step,nsim,true)
    save2DData(n_meas_l,depth_l,SA_ave_arr,SA_std_arr,string("data/230513/230513_",save_idx,"_SA"))
    save2DData(n_meas_l,depth_l,SC_ave_arr,SC_std_arr,string("data/230513/230513_",save_idx,"_SC"))
    save2DData(n_meas_l,depth_l,SAC_ave_arr,SAC_std_arr,string("data/230513/230513_",save_idx,"_SAC"))
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

# exp_1D erasure
let
    i = 1
    run(`clear`)
    n_Asites = 20
    n_Bsites = 100 + 10*i
    n_Csites = 20
    nsim = 100
    n_meas_start = 30
    n_meas_end = 70
    n_meas_step = 5
    depth_start = 0
    depth_end = 60
    depth_step = 6

    save_idx = 3
    n_meas_l = floor.(Int,collect(range(n_meas_start,stop=n_meas_end,step=n_meas_step)))
    n_meas_length = length(n_meas_l)
    depth_l = floor.(Int,collect(range(depth_start,stop=depth_end,step=depth_step)))
    depth_length = length(n_meas_l)
    SAB_ave_arr,SAB_std_arr,SBC_ave_arr,SBC_std_arr,SB_ave_arr,SB_std_arr,SABC_ave_arr,SABC_std_arr = scanNerasDepth1D(simBWEras1D,n_Asites,n_Bsites,n_Csites,n_meas_start,n_meas_end,n_meas_step,depth_start,depth_end,depth_step,nsim,true)
    save2DData(n_meas_l,depth_l,SAB_ave_arr,SAB_std_arr,string("data/240114/240114_",save_idx,"_SAB"))
    save2DData(n_meas_l,depth_l,SBC_ave_arr,SBC_std_arr,string("data/240114/240114_",save_idx,"_SBC"))
    save2DData(n_meas_l,depth_l,SB_ave_arr,SB_std_arr,string("data/240114/240114_",save_idx,"_SB"))
    save2DData(n_meas_l,depth_l,SABC_ave_arr,SABC_std_arr,string("data/240114/240114_",save_idx,"_SABC"))
end

# exp_1D erasure
let
    run(`clear`)
    for i = 1:5
        n_Asites = 20
        n_Bsites = 100 + 10*i
        n_Csites = 20
        nsim = 100
        n_meas_start = 30 + 3*i
        n_meas_end = 70 + 7*i
        n_meas_step = 5 + i
        depth_start = 20
        depth_end = 80
        depth_step = 6

        save_idx = 3+i
        n_meas_l = floor.(Int,collect(range(n_meas_start,stop=n_meas_end,step=n_meas_step)))
        n_meas_length = length(n_meas_l)
        depth_l = floor.(Int,collect(range(depth_start,stop=depth_end,step=depth_step)))
        depth_length = length(n_meas_l)
        SAB_ave_arr,SAB_std_arr,SBC_ave_arr,SBC_std_arr,SB_ave_arr,SB_std_arr,SABC_ave_arr,SABC_std_arr = scanNerasDepth1D(simBWEras1D,n_Asites,n_Bsites,n_Csites,n_meas_start,n_meas_end,n_meas_step,depth_start,depth_end,depth_step,nsim,true)
        save2DData(n_meas_l,depth_l,SAB_ave_arr,SAB_std_arr,string("data/240114/240114_",save_idx,"_SAB"))
        save2DData(n_meas_l,depth_l,SBC_ave_arr,SBC_std_arr,string("data/240114/240114_",save_idx,"_SBC"))
        save2DData(n_meas_l,depth_l,SB_ave_arr,SB_std_arr,string("data/240114/240114_",save_idx,"_SB"))
        save2DData(n_meas_l,depth_l,SABC_ave_arr,SABC_std_arr,string("data/240114/240114_",save_idx,"_SABC"))
    end
end

let 
    dx = 20
    dy = 20
    depth = 6
    dA = 1
    dC = 1
    idx_i = 10

    n_Asites = dA*dx
    n_Csites = dC*dx
    dB = dy-dA-dC
    n_Bsites = dB*dx
    reg = genProdInitABC(n_Asites,n_Bsites,n_Csites)
    apply2DBW(reg,dx,dy,depth)

    for i = 1:n_Bsites
        apply!(reg,sMZ(i+n_Asites,1))
    end
    for i = 1:dA
        if i != idx_i
            apply!(reg,sMZ(i,1))
        end
    end

    AsitesArr = collect(idx_i:idx_i)
    CsitesArr = collect(n_Asites+n_Bsites+1:n_Asites+n_Bsites+n_Csites)
    SA = entanglement_entropy(reg.stab,AsitesArr,Val(:rref))
    SC = entanglement_entropy(reg.stab,CsitesArr,Val(:rref))
    SAC = entanglement_entropy(reg.stab,vcat(AsitesArr,CsitesArr),Val(:rref))
    @show SA, SC, SAC
end

# 2D noisy BW
let
    run(`clear`)
    # for i in [10,20,30,40,50,60,70,80,90]
    # for i in [0]
    # for i in [10,20,30,40,50,60]
    for i in [0,1,2,3,4]
        dx = 10
        nsim = 100
        dy_start = 2
        dy_end = 20
        dy_step = 2
        noise = 0.01+0.0025*(i)
        depth = 5

        save_idx_start = 2800
        dy_l = floor.(Int,collect(range(dy_start,stop=dy_end,step=dy_step)))
        dist_ave_arr, dist_std_arr = scandY2D(sim2DNoisyBW,dx,dy_start,dy_end,dy_step,noise,depth,nsim,true)
        save1DData(dy_l,dist_ave_arr, dist_std_arr,string("data/250319/250319_",save_idx_start+i))
    end
end

let 
    for i in [0]
        dx = 10
        nsim = 100
        dy_start = 2
        dy_end = 20
        dy_step = 2
        noise = 0.01*i
        depth = 8
    
        save_idx_start = 300
        dy_l = floor.(Int,collect(range(dy_start,stop=dy_end,step=dy_step)))
        dist_ave_arr, dist_std_arr = scandY2D(sim2DNoisyBW,dx,dy_start,dy_end,dy_step,noise,depth,nsim,true)
        save1DData(dy_l,dist_ave_arr, dist_std_arr,string("data/250319/250319_",save_idx_start+i))
    end
end

let
    s = S"XZZ";
    ms = MixedStabilizer(s)
    new_state, idx, anticommute = project!(ms, P"IIY")
    @show new_state
    @show idx, anticommute
end

let 
    s = S"XXX
        IZZ";
    ms = MixedStabilizer(s)
    sv = stabilizerview(ms)
    p = genIStr(3)
    p[1] = (false, true)
    for gen in sv
        @show comm(gen, p)
    end
end

let 
    p = genIStr(3)
    p[1] = (false, true)
    reg = genProdInitABC(1,1,1)
    sv = stabilizerview(reg.stab)
    for gen in sv
        @show gen[1] == (false, true)
    end
    traceout!(reg,1)
    stabilizerview(reg.stab)
    @show expect(P"XII",reg.stab)==0
    apply!(reg,sMX(3,1))
    stabilizerview(reg.stab)
end