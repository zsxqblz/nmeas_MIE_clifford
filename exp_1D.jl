function simBWMeas1D(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas::Int64,depth::Int64)
    # choi
    # reg = genInitABC(n_Asites,n_Bsites,n_Csites)
    # applyBW(reg,n_Asites,n_Bsites,n_Csites,depth)
    # no choi
    reg = genProdInitABC(n_Asites,n_Bsites,n_Csites)
    applyProdBW(reg,n_Asites+n_Bsites+n_Csites,depth)

    randMeasB(reg,n_Asites,n_Bsites,n_Csites,n_meas)
    AsitesArr = collect(1:n_Asites)
    CsitesArr = collect(n_Asites+n_Bsites+1:n_Asites+n_Bsites+n_Csites)
    SA = entanglement_entropy(reg.stab,AsitesArr,Val(:rref))
    SC = entanglement_entropy(reg.stab,CsitesArr,Val(:rref))
    SAC = entanglement_entropy(reg.stab,vcat(AsitesArr,CsitesArr),Val(:rref))
    return SA, SC, SAC
    # return cmi(reg,n_Asites,n_Bsites,n_Csites)
end

function simBWEras1D(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas::Int64,depth::Int64)
    # choi
    # reg = genInitABC(n_Asites,n_Bsites,n_Csites)
    # applyBW(reg,n_Asites,n_Bsites,n_Csites,depth)
    # no choi
    reg = genProdInitABC(n_Asites,n_Bsites,n_Csites)
    applyProdBW(reg,n_Asites+n_Bsites+n_Csites,depth)

    randErasB(reg,n_Asites,n_Bsites,n_Csites,n_meas)
    AsitesArr = collect(1:n_Asites)
    CsitesArr = collect(n_Asites+n_Bsites+1:n_Asites+n_Bsites+n_Csites)
    SAB = entanglement_entropy(reg.stab,collect(1:n_Asites+n_Bsites),Val(:rref))
    SBC = entanglement_entropy(reg.stab,collect(n_Asites+1:n_Asites+n_Bsites+n_Csites),Val(:rref))
    SB = entanglement_entropy(reg.stab,collect(n_Asites+1:n_Asites+n_Bsites),Val(:rref))
    SABC = entanglement_entropy(reg.stab,collect(1:n_Asites+n_Bsites+n_Csites),Val(:rref))
    return SAB, SBC, SB, SABC
    # return cmi(reg,n_Asites,n_Bsites,n_Csites)
end

function genExp1D(sim,n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas::Int64,depth::Int64,nsim::Int64)
    # cmi_l = Vector{Int64}(undef,nsim)
    # for i = 1:nsim
    #     cmi_l[i] = sim(n_Asites,n_Bsites,n_Csites,n_meas,depth)
    # end
    # return cmi_l
    SA_l = Vector{Int64}(undef,nsim)
    SC_l = Vector{Int64}(undef,nsim)
    SAC_l = Vector{Int64}(undef,nsim)
    for i = 1:nsim
        SA, SC, SAC = sim(n_Asites,n_Bsites,n_Csites,n_meas,depth)
        SA_l[i] = SA
        SC_l[i] = SC
        SAC_l[i] = SAC
    end
    return SA_l, SC_l, SAC_l
end

function genErasExp1D(sim,n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas::Int64,depth::Int64,nsim::Int64)
    # cmi_l = Vector{Int64}(undef,nsim)
    # for i = 1:nsim
    #     cmi_l[i] = sim(n_Asites,n_Bsites,n_Csites,n_meas,depth)
    # end
    # return cmi_l
    SAB_l = Vector{Int64}(undef,nsim)
    SBC_l = Vector{Int64}(undef,nsim)
    SB_l = Vector{Int64}(undef,nsim)
    SABC_l = Vector{Int64}(undef,nsim)
    for i = 1:nsim
        SAB, SBC, SB, SABC = sim(n_Asites,n_Bsites,n_Csites,n_meas,depth)
        SAB_l[i] = SAB
        SBC_l[i] = SBC
        SB_l[i] = SB
        SABC_l[i] = SABC
    end
    return SAB_l, SBC_l, SB_l, SABC_l
end

function scanNmeas1D(sim,n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas_start::Int64,n_meas_end::Int64,n_meas_step::Int64,depth::Int64,nsim::Int64,showProg::Bool=false)
    n_meas_l = floor.(Int,collect(range(n_meas_start,stop=n_meas_end,step=n_meas_step)))
    n_meas_length = length(n_meas_l)
    # cmi_ave = Vector{Float64}(undef,n_meas_length)
    # cmi_std = Vector{Float64}(undef,n_meas_length)
    SA_ave = Vector{Float64}(undef,n_meas_length)
    SA_std = Vector{Float64}(undef,n_meas_length)
    SC_ave = Vector{Float64}(undef,n_meas_length)
    SC_std = Vector{Float64}(undef,n_meas_length)
    SAC_ave = Vector{Float64}(undef,n_meas_length)
    SAC_std = Vector{Float64}(undef,n_meas_length)
    if showProg
        @showprogress for (i,n_meas) in enumerate(n_meas_l)
            # cmi_l = genExp(sim,n_Asites,n_Bsites,n_Csites,n_meas,depth,nsim)
            # cmi_ave[i] = mean(cmi_l)
            # cmi_std[i] = std(cmi_l)
            SA_l, SC_l, SAC_l = genExp1D(sim,n_Asites,n_Bsites,n_Csites,n_meas,depth,nsim)
            SA_ave[i] = mean(SA_l)
            SA_std[i] = std(SA_l)
            SC_ave[i] = mean(SC_l)
            SC_std[i] = std(SC_l)
            SAC_ave[i] = mean(SAC_l)
            SAC_std[i] = std(SAC_l)
        end
    else
        for (i,n_meas) in enumerate(n_meas_l)
            # cmi_l = genExp(sim,n_Asites,n_Bsites,n_Csites,n_meas,depth,nsim)
            # cmi_ave[i] = mean(cmi_l)
            # cmi_std[i] = std(cmi_l)
            SA_l, SC_l, SAC_l = genExp1D(sim,n_Asites,n_Bsites,n_Csites,n_meas,depth,nsim)
            SA_ave[i] = mean(SA_l)
            SA_std[i] = std(SA_l)
            SC_ave[i] = mean(SC_l)
            SC_std[i] = std(SC_l)
            SAC_ave[i] = mean(SAC_l)
            SAC_std[i] = std(SAC_l)
        end
    end
    # return cmi_ave,cmi_std
    return SA_ave,SA_std,SC_ave,SC_std,SAC_ave,SAC_std
end

function scanNeras1D(sim,n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas_start::Int64,n_meas_end::Int64,n_meas_step::Int64,depth::Int64,nsim::Int64,showProg::Bool=false)
    n_meas_l = floor.(Int,collect(range(n_meas_start,stop=n_meas_end,step=n_meas_step)))
    n_meas_length = length(n_meas_l)
    # cmi_ave = Vector{Float64}(undef,n_meas_length)
    # cmi_std = Vector{Float64}(undef,n_meas_length)
    SAB_ave = Vector{Float64}(undef,n_meas_length)
    SAB_std = Vector{Float64}(undef,n_meas_length)
    SBC_ave = Vector{Float64}(undef,n_meas_length)
    SBC_std = Vector{Float64}(undef,n_meas_length)
    SB_ave = Vector{Float64}(undef,n_meas_length)
    SB_std = Vector{Float64}(undef,n_meas_length)
    SABC_ave = Vector{Float64}(undef,n_meas_length)
    SABC_std = Vector{Float64}(undef,n_meas_length)
    if showProg
        @showprogress for (i,n_meas) in enumerate(n_meas_l)
            # cmi_l = genExp(sim,n_Asites,n_Bsites,n_Csites,n_meas,depth,nsim)
            # cmi_ave[i] = mean(cmi_l)
            # cmi_std[i] = std(cmi_l)
            SAB_l, SBC_l, SB_l, SABC_l = genErasExp1D(sim,n_Asites,n_Bsites,n_Csites,n_meas,depth,nsim)
            SAB_ave[i] = mean(SAB_l)
            SAB_std[i] = std(SAB_l)
            SBC_ave[i] = mean(SBC_l)
            SBC_std[i] = std(SBC_l)
            SB_ave[i] = mean(SB_l)
            SB_std[i] = std(SB_l)
            SABC_ave[i] = mean(SABC_l)
            SABC_std[i] = std(SABC_l)
        end
    else
        for (i,n_meas) in enumerate(n_meas_l)
            # cmi_l = genExp(sim,n_Asites,n_Bsites,n_Csites,n_meas,depth,nsim)
            # cmi_ave[i] = mean(cmi_l)
            # cmi_std[i] = std(cmi_l)
            SAB_l, SBC_l, SB_l, SABC_l = genErasExp1D(sim,n_Asites,n_Bsites,n_Csites,n_meas,depth,nsim)
            SAB_ave[i] = mean(SAB_l)
            SAB_std[i] = std(SAB_l)
            SBC_ave[i] = mean(SBC_l)
            SBC_std[i] = std(SBC_l)
            SB_ave[i] = mean(SB_l)
            SB_std[i] = std(SB_l)
            SABC_ave[i] = mean(SABC_l)
            SABC_std[i] = std(SABC_l)
        end
    end
    # return cmi_ave,cmi_std
    return SAB_ave,SAB_std,SBC_ave,SBC_std,SB_ave,SB_std,SABC_ave,SABC_std
end

function scanNmeasDepth1D(sim,n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas_start::Int64,n_meas_end::Int64,n_meas_step::Int64,depth_start::Int64,depth_end::Int64,depth_step::Int64,nsim::Int64,showProg::Bool=false)
    n_meas_l = floor.(Int,collect(range(n_meas_start,stop=n_meas_end,step=n_meas_step)))
    n_meas_length = length(n_meas_l)
    depth_l = floor.(Int,collect(range(depth_start,stop=depth_end,step=depth_step)))
    depth_length = length(depth_l)
    # cmi_ave_arr = Array{Float64}(undef,n_meas_length,depth_length)
    # cmi_std_arr = Array{Float64}(undef,n_meas_length,depth_length)
    SA_ave_arr = Array{Float64}(undef,n_meas_length,depth_length)
    SA_std_arr = Array{Float64}(undef,n_meas_length,depth_length)
    SC_ave_arr = Array{Float64}(undef,n_meas_length,depth_length)
    SC_std_arr = Array{Float64}(undef,n_meas_length,depth_length)
    SAC_ave_arr = Array{Float64}(undef,n_meas_length,depth_length)
    SAC_std_arr = Array{Float64}(undef,n_meas_length,depth_length)

    if showProg
        @showprogress for (i,depth) in enumerate(depth_l)
            # cmi_ave_l, cmi_std_l = scanNmeas(sim,n_Asites,n_Bsites,n_Csites,n_meas_start,n_meas_end,n_meas_step,depth,nsim)
            # cmi_ave_arr[:,i] =  cmi_ave_l 
            # cmi_std_arr[:,i] =  cmi_std_l 
            SA_ave_l,SA_std_l,SC_ave_l,SC_std_l,SAC_ave_l,SAC_std_l = scanNmeas1D(sim,n_Asites,n_Bsites,n_Csites,n_meas_start,n_meas_end,n_meas_step,depth,nsim)
            SA_ave_arr[:,i] =  SA_ave_l 
            SA_std_arr[:,i] =  SA_std_l
            SC_ave_arr[:,i] =  SC_ave_l 
            SC_std_arr[:,i] =  SC_std_l 
            SAC_ave_arr[:,i] =  SAC_ave_l 
            SAC_std_arr[:,i] =  SAC_std_l 
        end
    else
        for (i,depth) in enumerate(depth_l)
            # cmi_ave_l, cmi_std_l = scanNmeas(sim,n_Asites,n_Bsites,n_Csites,n_meas_start,n_meas_end,n_meas_step,depth,nsim)
            # cmi_ave_arr[:,i] =  cmi_ave_l 
            # cmi_std_arr[:,i] =  cmi_std_l 
            SA_ave_l,SA_std_l,SC_ave_l,SC_std_l,SAC_ave_l,SAC_std_l = scanNmeas1D(sim,n_Asites,n_Bsites,n_Csites,n_meas_start,n_meas_end,n_meas_step,depth,nsim)
            SA_ave_arr[:,i] =  SA_ave_l 
            SA_std_arr[:,i] =  SA_std_l
            SC_ave_arr[:,i] =  SC_ave_l 
            SC_std_arr[:,i] =  SC_std_l 
            SAC_ave_arr[:,i] =  SAC_ave_l 
            SAC_std_arr[:,i] =  SAC_std_l 
        end
    end

    return SA_ave_arr,SA_std_arr,SC_ave_arr,SC_std_arr,SAC_ave_arr,SAC_std_arr
end

function scanNerasDepth1D(sim,n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas_start::Int64,n_meas_end::Int64,n_meas_step::Int64,depth_start::Int64,depth_end::Int64,depth_step::Int64,nsim::Int64,showProg::Bool=false)
    n_meas_l = floor.(Int,collect(range(n_meas_start,stop=n_meas_end,step=n_meas_step)))
    n_meas_length = length(n_meas_l)
    depth_l = floor.(Int,collect(range(depth_start,stop=depth_end,step=depth_step)))
    depth_length = length(depth_l)
    # cmi_ave_arr = Array{Float64}(undef,n_meas_length,depth_length)
    # cmi_std_arr = Array{Float64}(undef,n_meas_length,depth_length)
    SAB_ave_arr = Array{Float64}(undef,n_meas_length,depth_length)
    SAB_std_arr = Array{Float64}(undef,n_meas_length,depth_length)
    SBC_ave_arr = Array{Float64}(undef,n_meas_length,depth_length)
    SBC_std_arr = Array{Float64}(undef,n_meas_length,depth_length)
    SB_ave_arr = Array{Float64}(undef,n_meas_length,depth_length)
    SB_std_arr = Array{Float64}(undef,n_meas_length,depth_length)
    SABC_ave_arr = Array{Float64}(undef,n_meas_length,depth_length)
    SABC_std_arr = Array{Float64}(undef,n_meas_length,depth_length)

    if showProg
        @showprogress for (i,depth) in enumerate(depth_l)
            # cmi_ave_l, cmi_std_l = scanNmeas(sim,n_Asites,n_Bsites,n_Csites,n_meas_start,n_meas_end,n_meas_step,depth,nsim)
            # cmi_ave_arr[:,i] =  cmi_ave_l 
            # cmi_std_arr[:,i] =  cmi_std_l 
            SAB_ave_l,SAB_std_l,SBC_ave_l,SBC_std_l,SB_ave_l,SB_std_l,SABC_ave_l,SABC_std_l = scanNeras1D(sim,n_Asites,n_Bsites,n_Csites,n_meas_start,n_meas_end,n_meas_step,depth,nsim)
            SAB_ave_arr[:,i] =  SAB_ave_l 
            SAB_std_arr[:,i] =  SAB_std_l
            SBC_ave_arr[:,i] =  SBC_ave_l 
            SBC_std_arr[:,i] =  SBC_std_l 
            SB_ave_arr[:,i] =  SB_ave_l 
            SB_std_arr[:,i] =  SB_std_l 
            SABC_ave_arr[:,i] =  SABC_ave_l 
            SABC_std_arr[:,i] =  SABC_std_l 
        end
    else
        for (i,depth) in enumerate(depth_l)
            # cmi_ave_l, cmi_std_l = scanNmeas(sim,n_Asites,n_Bsites,n_Csites,n_meas_start,n_meas_end,n_meas_step,depth,nsim)
            # cmi_ave_arr[:,i] =  cmi_ave_l 
            # cmi_std_arr[:,i] =  cmi_std_l 
            SAB_ave_l,SAB_std_l,SBC_ave_l,SBC_std_l,SB_ave_l,SB_std_l,SABC_ave_l,SABC_std_l = scanNeras1D(sim,n_Asites,n_Bsites,n_Csites,n_meas_start,n_meas_end,n_meas_step,depth,nsim)
            SAB_ave_arr[:,i] =  SAB_ave_l 
            SAB_std_arr[:,i] =  SAB_std_l
            SBC_ave_arr[:,i] =  SBC_ave_l 
            SBC_std_arr[:,i] =  SBC_std_l 
            SB_ave_arr[:,i] =  SB_ave_l 
            SB_std_arr[:,i] =  SB_std_l 
            SABC_ave_arr[:,i] =  SABC_ave_l 
            SABC_std_arr[:,i] =  SABC_std_l 
        end
    end

    return SAB_ave_arr,SAB_std_arr,SBC_ave_arr,SBC_std_arr,SB_ave_arr,SB_std_arr,SABC_ave_arr,SABC_std_arr
end