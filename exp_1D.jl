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