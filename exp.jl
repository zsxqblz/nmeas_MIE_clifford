function expRndClifMeas(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas::Int64,nsim::Int64)
    cmi_l = Vector{Int64}(undef,nsim)
    for i = 1:nsim
        cmi_l[i] = simRndClifMeas(n_Asites,n_Bsites,n_Csites,n_meas)
    end
    return cmi_l
end

function expBWMeas(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas::Int64,depth::Int64,nsim::Int64)
    cmi_l = Vector{Int64}(undef,nsim)
    for i = 1:nsim
        cmi_l[i] = simBWMeas(n_Asites,n_Bsites,n_Csites,n_meas,depth)
    end
    return cmi_l
end

function expHPClifMeas(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas::Int64,nsim::Int64)
    cmi_l = Vector{Int64}(undef,nsim)
    for i = 1:nsim
        cmi_l[i] = simHPClifMeas(n_Asites,n_Bsites,n_Csites,n_meas)
    end
    return cmi_l
end

function expTMClifBMeas(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas::Int64,nsim::Int64)
    cmi_l = Vector{Int64}(undef,nsim)
    for i = 1:nsim
        cmi_l[i] = simTMClifBMeas(n_Asites,n_Bsites,n_Csites,n_meas)
    end
    return cmi_l
end

function expHPBWMeas(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas::Int64,depth::Int64,nsim::Int64)
    cmi_l = Vector{Int64}(undef,nsim)
    for i = 1:nsim
        cmi_l[i] = simHPBWMeas(n_Asites,n_Bsites,n_Csites,n_meas,depth)
    end
    return cmi_l
end

function expHPClifBMeas(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas::Int64,nsim::Int64)
    cmi_l = Vector{Int64}(undef,nsim)
    for i = 1:nsim
        cmi_l[i] = simHPClifBMeas(n_Asites,n_Bsites,n_Csites,n_meas)
    end
    return cmi_l
end

function expHPClifPMeas(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas::Int64,nsim::Int64)
    cmi_l = Vector{Int64}(undef,nsim)
    for i = 1:nsim
        cmi_l[i] = simHPClifPMeas(n_Asites,n_Bsites,n_Csites,n_meas)
    end
    return cmi_l
end

function expHPBWBMeas(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas::Int64,depth::Int64,nsim::Int64)
    cmi_l = Vector{Int64}(undef,nsim)
    for i = 1:nsim
        cmi_l[i] = simHPBWBMeas(n_Asites,n_Bsites,n_Csites,n_meas,depth)
    end
    return cmi_l
end

function expHPDiffBMeas(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas::Int64,depth::Int64,nsim::Int64)
    cmi_l = Vector{Int64}(undef,nsim)
    for i = 1:nsim
        cmi_l[i] = simHPDiffBMeas(n_Asites,n_Bsites,n_Csites,n_meas,depth)
    end
    return cmi_l
end

function expTwoClifMeas(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas::Int64,nsim::Int64)
    cmi_l = Vector{Int64}(undef,nsim)
    for i = 1:nsim
        cmi_l[i] = simTwoClifMeas(n_Asites,n_Bsites,n_Csites,n_meas)
    end
    return cmi_l
end

function expTwoClifBMeas(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas::Int64,nsim::Int64)
    cmi_l = Vector{Int64}(undef,nsim)
    for i = 1:nsim
        cmi_l[i] = simTwoClifBMeas(n_Asites,n_Bsites,n_Csites,n_meas)
    end
    return cmi_l
end

function scanNmeasRndClif(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas_start::Int64,n_meas_end::Int64,n_meas_step::Int64,nsim::Int64,showProg::Bool=false)
    n_meas_l = floor.(Int,collect(range(n_meas_start,stop=n_meas_end,step=n_meas_step)))
    n_meas_length = length(n_meas_l)
    cmi_ave = Vector{Float64}(undef,n_meas_length)
    cmi_std = Vector{Float64}(undef,n_meas_length)
    if showProg
        @showprogress for (i,n_meas) in enumerate(n_meas_l)
            cmi_l = expRndClifMeas(n_Asites,n_Bsites,n_Csites,n_meas,nsim)
            cmi_ave[i] = mean(cmi_l)
            cmi_std[i] = std(cmi_l)
        end
    else
        for (i,n_meas) in enumerate(n_meas_l)
            cmi_l = expRndClifMeas(n_Asites,n_Bsites,n_Csites,n_meas,nsim)
            cmi_ave[i] = mean(cmi_l)
            cmi_std[i] = std(cmi_l)
        end
    end
    return cmi_ave,cmi_std
end

function scanNmeasBW(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas_start::Int64,n_meas_end::Int64,n_meas_step::Int64,depth::Int64,nsim::Int64,showProg::Bool=false)
    n_meas_l = floor.(Int,collect(range(n_meas_start,stop=n_meas_end,step=n_meas_step)))
    n_meas_length = length(n_meas_l)
    cmi_ave = Vector{Float64}(undef,n_meas_length)
    cmi_std = Vector{Float64}(undef,n_meas_length)
    if showProg
        @showprogress for (i,n_meas) in enumerate(n_meas_l)
            cmi_l = expBWMeas(n_Asites,n_Bsites,n_Csites,n_meas,depth,nsim)
            cmi_ave[i] = mean(cmi_l)
            cmi_std[i] = std(cmi_l)
        end
    else
        for (i,n_meas) in enumerate(n_meas_l)
            cmi_l = expBWMeas(n_Asites,n_Bsites,n_Csites,n_meas,depth,nsim)
            cmi_ave[i] = mean(cmi_l)
            cmi_std[i] = std(cmi_l)
        end
    end
    return cmi_ave,cmi_std
end

function scanNmeasHP(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas_start::Int64,n_meas_end::Int64,n_meas_step::Int64,nsim::Int64,showProg::Bool=false)
    n_meas_l = floor.(Int,collect(range(n_meas_start,stop=n_meas_end,step=n_meas_step)))
    n_meas_length = length(n_meas_l)
    cmi_ave = Vector{Float64}(undef,n_meas_length)
    cmi_std = Vector{Float64}(undef,n_meas_length)
    if showProg
        @showprogress for (i,n_meas) in enumerate(n_meas_l)
            cmi_l = expHPClifMeas(n_Asites,n_Bsites,n_Csites,n_meas,nsim)
            cmi_ave[i] = mean(cmi_l)
            cmi_std[i] = std(cmi_l)
        end
    else
        for (i,n_meas) in enumerate(n_meas_l)
            cmi_l = expHPClifMeas(n_Asites,n_Bsites,n_Csites,n_meas,nsim)
            cmi_ave[i] = mean(cmi_l)
            cmi_std[i] = std(cmi_l)
        end
    end
    return cmi_ave,cmi_std
end

function scanNmeasTMB(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas_start::Int64,n_meas_end::Int64,n_meas_step::Int64,nsim::Int64,showProg::Bool=false)
    n_meas_l = floor.(Int,collect(range(n_meas_start,stop=n_meas_end,step=n_meas_step)))
    n_meas_length = length(n_meas_l)
    cmi_ave = Vector{Float64}(undef,n_meas_length)
    cmi_std = Vector{Float64}(undef,n_meas_length)
    if showProg
        @showprogress for (i,n_meas) in enumerate(n_meas_l)
            cmi_l = expTMClifBMeas(n_Asites,n_Bsites,n_Csites,n_meas,nsim)
            cmi_ave[i] = mean(cmi_l)
            cmi_std[i] = std(cmi_l)
        end
    else
        for (i,n_meas) in enumerate(n_meas_l)
            cmi_l = expTMClifBMeas(n_Asites,n_Bsites,n_Csites,n_meas,nsim)
            cmi_ave[i] = mean(cmi_l)
            cmi_std[i] = std(cmi_l)
        end
    end
    return cmi_ave,cmi_std
end

function scanNmeasHPB(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas_start::Int64,n_meas_end::Int64,n_meas_step::Int64,nsim::Int64,showProg::Bool=false)
    n_meas_l = floor.(Int,collect(range(n_meas_start,stop=n_meas_end,step=n_meas_step)))
    n_meas_length = length(n_meas_l)
    cmi_ave = Vector{Float64}(undef,n_meas_length)
    cmi_std = Vector{Float64}(undef,n_meas_length)
    if showProg
        @showprogress for (i,n_meas) in enumerate(n_meas_l)
            cmi_l = expHPClifBMeas(n_Asites,n_Bsites,n_Csites,n_meas,nsim)
            cmi_ave[i] = mean(cmi_l)
            cmi_std[i] = std(cmi_l)
        end
    else
        for (i,n_meas) in enumerate(n_meas_l)
            cmi_l = expHPClifBMeas(n_Asites,n_Bsites,n_Csites,n_meas,nsim)
            cmi_ave[i] = mean(cmi_l)
            cmi_std[i] = std(cmi_l)
        end
    end
    return cmi_ave,cmi_std
end

function scanNmeasHPP(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas_start::Int64,n_meas_end::Int64,n_meas_step::Int64,nsim::Int64,showProg::Bool=false)
    n_meas_l = floor.(Int,collect(range(n_meas_start,stop=n_meas_end,step=n_meas_step)))
    n_meas_length = length(n_meas_l)
    cmi_ave = Vector{Float64}(undef,n_meas_length)
    cmi_std = Vector{Float64}(undef,n_meas_length)
    if showProg
        @showprogress for (i,n_meas) in enumerate(n_meas_l)
            cmi_l = expHPClifPMeas(n_Asites,n_Bsites,n_Csites,n_meas,nsim)
            cmi_ave[i] = mean(cmi_l)
            cmi_std[i] = std(cmi_l)
        end
    else
        for (i,n_meas) in enumerate(n_meas_l)
            cmi_l = expHPClifPMeas(n_Asites,n_Bsites,n_Csites,n_meas,nsim)
            cmi_ave[i] = mean(cmi_l)
            cmi_std[i] = std(cmi_l)
        end
    end
    return cmi_ave,cmi_std
end

function scanNmeasHPBBW(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas_start::Int64,n_meas_end::Int64,n_meas_step::Int64,depth::Int64,nsim::Int64,showProg::Bool=false)
    n_meas_l = floor.(Int,collect(range(n_meas_start,stop=n_meas_end,step=n_meas_step)))
    n_meas_length = length(n_meas_l)
    cmi_ave = Vector{Float64}(undef,n_meas_length)
    cmi_std = Vector{Float64}(undef,n_meas_length)
    if showProg
        @showprogress for (i,n_meas) in enumerate(n_meas_l)
            cmi_l = expHPBWBMeas(n_Asites,n_Bsites,n_Csites,n_meas,depth,nsim)
            cmi_ave[i] = mean(cmi_l)
            cmi_std[i] = std(cmi_l)
        end
    else
        for (i,n_meas) in enumerate(n_meas_l)
            cmi_l = expHPBWBMeas(n_Asites,n_Bsites,n_Csites,n_meas,depth,nsim)
            cmi_ave[i] = mean(cmi_l)
            cmi_std[i] = std(cmi_l)
        end
    end
    return cmi_ave,cmi_std
end

function scanNmeasHPBDiff(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas_start::Int64,n_meas_end::Int64,n_meas_step::Int64,depth::Int64,nsim::Int64,showProg::Bool=false)
    n_meas_l = floor.(Int,collect(range(n_meas_start,stop=n_meas_end,step=n_meas_step)))
    n_meas_length = length(n_meas_l)
    cmi_ave = Vector{Float64}(undef,n_meas_length)
    cmi_std = Vector{Float64}(undef,n_meas_length)
    if showProg
        @showprogress for (i,n_meas) in enumerate(n_meas_l)
            cmi_l = expHPDiffBMeas(n_Asites,n_Bsites,n_Csites,n_meas,depth,nsim)
            cmi_ave[i] = mean(cmi_l)
            cmi_std[i] = std(cmi_l)
        end
    else
        for (i,n_meas) in enumerate(n_meas_l)
            cmi_l = expHPDiffBMeas(n_Asites,n_Bsites,n_Csites,n_meas,depth,nsim)
            cmi_ave[i] = mean(cmi_l)
            cmi_std[i] = std(cmi_l)
        end
    end
    return cmi_ave,cmi_std
end

function scanNmeasHPBW(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas_start::Int64,n_meas_end::Int64,n_meas_step::Int64,depth::Int64,nsim::Int64,showProg::Bool=false)
    n_meas_l = floor.(Int,collect(range(n_meas_start,stop=n_meas_end,step=n_meas_step)))
    n_meas_length = length(n_meas_l)
    cmi_ave = Vector{Float64}(undef,n_meas_length)
    cmi_std = Vector{Float64}(undef,n_meas_length)
    if showProg
        @showprogress for (i,n_meas) in enumerate(n_meas_l)
            cmi_l = expHPBWMeas(n_Asites,n_Bsites,n_Csites,n_meas,depth,nsim)
            cmi_ave[i] = mean(cmi_l)
            cmi_std[i] = std(cmi_l)
        end
    else
        for (i,n_meas) in enumerate(n_meas_l)
            cmi_l = expHPBWMeas(n_Asites,n_Bsites,n_Csites,n_meas,depth,nsim)
            cmi_ave[i] = mean(cmi_l)
            cmi_std[i] = std(cmi_l)
        end
    end
    return cmi_ave,cmi_std
end


function scanNmeasTwo(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas_start::Int64,n_meas_end::Int64,n_meas_step::Int64,nsim::Int64,showProg::Bool=false)
    n_meas_l = floor.(Int,collect(range(n_meas_start,stop=n_meas_end,step=n_meas_step)))
    n_meas_length = length(n_meas_l)
    cmi_ave = Vector{Float64}(undef,n_meas_length)
    cmi_std = Vector{Float64}(undef,n_meas_length)
    if showProg
        @showprogress for (i,n_meas) in enumerate(n_meas_l)
            cmi_l = expTwoClifMeas(n_Asites,n_Bsites,n_Csites,n_meas,nsim)
            cmi_ave[i] = mean(cmi_l)
            cmi_std[i] = std(cmi_l)
        end
    else
        for (i,n_meas) in enumerate(n_meas_l)
            cmi_l = expTwoClifMeas(n_Asites,n_Bsites,n_Csites,n_meas,nsim)
            cmi_ave[i] = mean(cmi_l)
            cmi_std[i] = std(cmi_l)
        end
    end
    return cmi_ave,cmi_std
end

function scanNmeasTwoB(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas_start::Int64,n_meas_end::Int64,n_meas_step::Int64,nsim::Int64,showProg::Bool=false)
    n_meas_l = floor.(Int,collect(range(n_meas_start,stop=n_meas_end,step=n_meas_step)))
    n_meas_length = length(n_meas_l)
    cmi_ave = Vector{Float64}(undef,n_meas_length)
    cmi_std = Vector{Float64}(undef,n_meas_length)
    if showProg
        @showprogress for (i,n_meas) in enumerate(n_meas_l)
            cmi_l = expTwoClifBMeas(n_Asites,n_Bsites,n_Csites,n_meas,nsim)
            cmi_ave[i] = mean(cmi_l)
            cmi_std[i] = std(cmi_l)
        end
    else
        for (i,n_meas) in enumerate(n_meas_l)
            cmi_l = expTwoClifBMeas(n_Asites,n_Bsites,n_Csites,n_meas,nsim)
            cmi_ave[i] = mean(cmi_l)
            cmi_std[i] = std(cmi_l)
        end
    end
    return cmi_ave,cmi_std
end


function save1DData(scan_l,cmi_ave_l,cmi_std_l,file_name)
    df = DataFrame()
    df.scan_l = scan_l
    df.cmi_ave_l = cmi_ave_l
    df.cmi_std_l = cmi_std_l
    CSV.write(file_name, df)
end

function scanNmeasDepthBW(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas_start::Int64,n_meas_end::Int64,n_meas_step::Int64,depth_start::Int64,depth_end::Int64,depth_step::Int64,nsim::Int64,showProg::Bool=false)
    n_meas_l = floor.(Int,collect(range(n_meas_start,stop=n_meas_end,step=n_meas_step)))
    n_meas_length = length(n_meas_l)
    depth_l = floor.(Int,collect(range(depth_start,stop=depth_end,step=depth_step)))
    depth_length = length(depth_l)
    cmi_ave_arr = Array{Float64}(undef,n_meas_length,depth_length)
    cmi_std_arr = Array{Float64}(undef,n_meas_length,depth_length)

    if showProg
        @showprogress for (i,depth) in enumerate(depth_l)
            cmi_ave_l, cmi_std_l = scanNmeasBW(n_Asites,n_Bsites,n_Csites,n_meas_start,n_meas_end,n_meas_step,depth,nsim)
            cmi_ave_arr[:,i] =  cmi_ave_l 
            cmi_std_arr[:,i] =  cmi_std_l 
        end
    else
        for (i,depth) in enumerate(depth_l)
            cmi_ave_l, cmi_std_l = scanNmeasBW(n_Asites,n_Bsites,n_Csites,n_meas_start,n_meas_end,n_meas_step,depth,nsim)
            cmi_ave_arr[:,i] =  cmi_ave_l 
            cmi_std_arr[:,i] =  cmi_std_l 
        end
    end

    return cmi_ave_arr,cmi_std_arr
end

function scanNmeasDepthHPBBW(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas_start::Int64,n_meas_end::Int64,n_meas_step::Int64,depth_start::Int64,depth_end::Int64,depth_step::Int64,nsim::Int64,showProg::Bool=false)
    n_meas_l = floor.(Int,collect(range(n_meas_start,stop=n_meas_end,step=n_meas_step)))
    n_meas_length = length(n_meas_l)
    depth_l = floor.(Int,collect(range(depth_start,stop=depth_end,step=depth_step)))
    depth_length = length(depth_l)
    cmi_ave_arr = Array{Float64}(undef,n_meas_length,depth_length)
    cmi_std_arr = Array{Float64}(undef,n_meas_length,depth_length)

    if showProg
        @showprogress for (i,depth) in enumerate(depth_l)
            cmi_ave_l, cmi_std_l = scanNmeasHPBBW(n_Asites,n_Bsites,n_Csites,n_meas_start,n_meas_end,n_meas_step,depth,nsim)
            cmi_ave_arr[:,i] =  cmi_ave_l 
            cmi_std_arr[:,i] =  cmi_std_l 
        end
    else
        for (i,depth) in enumerate(depth_l)
            cmi_ave_l, cmi_std_l = scanNmeasHPBBW(n_Asites,n_Bsites,n_Csites,n_meas_start,n_meas_end,n_meas_step,depth,nsim)
            cmi_ave_arr[:,i] =  cmi_ave_l 
            cmi_std_arr[:,i] =  cmi_std_l 
        end
    end

    return cmi_ave_arr,cmi_std_arr
end

function scanNmeasDepthHPBDiff(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas_start::Int64,n_meas_end::Int64,n_meas_step::Int64,depth_start::Int64,depth_end::Int64,depth_step::Int64,nsim::Int64,showProg::Bool=false)
    n_meas_l = floor.(Int,collect(range(n_meas_start,stop=n_meas_end,step=n_meas_step)))
    n_meas_length = length(n_meas_l)
    depth_l = floor.(Int,collect(range(depth_start,stop=depth_end,step=depth_step)))
    depth_length = length(depth_l)
    cmi_ave_arr = Array{Float64}(undef,n_meas_length,depth_length)
    cmi_std_arr = Array{Float64}(undef,n_meas_length,depth_length)

    if showProg
        @showprogress for (i,depth) in enumerate(depth_l)
            cmi_ave_l, cmi_std_l = scanNmeasHPBDiff(n_Asites,n_Bsites,n_Csites,n_meas_start,n_meas_end,n_meas_step,depth,nsim)
            cmi_ave_arr[:,i] =  cmi_ave_l 
            cmi_std_arr[:,i] =  cmi_std_l 
        end
    else
        for (i,depth) in enumerate(depth_l)
            cmi_ave_l, cmi_std_l = scanNmeasHPBDiff(n_Asites,n_Bsites,n_Csites,n_meas_start,n_meas_end,n_meas_step,depth,nsim)
            cmi_ave_arr[:,i] =  cmi_ave_l 
            cmi_std_arr[:,i] =  cmi_std_l 
        end
    end

    return cmi_ave_arr,cmi_std_arr
end

function scanNmeasDepthHPBW(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas_start::Int64,n_meas_end::Int64,n_meas_step::Int64,depth_start::Int64,depth_end::Int64,depth_step::Int64,nsim::Int64,showProg::Bool=false)
    n_meas_l = floor.(Int,collect(range(n_meas_start,stop=n_meas_end,step=n_meas_step)))
    n_meas_length = length(n_meas_l)
    depth_l = floor.(Int,collect(range(depth_start,stop=depth_end,step=depth_step)))
    depth_length = length(depth_l)
    cmi_ave_arr = Array{Float64}(undef,n_meas_length,depth_length)
    cmi_std_arr = Array{Float64}(undef,n_meas_length,depth_length)

    if showProg
        @showprogress for (i,depth) in enumerate(depth_l)
            cmi_ave_l, cmi_std_l = scanNmeasHPBW(n_Asites,n_Bsites,n_Csites,n_meas_start,n_meas_end,n_meas_step,depth,nsim)
            cmi_ave_arr[:,i] =  cmi_ave_l 
            cmi_std_arr[:,i] =  cmi_std_l 
        end
    else
        for (i,depth) in enumerate(depth_l)
            cmi_ave_l, cmi_std_l = scanNmeasHPBW(n_Asites,n_Bsites,n_Csites,n_meas_start,n_meas_end,n_meas_step,depth,nsim)
            cmi_ave_arr[:,i] =  cmi_ave_l 
            cmi_std_arr[:,i] =  cmi_std_l 
        end
    end

    return cmi_ave_arr,cmi_std_arr
end

function save2DData(scanx_l,scany_l,cmi_ave_arr,cmi_std_arr,file_name)
    df_scanx = DataFrame()
    df_scany = DataFrame()
    df_data = DataFrame()
    df_scanx.scanx_l = scanx_l
    df_scany.scany_l = scany_l
    df_data.cmi_ave_l = collect(Iterators.flatten(cmi_ave_arr)) 
    df_data.cmi_std_l = collect(Iterators.flatten(cmi_std_arr))

    CSV.write(file_name*"_scanx.csv", df_scanx)
    CSV.write(file_name*"_scany.csv", df_scany)
    CSV.write(file_name*"_data.csv", df_data)
end