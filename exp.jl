function expRndClifMeas(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas::Int64,nsim::Int64)
    cmi_l = Vector{Int64}(undef,nsim)
    for i = 1:nsim
        cmi_l[i] = simRndClifMeas(n_Asites,n_Bsites,n_Csites,n_meas)
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

function expHPClifBMeas(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas::Int64,nsim::Int64)
    cmi_l = Vector{Int64}(undef,nsim)
    for i = 1:nsim
        cmi_l[i] = simHPClifBMeas(n_Asites,n_Bsites,n_Csites,n_meas)
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
