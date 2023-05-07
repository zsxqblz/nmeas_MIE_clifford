function vertexToId2D(dx::Int64,dy::Int64,x::Int64,y::Int64)
    return dx*(y-1) + x
end

function idToVertex2D(dx::Int64,dy::Int64,id::Int64)
    x = id % dx
    y = (id-x)÷dx + 1
    return x,y
end

function apply2DBW(reg::Register,dx::Int64,dy::Int64,depth::Int64)
    n_Bsites = dx*dy
    n_Edge = dx*(dy-1) + (dx-1)*dy
    for i = 1:n_Bsites*depth
        edge_start = rand(1:n_Edge)
        # edge along y
        if edge_start <= dx*(dy-1)
            edge_end_x, edge_end_y = idToVertex2D(dx,dy,edge_start)
            edge_end_y = edge_end_y + 1
            edge_end = vertexToId2D(dx,dy,edge_end_x,edge_end_y)
        # edge along x
        else
            edge_start = edge_start - dx*(dy-1)
            edge_start = edge_start + (edge_start-1)÷(dx-1)
            edge_end_x, edge_end_y = idToVertex2D(dx,dy,edge_start)
            edge_end_y = edge_end_y + 1
            edge_end = vertexToId2D(dx,dy,edge_end_x,edge_end_y)
        end
        apply!(reg.stab,random_clifford(2),[edge_start+1,edge_end+1])
    end
end

function sim2DBW(dx::Int64,dy::Int64,n_meas::Int64,depth::Int64)
    n_Asites = 1
    n_Csites = 1
    n_Bsites = dx*dy
    reg = genInitABC(n_Asites,n_Bsites,n_Csites)
    apply2DBW(reg,dx,dy,depth)
    randMeasB(reg,n_Asites,n_Bsites,n_Csites,n_meas)
    return cmi(reg,n_Asites,n_Bsites,n_Csites)
end

function gen2DExp(sim,dx::Int64,dy::Int64,n_meas::Int64,depth::Int64,nsim::Int64)
    cmi_l = Vector{Int64}(undef,nsim)
    for i = 1:nsim
        cmi_l[i] = sim(dx,dy,n_meas,depth)
    end
    return cmi_l
end

function scanNmeas2D(sim,dx::Int64,dy::Int64,n_meas_start::Int64,n_meas_end::Int64,n_meas_step::Int64,depth::Int64,nsim::Int64,showProg::Bool=false)
    n_meas_l = floor.(Int,collect(range(n_meas_start,stop=n_meas_end,step=n_meas_step)))
    n_meas_length = length(n_meas_l)
    cmi_ave = Vector{Float64}(undef,n_meas_length)
    cmi_std = Vector{Float64}(undef,n_meas_length)
    if showProg
        @showprogress for (i,n_meas) in enumerate(n_meas_l)
            cmi_l = gen2DExp(sim,dx,dy,n_meas,depth,nsim)
            cmi_ave[i] = mean(cmi_l)
            cmi_std[i] = std(cmi_l)
        end
    else
        for (i,n_meas) in enumerate(n_meas_l)
            cmi_l = gen2DExp(sim,dx,dy,n_meas,depth,nsim)
            cmi_ave[i] = mean(cmi_l)
            cmi_std[i] = std(cmi_l)
        end
    end
    return cmi_ave,cmi_std
end


function scanNmeasDepth2D(sim,dx::Int64,dy::Int64,n_meas_start::Int64,n_meas_end::Int64,n_meas_step::Int64,depth_start::Int64,depth_end::Int64,depth_step::Int64,nsim::Int64,showProg::Bool=false)
    n_meas_l = floor.(Int,collect(range(n_meas_start,stop=n_meas_end,step=n_meas_step)))
    n_meas_length = length(n_meas_l)
    depth_l = floor.(Int,collect(range(depth_start,stop=depth_end,step=depth_step)))
    depth_length = length(depth_l)
    cmi_ave_arr = Array{Float64}(undef,n_meas_length,depth_length)
    cmi_std_arr = Array{Float64}(undef,n_meas_length,depth_length)

    if showProg
        @showprogress for (i,depth) in enumerate(depth_l)
            cmi_ave_l, cmi_std_l = scanNmeas2D(sim,dx,dy,n_meas_start,n_meas_end,n_meas_step,depth,nsim)
            cmi_ave_arr[:,i] =  cmi_ave_l 
            cmi_std_arr[:,i] =  cmi_std_l 
        end
    else
        for (i,depth) in enumerate(depth_l)
            cmi_ave_l, cmi_std_l = scanNmeas2D(sim,dx,dy,n_meas_start,n_meas_end,n_meas_step,depth,nsim)
            cmi_ave_arr[:,i] =  cmi_ave_l 
            cmi_std_arr[:,i] =  cmi_std_l 
        end
    end

    return cmi_ave_arr,cmi_std_arr
end