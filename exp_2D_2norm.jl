function vertexToId2D(dx::Int64,dy::Int64,x::Int64,y::Int64)
    return dx*y + x + 1
end

function idToVertex2D(dx::Int64,dy::Int64,id::Int64)
    x = (id-1) % dx
    y = (id-1-x) ÷ dx
    return x,y
end

function apply2DNoisyBW(reg::Register,dx::Int64,dy::Int64,dA::Int64, dC::Int64,depth::Int64,noise::Float64)
    n_Bsites = dx*dy
    n_Edge = dx*(dy-1) + (dx-1)*dy
    # for i = 1:n_Bsites*depth
    #     edge_start = rand(1:n_Edge)
    for i = 1:depth
        for edge_start = 1:n_Edge
            edge_start_x, edge_start_y = idToVertex2D(dx,dy,edge_start)
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
                edge_end_x = edge_end_x + 1
                edge_end = vertexToId2D(dx,dy,edge_end_x,edge_end_y)
            end
            apply!(reg.stab,random_clifford(2),[edge_start,edge_end])

            rnd = rand()
            if rnd < noise && edge_start_y > i-1 && edge_start_y <= dy-i+1
                traceout!(reg, edge_start)
            end

            rnd = rand()
            if rnd < noise && edge_end_y > i-1 && edge_end_y <= dy-i+1
                traceout!(reg, edge_end)
            end
        end
    end
end

function sim2DNoisyBW(dx::Int64,dy::Int64,dA::Int64,dC::Int64,noise::Float64,depth::Int64)
    n_Asites = dA*dx
    n_Csites = dC*dx
    dB = dy-dA-dC
    n_Bsites = dB*dx
    reg = genProdInitABC(n_Asites,n_Bsites,n_Csites)
    apply2DNoisyBW(reg,dx,dy,dA,dC,depth,noise)
    # num_anticommute = randMeasBProb(reg,n_Asites,n_Bsites,n_Csites,n_Bsites)
    num_anticommute = randMeasBCProb(reg,n_Asites,n_Bsites,n_Csites,dx,dy,dC)

    AsitesArr = collect(1:n_Asites)
    # BsitesArr = collect(n_Asites+1:n_Asites+n_Bsites)
    # CsitesArr = collect(n_Asites+n_Bsites+1:n_Asites+n_Bsites+n_Csites)
    CsiteID = vertexToId2D(dx, dy, Int(floor(dx/2)), dy-dC)
    CsitesArr = [CsiteID]

    SA = entanglement_entropy(reg.stab,AsitesArr,Val(:rref))
    SC = entanglement_entropy(reg.stab,CsitesArr,Val(:rref))
    SAC = entanglement_entropy(reg.stab,vcat(AsitesArr,CsitesArr),Val(:rref))
    
    # prob_sq = 2.0^(-2*num_anticommute)
    # dist_sq = 2.0^(-SAC) - 2.0^(-n_Asites-n_Csites)
    # dist_sq = 2.0^(-SAC)

    # return 2.0^(-2*num_anticommute + dx*(dA+dC) + 2*dx*dB) * (2.0^(-SAC) - 2.0^(-SA-SC))

    # return 2.0^(-2*num_anticommute + dx*dA+1 + 2*dx*dB) * (2.0^(-SAC) - 2.0^(-SA-1))
    return 2.0^(-2*num_anticommute + dx*dA+1 + 2*(dx*(dB+dC)-1)) * (2.0^(-SAC) - 2.0^(-SA-1))

    # return prob_sq*dist_sq
    # return prob
    # return dist_sq
end

function gen2DExp(sim,dx::Int64,dy::Int64,dA::Int64,dC::Int64,noise::Float64,depth::Int64,nsim::Int64)
    # cmi_l = Vector{Int64}(undef,nsim)
    # for i = 1:nsim
    #     cmi_l[i] = sim(dx,dy,n_meas,depth)
    # end
    # return cmi_l
    dist_l = Vector{Float64}(undef,nsim)
    for i = 1:nsim
        dist = sim(dx,dy,dA,dC,noise,depth)
        dist_l[i] = dist
    end
    return dist_l
end



function scandY2D(sim,dx::Int64,dy_start::Int64,dy_end::Int64,dy_step::Int64,noise::Float64,depth::Int64,nsim::Int64,showProg::Bool=false)
    dA = depth
    dC = depth
    # dA = 1
    # dC = 1

    dy_l = floor.(Int,collect(range(dy_start+2*depth,stop=dy_end+2*depth,step=dy_step)))
    dy_length = length(dy_l)
    dist_ave_arr = Array{Float64}(undef,dy_length)
    dist_std_arr = Array{Float64}(undef,dy_length)

    if showProg
        @showprogress for (i,dy) in enumerate(dy_l)
            dist_l = gen2DExp(sim,dx,dy,dA,dC,noise,depth,nsim)
            dist_ave_arr[i] = mean(dist_l)
            dist_std_arr[i] = std(dist_l)
        end
    else
        for (i,depth) in enumerate(depth_l)
            dist_l = gen2DExp(sim,dx,dy,dA,dC,noise,depth,nsim)
            dist_ave_arr[i] = mean(dist_l)
            dist_std_arr[i] = std(dist_l)
        end
    end

    return dist_ave_arr, dist_std_arr
end