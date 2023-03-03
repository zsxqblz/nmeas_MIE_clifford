function genIStr(nsites::Int64)
    p = P"I"
    for i in 2:nsites
        p = p ⊗ P"I"
    end
    return p
end

function genIOp(nsites::Int64)
    p = zero(CliffordOperator,nsites)
    for i in 1:nsites
        p[i,i] = (true,false)
        p[i+nsites,i] = (false,true)
    end
    return p
end

function genBWClif(nsites::Int64,depth::Int64)
    U = genIOp(nsites)
    for i = 1:depth
        parity = i%2
        nbrick = floor(Int,(nsites-parity)/2)
        for j = 1:nbrick
            if j == 1
                UBW = random_clifford(2)
            else
                UBW = UBW ⊗ random_clifford(2)
            end
        end

        if parity == 1
            UBW = genIOp(1) ⊗ UBW
        end

        if (nsites-parity)%2 == 1
            UBW = UBW ⊗ genIOp(1) 
        end

        U = U * UBW 
    end
    return U
end

function reverseOp(UL::CliffordOperator,nsites::Int64)
    UR = copy(UL)
    for i = 1:nsites
        for j = 1:nsites
            UR[i,j] = UL[nsites+1-i,nsites+1-j]
            UR[i+nsites,j] = UL[2*nsites+1-i,nsites+1-j]
        end
    end
    return UR
end

function genBellState(nsites::Int64)
    s = zero(Stabilizer,2*nsites)
    for i = 1:nsites
        s[i,i] = (true,false)
        s[i,i+nsites] = (true,false)
        s[i+nsites,i] = (false,true)
        s[i+nsites,i+nsites] = (false,true)
    end
    return s
end

function genSBellState(nsites::Int64)
    s = zero(Stabilizer,2*nsites)
    for i = 1:nsites
        s[i,i] = (true,false)
        s[i,2*nsites+1-i] = (true,false)
        s[i+nsites,i] = (false,true)
        s[i+nsites,2*nsites+1-i] = (false,true)
    end
    return s
end

function genProdZState(nsites::Int64)
    s = zero(Stabilizer,nsites)
    for i = 1:nsites
        s[i,i] = (false,true)
    end
    return s
end

function genInitABC(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64)
    sA = genBellState(n_Asites)
    sC = genBellState(n_Csites)
    sB = genProdZState(n_Bsites-n_Asites-n_Csites)
    boolStr = Vector{Bool}(undef,n_Asites+n_Bsites+n_Csites)
    reg = Register(MixedDestabilizer(sA⊗sB⊗sC), boolStr)
    return reg
end

function genInitHP(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64)
    sA = genBellState(n_Asites)
    sC = genBellState(n_Csites)
    sB = genSBellState(floor(Int,(n_Bsites-n_Asites-n_Csites)/2))
    boolStr = Vector{Bool}(undef,n_Asites+n_Bsites+n_Csites)
    reg = Register(MixedDestabilizer(sA⊗sB⊗sC), boolStr)
    return reg
end

function applyRandClif(reg::Register,n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64)
    IA = genIOp(n_Asites)
    IC = genIOp(n_Csites)
    UB = random_clifford(n_Bsites)
    U = IA ⊗ UB ⊗ IC
    apply!(reg,U)
end

function applyBW(reg::Register,n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,depth::Int64)
    IA = genIOp(n_Asites)
    IC = genIOp(n_Csites)
    UB = genBWClif(n_Bsites,depth)
    U = IA ⊗ UB ⊗ IC
    apply!(reg,U)
end

function applyHPClif(reg::Register,n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64)
    half_Bsites = floor(Int,n_Bsites/2)
    IA = genIOp(n_Asites)
    IC = genIOp(n_Csites)
    UBL = random_clifford(half_Bsites)
    UBR = reverseOp(UBL,half_Bsites)
    U = IA ⊗ UBL ⊗ UBR ⊗ IC
    apply!(reg,U)
end

function applyHPBW(reg::Register,n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,depth::Int64)
    half_Bsites = floor(Int,n_Bsites/2)
    IA = genIOp(n_Asites)
    IC = genIOp(n_Csites)
    UBL = genBWClif(half_Bsites,depth)
    UBR = reverseOp(UBL,half_Bsites)
    U = IA ⊗ UBL ⊗ UBR ⊗ IC
    apply!(reg,U)
end

function applyTwoRandClif(reg::Register,n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64)
    half_Bsites = floor(Int,n_Bsites/2)
    IA = genIOp(n_Asites)
    IC = genIOp(n_Csites)
    UBL = random_clifford(half_Bsites)
    UBR = random_clifford(half_Bsites)
    U = IA ⊗ UBL ⊗ UBR ⊗ IC
    apply!(reg,U)
end

function randMeasB(reg::Register,n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,nmeas::Int64)
    n_sites = n_Asites+n_Bsites+n_Csites
    meas_idx = sample(collect(1:n_Bsites),nmeas,replace=false)
    for i in meas_idx
        apply!(reg,sMZ(i+n_Asites,1))
    end
end

function randBellMeasB(reg::Register,n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,nmeas::Int64)
    n_sites = n_Asites+n_Bsites+n_Csites
    h_Bsites = floor(Int,n_Bsites/2)
    hmeas = floor(Int,nmeas/2)
    meas_idx = sample(collect(1:h_Bsites),hmeas,replace=false)
    for i in meas_idx
        zz_op = single_z(n_sites,i+n_Asites)
        zz_op[n_Asites+n_Bsites+1-i] = (false,true)
        projectrand!(reg,zz_op)
    end
end

# function applyRandPauli(reg::Register,n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,pauli_prob::Float64)
#     half_Bsites = floor(Int,n_Bsites/2)
#     IA = genIOp(n_Asites)
#     IC = genIOp(n_Csites)
#     IBL = genIOp(half_Bsites)
#     # UBR = random_clifford(half_Bsites)
#     UBR = copy(IBL)
#     pauli_rand = rand(Float64,half_Bsites)
#     for i = 1:half_Bsites
#         if pauli_rand[i] < pauli_prob
#             pauli = mod(rand(Int64),3)
#             if pauli == 0
#                 UBR[i,i]
#         end
#     end
#     U = IA ⊗ IBL ⊗ UBR ⊗ IC
#     apply!(reg,U)
# end

function cmi(reg::Register,n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64)
    SAB = entanglement_entropy(reg.stab,collect(1:n_Asites+n_Bsites),Val(:rref))
    SBC = entanglement_entropy(reg.stab,collect(n_Asites+1:n_Asites+n_Bsites+n_Csites),Val(:rref))
    SB = entanglement_entropy(reg.stab,collect(n_Asites+1:n_Asites+n_Bsites),Val(:rref))
    SABC = entanglement_entropy(reg.stab,collect(1:n_Asites+n_Bsites+n_Csites),Val(:rref))
    return SAB + SBC - SB - SABC 
end

function simRndClifMeas(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas::Int64)
    reg = genInitABC(n_Asites,n_Bsites,n_Csites)
    applyRandClif(reg,n_Asites,n_Bsites,n_Csites)
    randMeasB(reg,n_Asites,n_Bsites,n_Csites,n_meas)
    return cmi(reg,n_Asites,n_Bsites,n_Csites)
end

function simBWMeas(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas::Int64,depth::Int64)
    reg = genInitABC(n_Asites,n_Bsites,n_Csites)
    applyBW(reg,n_Asites,n_Bsites,n_Csites,depth)
    randMeasB(reg,n_Asites,n_Bsites,n_Csites,n_meas)
    return cmi(reg,n_Asites,n_Bsites,n_Csites)
end

function simHPClifMeas(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas::Int64)
    reg = genInitHP(n_Asites,n_Bsites,n_Csites)
    applyHPClif(reg,n_Asites,n_Bsites,n_Csites)
    randMeasB(reg,n_Asites,n_Bsites,n_Csites,n_meas)
    return cmi(reg,n_Asites,n_Bsites,n_Csites)
end

function simTMClifBMeas(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas::Int64)
    reg = genInitABC(n_Asites,n_Bsites,n_Csites)
    applyHPClif(reg,n_Asites,n_Bsites,n_Csites)
    randBellMeasB(reg,n_Asites,n_Bsites,n_Csites,n_meas)
    return cmi(reg,n_Asites,n_Bsites,n_Csites)
end

function simHPBWMeas(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas::Int64,depth::Int64)
    reg = genInitHP(n_Asites,n_Bsites,n_Csites)
    applyHPBW(reg,n_Asites,n_Bsites,n_Csites,depth)
    randMeasB(reg,n_Asites,n_Bsites,n_Csites,n_meas)
    return cmi(reg,n_Asites,n_Bsites,n_Csites)
end

function simHPClifBMeas(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas::Int64)
    reg = genInitHP(n_Asites,n_Bsites,n_Csites)
    applyHPClif(reg,n_Asites,n_Bsites,n_Csites)
    randBellMeasB(reg,n_Asites,n_Bsites,n_Csites,n_meas)
    return cmi(reg,n_Asites,n_Bsites,n_Csites)
end

function simHPBWBMeas(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas::Int64,depth::Int64)
    reg = genInitHP(n_Asites,n_Bsites,n_Csites)
    applyHPBW(reg,n_Asites,n_Bsites,n_Csites,depth)
    randBellMeasB(reg,n_Asites,n_Bsites,n_Csites,n_meas)
    return cmi(reg,n_Asites,n_Bsites,n_Csites)
end

function simTwoClifMeas(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas::Int64)
    reg = genInitHP(n_Asites,n_Bsites,n_Csites)
    applyTwoRandClif(reg,n_Asites,n_Bsites,n_Csites)
    randMeasB(reg,n_Asites,n_Bsites,n_Csites,n_meas)
    return cmi(reg,n_Asites,n_Bsites,n_Csites)
end

function simTwoClifBMeas(n_Asites::Int64,n_Bsites::Int64,n_Csites::Int64,n_meas::Int64)
    reg = genInitHP(n_Asites,n_Bsites,n_Csites)
    applyTwoRandClif(reg,n_Asites,n_Bsites,n_Csites)
    randBellMeasB(reg,n_Asites,n_Bsites,n_Csites,n_meas)
    return cmi(reg,n_Asites,n_Bsites,n_Csites)
end