using Yao
using BitBasis
using LinearAlgebra
using Statistics
using QuAlgorithmZoo: Sequence

abstract type AbstractModel{D} end

nspin(model::AbstractModel) = prod(model.size)


struct Heisenberg{D} <: AbstractModel{D}
    size::NTuple{D,Int}
    periodic::Bool
    Heisenberg(size::Int...;periodic::Bool) = new{length(size)}(size,periodic)
end

heisenberg_ij(nbit::Int, m::Int, n::Int = m + 1) = put(nbit,m=>X) * put(nbit,n=>X) + put(nbit,m=>Y) * put(nbit,n=>Y) + put(nbit,m=>Z) * put(nbit,n=>Z)
#heisenberg_iijj(nbit::Int, m::Int, n::Int = m + 2) = put(nbit,m=>X) * put(nbit,n=>X) + put(nbit,m=>Y) * put(nbit,n=>Y) + put(nbit,m=>Z) * put(nbit,n=>Z)

function get_bonds(model::AbstractModel{1})
    nbit, = model.size
    [(i,i%nbit+1) for i in 1:(model.periodic ? nbit : nbit-1)]
end

function get_nextnear_bonds(model::AbstractModel{1})
    nbit, = model.size
    [(i,i%nbit+2) for i in 1:(model.periodic ? nbit : nbit-2)]
end

function get_bonds(model::AbstractModel{2})
    m,n = model.size
    cis = LinearIndices(model.size)
    bonds = Tuple{Int,Int,Float64}[]
    for i = 1:m, j = 1:n
        (i!=m || model.periodic) && push!(bonds,(cis[i,j],cis[i%m+1,j]))
        (j!=n || model.periodic) && push!(bonds,(cis[i,j],cis[i,j%n+1]))
    end
    bonds
end

function hamiltionian(model::Heisenberg)
    nbit = nspin(model)
    nearest = sum(x->heisenberg_ij(nbit,x[1],x[2]),get_bonds(model))*0.25;
    #nextnear = sum(x->heisenberg_ij(nbit,x[1],x[2]),get_nextnear_bonds(model))*0.1;
    return nearest #+ nextnear
end


function gensample(circuit,operator;nbatch=1024)
    mblocks = collect_blocks(Measure,circuit)
    for m in mblocks
        m.operator = operator
    end
    reg = zero_state(nqubits(circuit);nbatch=nbatch)
    reg |> circuit
    mblocks
end

function ising_energy(circuit,bonds,basis;nbatch = nbatch)
    mblocks = gensample(circuit,basis;nbatch = nbatch)
    nspin = length(mblocks)
    local eng = 0.0
    for (a,b) in bonds
        eng += mean(mblocks[a].results.*mblocks[b].results)
    end
    eng/=4
end

function energy(circuit,model::Heisenberg;nbatch = 1024)
    bonds = get_bonds(model)
    sum(basis->ising_energy(circuit,bonds,basis;nbatch = nbatch),[X,Y,Z])
end


#hei_model = Heisenberg(4;periodic = false)

function cos_fit_min(x::Array{Float64,1},y::Array{Float64,1})
#   A = std(y);
#    b = mean(y);
#    w = 1.0;
#    ph = 0.0;
#    x0 = vec(x);
#    p = [A,w,ph,b];
#    fun(x0,p) = p[1].*sin.(p[2].*x0 .+ p[3]) .+ p[4];
#    fit = curve_fit(fun,x0,y0,p);
#    return fit.param[1] > 0 ? mod2pi(3/2*π - fit.param[3]) : mod2pi(π/2 - fit.param[3])
    b = (y[1] + y[3])/2
    ϕ = atan((y[1] - b)/(y[2] - b)) - x[1]
    A = (y[1] - b)/sin(x[1] + ϕ)
    if A < 0
        A = -A
        ϕ = mod2pi(ϕ+π)
    else
        ϕ = mod2pi(ϕ)
    end
    return A,ϕ,b
end

function optimalr(A,ϕ,b)
    A = median(A)
    ϕ = median(ϕ)
    b = median(b)
    #A = mean(A)
    #ϕ = mean(ϕ)
    #b = mean(b)
    return mod2pi(3/2*π - ϕ)
end

function fidelity(circuit,nphy,nex,VG)
    psi = zero_state(nphy+nex+1,)
end
