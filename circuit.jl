using Yao
using LinearAlgebra
using Statistics:mean
using Random

rng = Random.GLOBAL_RNG
include("non-unitary.jl")

U(θ) = matblock(ComplexF64[1 0 0 0; 0 1/sqrt(2)*exp(im*θ) 1/sqrt(2) 0;0 -exp(im*θ)/sqrt(2) 1/sqrt(2) 0;0 0 0 1])
S(ρ) = [1 0 0 0; 0 ρ 0 0; 0 0 1 0; 0 0 0 1]


function initialblock(nphy,nex)
        ntot = nphy + nex
        nan = 1
        chain(ntot + nan, put(ntot + nan, i=>X) for i in nphy+1:nphy+nex);
end

function R(nphy::Int,nex::Int,seqlayer::Int,seqR::Int,θ,ρ)
        ntot = nphy + nex
        nan = 1
        layerbot = nphy + seqlayer
        bot = layerbot - seqR + 1
        mid = bot - 1
        top = ntot + nan
         chain(
         ntot + nan,
         put(ntot + nan,(mid, bot)=>U(0)),
         lcu(ntot + nan,top,mid,bot,S(ρ)),
         put(ntot + nan,(mid,bot)=>U(θ))
)
end

function layer(nphy::Int,nex::Int,seqlayer::Int,θ,ρ)
        ntot = nphy + nex
        nan = 1
        layerbot = nphy + seqlayer
        layertop = seqlayer
        chain(ntot+nan,
        chain(ntot+nan,R(nphy,nex,seqlayer,seqR,θ,ρ)  for seqR in 1:nphy),
        put(ntot + nan, (layertop=>ConstGate.P1)))
end

function fcircuit(nphy::Int,nex::Int,θ,ρ)
        ntot = nphy + nex
        nan = 1
        chain(ntot + nan, initialblock(nphy,nex),
        chain(ntot + nan,
        layer(nphy,nex,seqlayer,θ[seqlayer],ρ[seqlayer]) for seqlayer in 1:nex))
end

function allmeasure(nphy,nex)
        ntot = nphy + nex
        nan = 1
        chain(ntot + nan,
        Measure{ntot + nan,1,AbstractBlock,typeof(rng)}(rng,Z,(i,),0,false) for i in (nex+1) : (nex+nphy))
end

tcircuit(nphy,nex,θ::Array{Float64,1},ρ::Array{Float64,1}) = chain(nphy + nex + 1,
        fcircuit(nphy,nex,θ,ρ),
        allmeasure(nphy,nex)
        )
