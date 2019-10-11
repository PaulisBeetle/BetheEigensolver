using Yao
using LinearAlgebra
using Statistics:mean
using Random

rng = Random.GLOBAL_RNG

U(θ) = matblock(ComplexF64[1 0 0 0; 0 1/sqrt(2)*exp(im*θ) 1/sqrt(2) 0;0 -exp(im*θ)/sqrt(2) 1/sqrt(2) 0;0 0 0 1])
S(ρ) = [1 0 0 0; 0 ρ 0 0; 0 0 1 0; 0 0 0 1]
T(ρ,ϵ) = matblock(ComplexF64[sin(ϵ*S(ρ)) cos(ϵ*S(ρ)); -cos(ϵ*S(ρ)) sin(ϵ*S(ρ))])

function initialblock(nphy,nex)
        ntot = nphy + nex
        nan = nex * nphy
        chain(ntot + nan, put(ntot + nan, i=>X) for i in nphy+1:nphy+nex);
end

function R(nphy,nex,seqlayer,seqR,θ,ρ,ϵ)
        ntot = nphy + nex
        nan = nex * nphy
        layerbot = ntot - nex + seqlayer
        bot = layerbot - seqR + 1
        mid = bot - 1
        top = ntot + (seqlayer - 1) * nphy + seqR
         chain(
         ntot+nan,
         put(ntot + nan,(mid, bot)=>U(0)),
         put(ntot + nan ,(top,mid,bot)=>T(ρ,ϵ)),
         put(ntot + nan,(mid,bot)=>U(θ)),
         Measure{ntot+nan,1,AbstractBlock,typeof(rng)}(rng,Z,(top,),0,false)
)
end

function layer(nphy,nex,seqlayer,θ,ρ,ϵ)
        ntot = nphy + nex
        nan = nphy * nex
        layerbot = ntot - nex + seqlayer
        layertop = layerbot - nphy
        nmlayer = chain(ntot + nan,
        R(nphy,nex,seqlayer,seqR,θ,ρ,ϵ)  for seqR in 1:nphy)
        chain(ntot + nan, nmlayer, Measure{ntot+nan,1,AbstractBlock,typeof(rng)}(rng,Z,(layertop,),0,false))
end

function fcircuit(nphy,nex,θ,ρ,ϵ)
        ntot = nphy + nex
        nan = nphy * nex
        chain(ntot + nan, initialblock(nphy,nex),
        chain(ntot + nan,
        layer(nphy,nex,seqlayer,θ[seqlayer],ρ[seqlayer],ϵ) for seqlayer in 1:nex))
end

function allmeasure(nphy,nex)
        ntot = nphy + nex
        nan = nphy * nex
        chain(ntot + nan,
        Measure{ntot + nan,1,AbstractBlock,typeof(rng)}(rng,Z,(i,),0,false) for i in (nex+1) : (nex+nphy))
end

tcircuit(nphy,nex,θ::Array{Float64,1},ρ::Array{Float64,1},ϵ) = chain(nphy + nex + nphy * nex,
        fcircuit(nphy,nex,θ,ρ,ϵ),
        allmeasure(nphy,nex)
        )
