using Yao
using BitBasis
using LinearAlgebra
using Statistics
using QuAlgorithmZoo: Sequence,Adam
using Plots
using Evolutionary

include("circuit.jl")
include("ground.jl")

nphy = 6;
nex = 3;
model = Heisenberg(nphy;periodic = false)
h = hamiltionian(model)
res = eigen(mat(h)|>Matrix)
EG = res.values[1]/nspin(model)
@show EG
VG = res.vectors[:,1]
mfcircuit = fcircuit(nphy,nex,θ,ρ,ϵ)
mmeasure = allmeasure(nphy,nex)
mtcircuit = chain(nphy+nex+nphy*nex,mfcircuit,mmeasure)


ρ = Float64[1.,1.,1.]
θ = Float64[0.0,0.0,0.0]
ϵ = 0.5

function cmaes_train(circuit,model,nphy,nex;maxiter = 200,α = 0.3,nbatch=1024)
        loss = energy(cir)
        ci = CMAESIter()
end
