using Yao
using BitBasis
using LinearAlgebra
using Statistics
using QuAlgorithmZoo: Sequence,Adam
using Plots
using Evolutionary

include("circuit.jl")
include("ground.jl")

function mydispatch!(circuit,x)
        
end


function cmaes_train(mytcircuit,myfcircuit,model;VG = nothing, maxiter = 200,α = 0.3,nbatch=1024)
        function loss(x)

                energy(tcircuit,model,nbatch=nbatch/nspin(model))
        end

        ci = CMAESIter(,num_offsprings = 50,num_parents = 10)
        loss_history = Float64[]
        for (i,info) in enumerate(ci)
                push!(loss_history,best(info)[2])
                if i%10 == 0
                        print("Iter $i, E/N = $(loss_history[end])")
                        fid = fidelity(fcircuit,VG)
                        println(VG isa Nothing ? "" : ",fidelity = $(fid)")
                end
                if i>=maxiter
                        return loss_history
                else i+=1
                end
        end
end


nphy = 6;
nex = 3;
model = Heisenberg(nphy;periodic = false)
h = hamiltionian(model)
res = eigen(mat(h)|>Matrix)
EG = res.values[1]/nspin(model)
@show EG
VG = res.vectors[:,1]


ρ = Float64[1.,1.,1.]
θ = Float64[0.0,0.0,0.0]
mfcircuit = fcircuit(nphy,nex,θ,ρ)
mycircuit = tcircuit(nphy,nex,θ,ρ)
