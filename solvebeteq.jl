using Plots


function  groundofaiferro(N::Int64,sz=0,stp = 100)
    r = Int(round(N/2 - sz))
    z = fill(0.0,r)
    IA = [ N/4-i+1/2-sz/2 for i in 1:r]
    for i in 1:stp
        for j in 1:r
            z[j] = tan(π/N*IA[j]+1/(2*N)*sum(x->2*atan((z[j]-x)/2),z))
        end
    end
    plot(IA/N, z)
    k = π .- 2*atan.(z)
    plot(IA/N, k)
    energy = sum(x->-2. / (1. + x^2),z)/N
end

N = 32
E = Float64[]
for i in 1:N/2
    push!(E,groundofaiferro(N,N/2-i))
end

scatter([1:N/2]/N,E)
