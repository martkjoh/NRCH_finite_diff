 # Timesteps, gridpoints, length
const M = 2_000_000
const N = 500
const L = N/5.
const c = .01 # dt/dx^4 = c
# P rint progress bar
pr = true
# Folder for data
write_folder = "data/assym/test/"
rm(write_folder, recursive=true, force=true)
mkdir(write_folder[1:end-1])

include("../../numerics.jl")

D = 1e-5
u = 10.
α = 15
φ1, φ2  = 0, 0

# circ = -(1/√6 - 1e-3)
# φ1, φ2  = circ, circ
# φ1, φ2  = -.6, -.3
param = (u, α, D, φ1, φ2)

o = ones(N÷4)/√2
l = LinRange(-1/√2, 1/√2, N÷4)
# φ0 = [cat(o, -l, -o, l, dims=1) cat(-l, -o, l, o, dims=1) ]
φ0 = [cat(o, o, -o, -o, dims=1) cat(o, -o, -o, o, dims=1) ]

# o = ones(N÷2)/√2
# φ0 = [cat(o, -o, dims=1) cat(o, -o, dims=1) ]

# o1 = ones(50)/√2
# o2 = ones(100)/√2
# φ0 = [cat(o1, o1, -o2, -o2, dims=1) cat(o1, -o1, -o2, o2, dims=1) ]
# bφ = sum(φ0, dims=1)/N
# φ0[:, 1] .-= bφ[1]
# φ0[:, 2] .-= bφ[2]

# @time run_euler(param; step="C4", φ0=φ0);


# φ1, φ2  = -.41, -.41
# # param = (u, α, D, φ1, φ2)

@time run_euler(param; step="C4", init=1);