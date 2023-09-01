 # Timesteps, gridpoints, length
const M = 10_000_000
const N = 200
const L = N/10.
const c = .01 # dt/dx^4 = c
# P rint progress bar
pr = true
# Folder for data
write_folder = "data/assym/sep3/"
rm(write_folder, recursive=true, force=true)
mkdir(write_folder[1:end-1])

include("../../numerics.jl")


using Base.Threads

u = 10
D = 1e-5 
φ1, φ2 = 0.45, 0.45
α = 0. 
param = (u, α, D, φ1, φ2)
o = ones(N÷4)/√2
φ0 = [cat(o, o, -o, -o, dims=1) cat(o, -o, -o, o, dims=1) ]

@time run_euler(param, step="C4", φ0=φ0)

 