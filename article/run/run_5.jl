using Base.Threads

# Timesteps, gridpoints, length
const M = 5_000_000
const N = 200
const L = N/5.
const c = .01 # dt/dx^4 = c
const frames = 1000
# Folder for data
write_folder = "article/data/5/"

include("../../numerics.jl")


u = 10
D = 1e-5 
φ1, φ2 = -.2, -.1
α = 2.
param = (u, α, D, φ1, φ2)
o = ones(N÷4)/√2
φ0 = [cat(o, o, -o, -o, dims=1) cat(o, -o, -o, o, dims=1) ]

@time run_euler(param, step="C4", φ0=φ0, name_app="(VID6)")
