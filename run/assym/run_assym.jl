 # Timesteps, gridpoints, length
const M = 1_000_000
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


u = 10.
D = 1e-5
φ1 = -0
φ2 = 0
α = 15.
param = (u, α, D, φ1, φ2)
o = ones(N÷4)/√2
l = LinRange(-1/√2, 1/√2, N÷4)
φ0 = [cat(o, -l, -o, l, dims=1) cat(-l, -o, l, o, dims=1) ]
φ0 = [cat(o, o, -o, -o, dims=1) cat(o, -o, -o, o, dims=1) ]

@time run_euler(param; φ0=φ0, step="C4");

