 # Timesteps, gridpoints, length
const M = 10_000_000
const N = 500
const L = N/5.
const c = .01 # dt/dx^4 = c
const frames = 1000
# P rint progress bar
pr = true
# Folder for data
write_folder = "data/test/"
rm(write_folder, recursive=true, force=true)
mkdir(write_folder[1:end-1])

include("../../numerics.jl")


u = 10.
D = 1e-5
bφ1 = 0.
bφ2 = 0.
α = 20.
param = (u, α, D, bφ1, bφ2)
@time run_euler(param; step="SO2", init=1);
