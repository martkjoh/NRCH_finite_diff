 # Timesteps, gridpoints, length
const M = 100_000
const N = 50
const L = N/5.
const c = .01 # dt/dx^4 = c
# P rint progress bar
pr = true
# Folder for data
write_folder = "data/test/"
rm(write_folder, recursive=true, force=true)
mkdir(write_folder[1:end-1])

include("numerics.jl")


u = 10.
D = 1e-5
bφ = -.7
α = 4.5
param = (u, α, D, bφ)
@time run_euler(param);
