 # Timesteps, gridpoints, length
const M = 200_000_000
const N = 500
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
bφ = -.8
α = 4.
param = (u, α, D, bφ)
@time run_euler(param);
