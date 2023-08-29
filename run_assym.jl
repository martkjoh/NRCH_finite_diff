 # Timesteps, gridpoints, length
const M = 10_000_000
const N = 500
const L = N/5.
const c = .01 # dt/dx^4 = c
# P rint progress bar
pr = true
# Folder for data
write_folder = "assym/data/test/"
rm(write_folder, recursive=true, force=true)
mkdir(write_folder[1:end-1])

include("numerics.jl")


u = 10.
D = 1e-5
φ1 = 0.
φ2 = 0
α = 10.
param = (u, α, D, φ1, φ2)
@time run_euler(param; init=5, name_app="", step="C4");
