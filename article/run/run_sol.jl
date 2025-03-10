# Timesteps, gridpoints, length
const M = 10_000_000
const N = 200
const L = N/5.
const c = .01 # dt/dx^4 = c
const frames = 100
# P rint progress bar
pr = true
# Folder for data
write_folder = "article/data/sol/"
rm(write_folder, recursive=true, force=true)
mkdir(write_folder[1:end-1])

include("../../numerics.jl")


u = 40.
D = 1e-7
α = 0

bφ = -.8
param = (u, α, D, bφ, 0)
@time run_euler(param; init=2, name_app="1");

