 # Timesteps, gridpoints, length
const M = 5_000_000
const N = 200
const L = N/10.
const c = .01 # dt/dx^4 = c
const frames = 10_000
# P rint progress bar
pr = true
# Folder for data
write_folder = "data/assym/sep/"
rm(write_folder, recursive=true, force=true)
mkdir(write_folder[1:end-1])

include("../../numerics.jl")


using Base.Threads

us = [1., 10., 100., 1000.]
D = 1e-5 
φ1 = - 1/sqrt(2)
φ2 = 0
α = 0.

@time @threads for u in us
    param = (u, α, D, φ1, φ2)
    @time run_euler(param, step="C4", init=4)
end
