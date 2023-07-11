# Timesteps, gridpoints, length
const M = 2_000_000
const N = 200
const L = 10.
const c = .1
# P rint progress bar
pr = true
# Folder for data
write_folder = "data/sym/"
rm(write_folder, recursive=true, force=true)
mkdir(write_folder[1:end-1])

include("numerics.jl")


u = 10.
D = 4e-6
bφ = -.707
α = 4.5
param = (u, α, D, bφ)
@time run_euler(param);



# using Base.Threads
# u = 10.
# D = 0.0001

# αs = LinRange(0, 6, 2)
# φs = [-.8, -0.71, -.6, -.5]
# αφ = [(α,φ) for α in αs for φ in φs]

# @time @threads for (α, bφ) in αφ
#     param = (u, α, D, bφ)
#     @time run_euler(param)
# end

