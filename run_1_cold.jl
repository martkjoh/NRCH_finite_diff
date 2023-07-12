# Timesteps, gridpoints, length
const M = 200_000_000
const N = 500
const L = N/5.
const c = .01
# P rint progress bar
pr = true
# Folder for data
write_folder = "data/sym_long_cold/"
rm(write_folder, recursive=true, force=true)
mkdir(write_folder[1:end-1])

include("numerics.jl")


using Base.Threads
u = 10.
D = 1e-12

αs = LinRange(0, 6, 4)
φs = [-.8, -0.7, -.6, -.2]
αφ = [(α,φ) for α in αs for φ in φs]

@time @threads for (α, bφ) in αφ
    param = (u, α, D, bφ)
    @time run_euler(param)
end

