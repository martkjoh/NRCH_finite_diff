# Timesteps, gridpoints, length
const M = 5_000_000
const N = 500
const L = N/5.
const c = .01 # dt/dx^4 = c
# P rint progress bar
pr = true
# Folder for data
write_folder = "assym/data/long/"
rm(write_folder, recursive=true, force=true)
mkdir(write_folder[1:end-1])

include("../../numerics.jl")
using Base.Threads

u = 10.
D = 1e-5

n = 6
m = 4
l = -LinRange(0, 1, n)
φ1s = cat( l, l, l, dims=1 )
φ2s = cat( 0 .* l, 1/2 .* l, l, dims=1 )
αs = LinRange(0, 15, m)
αφ = [(α, φi[1], φi[2]) for α in αs for φi in zip(φ1s, φ2s)]

@time @threads for (α, φ1, φ2) in αφ
    param = (u, α, D, φ1, φ2)
    @time run_euler(param, step="C4")
end

