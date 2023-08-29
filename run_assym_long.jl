# Timesteps, gridpoints, length
const M = 10_000_000
const N = 500
const L = N/5.
const c = .01 # dt/dx^4 = c
# P rint progress bar
pr = true
# Folder for data
write_folder = "assym/data/long/"
rm(write_folder, recursive=true, force=true)
mkdir(write_folder[1:end-1])

include("numerics.jl")


using Base.Threads
u = 10.
D = 1e-5

n = 4
φ1s = cat( - LinRangframese(0, 1, n), -  LinRange(0, 1, n), dims=1 )
φ2s = cat( zeros(n), -  LinRange(0, 1, n), dims=1 )
αs = [0. 2. 4.]
αφ = [(α, φ1, φ2) for α in αs for φ1 in φ1s for φ2 in φ2s]

@time @threads for (α, φ1, φ2) in αφ
    param = (u, α, D, φ1, φ2)
    @time run_euler(param, step="C4")
end

