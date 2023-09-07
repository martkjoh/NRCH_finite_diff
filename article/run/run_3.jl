using Base.Threads

# Timesteps, gridpoints, length
const M = 100_000_000
const N = 500
const L = N/5.
const c = .01 # dt/dx^4 = c
const frames = 10_000
# Folder for data
write_folder = "article/data/3/"
include("../../numerics.jl")

u = 40.
D = 1e-7
α = 0
φs = [-.95, -.9, -.8]
n, = axes(φs)

@time @threads for i in n
    bφ = φs[i]
    param = (u, α, D, bφ, 0., 0)
    @time run_euler(param; init=2, name_app="(VID0)");
end
