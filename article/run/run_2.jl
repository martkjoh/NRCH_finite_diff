using Base.Threads

# Timesteps, gridpoints, length
const M = 400_000_000
const N = 500
const L = N/5.
const c = .01 # dt/dx^4 = c
const frames = 4_000
# Folder for data
write_folder = "article/data/2/"
include("../../numerics.jl")

u = 10.
D = 2e-4
bφ = -.2
α = 1.
param = (u, α, D, bφ, 0., 0)

inits = [0, 0, 0, 1]
name_apps = ["1", "2", "3", "4"]
n, = axes(inits)

@time @threads for i in n
    init, name_app = inits[i], name_apps[i]
    @time run_euler(param; init=init, name_app=name_app)
end

