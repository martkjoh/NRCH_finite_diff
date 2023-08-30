# Timesteps, gridpoints, length
const M = 10_000_000
const N = 500
const L = N/5.
const c = .01 # dt/dx^4 = c
# P rint progress bar
pr = true
# Folder for data
write_folder = "data/traveling2/"
rm(write_folder, recursive=true, force=true)
mkdir(write_folder[1:end-1])

include("../../numerics.jl")


using Base.Threads
u = 10.
D = 1e-5
bφ = 0
α = 1.
param = (u, α, D, bφ)

inits = [false, false, false, true]
name_apps = ["1", "2", "3", "4"]

@time @threads for (init, name_app) in collect(zip(inits, name_apps))
    @time run_euler(param; init=init, name_app=name_app)
end

