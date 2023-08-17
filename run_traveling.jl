# Timesteps, gridpoints, length
const M = 25_000_000
const N = 500
const L = N/5.
const c = .01 # dt/dx^4 = c
# P rint progress bar
pr = true
# Folder for data
write_folder = "data/traveling/"
rm(write_folder, recursive=true, force=true)
mkdir(write_folder[1:end-1])

include("numerics.jl")


using Base.Threads
u = 10.
D = 1e-5
bφ = -.2
α = .2
param = (u, α, D, bφ)

circs = [false, false, false, true]
name_apps = ["1", "2", "3", "4"]

@time @threads for (circ, name_app) in collect(zip(circs, name_apps))
    @time run_euler(param; circ=circ, name_app=name_app)
end

