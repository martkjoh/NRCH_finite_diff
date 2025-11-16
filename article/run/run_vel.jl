using Base.Threads

# Timesteps, gridpoints, length
const M = 1_000_000
const N = 1000
const L = N/5.
const c = .01 # dt/dx^4 = c
const frames = 100
# Folder for data
include("../../numerics.jl")

u = 10
D = 0e-5

SS = [50, 75, 100, 125, 150, 200, 250, 300,]

@threads for SIZE in SS
    write_folder = "article/data/vel/$(SIZE)/"
    rm(write_folder, recursive=true, force=true)

    for α in LinRange(0., 2., 20)
        φ0 = [ cat(-ones(SIZE), ones(SIZE), -ones(N-SIZE-SIZE), dims=1) cat(-ones(SIZE÷2), ones(SIZE), -ones(N-(SIZE÷2)-SIZE), dims=1) ] ./ √2
        φ1, φ2 = mean(φ0[:, 1]), mean(φ0[:, 2])
        φ0[:, 1] .= (φ0[:, 1] .- φ1)
        φ0[:, 2] .= (φ0[:, 2] .- φ2)
        param = (u, α, D, φ1, φ2)
    
        @time run_euler(param, step="C4", φ0=φ0, write_folder=write_folder)
    end
end
