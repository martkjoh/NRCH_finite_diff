using Base.Threads

# Timesteps, gridpoints, length
const M = 400_000_000
const N = 500
const L = N/5.
const c = .01 # dt/dx^4 = c
const frames = 4_000
# Folder for data
write_folder = "article/data/4/"
include("../../numerics.jl")

D = 2e-4
u = 10.
α = 15
φ1, φ2  = 0, 0

φ1s = [-0.2, -0.8, -0.4]
φ2s = [-0.2, 0., -0.2]
αs = [10., 15., 5.]
nums = ["7", "8", "9"]
n, = axes(nums)

@time @threads for i in n
    α, bφ1, bφ2, num = αs[i], φ1s[i], φ2s[i], nums[i]
    param = (u, α, D, bφ1, bφ2)
    name = "(VID"*num*")"
    @time run_euler(param, step="C4", name_app=name)
end
