using Base.Threads

# Timesteps, gridpoints, length
const M = 400_000_000
const N = 500
const L = N/5.
const c = .01 # dt/dx^4 = c
# const frames = 10_000
const frames = 4_000
# Folder for data
write_folder = "article/data/1/"
include("../../numerics.jl")

u = 10.
D = 2e-4
bφ2 = 0.

φs = [-0.5, -0.7, -0.8, -0.5]
αs = [0., 2., 4., 6.,]
nums = ["1", "3", "4", "5"]
n, = axes(nums)

@time @threads for i in n
    α, bφ1, num = αs[i], φs[i], nums[i]
    param = (u, α, D, bφ1, bφ2)
    name = "(VID"*num*")"
    @time run_euler(param, name_app=name)
end
