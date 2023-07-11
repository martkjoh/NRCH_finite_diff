using DelimitedFiles
using BenchmarkTools
using Random
import KernelAbstractions.Extras.LoopInfo: @unroll

const dx = L / N
const dt = round(c * (dx)^4; sigdigits=6)
const frames = 1_000
const skip = div(M, frames)
const rng = MersenneTwister(1234)

print("T    = ", round(M*dt; sigdigits=6), '\n')
print("dT   = ", dt, '\n')
param_names = ["u, -r", "a", "D", "phi", "N", "L", "dt"]


# Finite difference coeficcients, https://en.wikipedia.org/wiki/Finite_difference_coefficient
# const co = [1/2] / dx
# const co = [2/3, -1/12] / dx
# const co = [3/4, -3/20, 1/60] / dx
const co = [ 4/5, -1/5, 4/105, -1/280 ] / dx

const order = length(co)

@inline ind(i) = mod(i-1, N)+1

@inline function ∇(A, i)
	d = 0.
	@inbounds @unroll for j in 1:order
		d += co[j]*( A[ind(i+j)] - A[ind(i-j)] ) 
	end
	return d
end

@inline function ∇²(A, i)
	d = 0.
	@inbounds @unroll for j in 1:order
		d += co[j]*( ∇(A, i+j) - ∇(A, i-j)  ) 
	end
	return d
end

function euler!(φ, μ, δφ, ξ, param_r)
    u, α, σ = param_r
    @inbounds for i in 1:N
        @views ruφ² = u * (-1 + (φ[i, 1]^2 + φ[i, 2]^2 ))
        @views μ[i, 1] = ruφ² * φ[i, 1] - ∇²(φ[:, 1], i) + α * φ[i, 2]
        @views μ[i, 2] = ruφ² * φ[i, 2] - ∇²(φ[:, 2], i) - α * φ[i, 1]
    end 
    randn!(rng, ξ)
    ξ .*= σ
    @inbounds for i in 1:N
        @views δφ[i, 1] = ( ∇²(μ[:, 1], i) - ∇(ξ[:, 1], i) ) * dt
        @views δφ[i, 2] = ( ∇²(μ[:, 2], i) - ∇(ξ[:, 2], i) ) * dt
    end
end

function loop!(φt,  φ, μ, δφ, ξ, param_r)
    for i in axes(φt, 1)[2:end]
        for _ in 1:skip
            euler!(φ, μ, δφ, ξ, param_r)
            φ .+= δφ
        end
        check(φ, i)
        φt[i,:,:] .= φ
    end
    print('\n')
end

function run_euler(param)
    u, α, D, bφ = param
    σ = sqrt(2 * D / dt / dx)

    x = LinRange(0, L-dx, N)
    φ = .1 * [sin.(2*pi*x/L) cos.(2*pi*x/L)]

    param_r = (u, α, σ)

    φt = zeros(frames, N, 2)
    μ = zeros(N, 2)
    δφ = zeros(N, 2)
    ξ = zeros(N, 2)
    
    φ[:,1] .+= bφ
    φt[1,:,:] .= φ

    loop!(φt, φ, μ, δφ, ξ, param_r)

    param_write = (u, α, D, bφ, N, L, dt)
    
    write_file(φt, param_write)
end

##############
# Utillities #
##############

function write_file(φt, param)
    filename = join(
        param_names[i] * '=' * string(param[i]) * "_"
        for i in range(1, length(param_names))
        )[1:end-1]
    w = reshape(permutedims(φt, (3, 2, 1)), (frames*2*N))
    writedlm(write_folder*filename*".txt", w)
end

function check(φ, i)
    n = frames//10
    if any(isnan, φ)
        throw(ErrorException("Error, NaN detected" ))
    end

    if (div(i,n)) - div(i-1,n) == 1 && pr
        print("\r"*"|"^div(i,n))
    end
end

