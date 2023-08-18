using DelimitedFiles
using BenchmarkTools
using Random

const dx = L / N
const dt = round(c * (dx)^4; sigdigits=6)
const frames = 10_000
const skip = div(M, frames)

print("T    = ", round(M*dt; sigdigits=6), '\n')
print("dT   = ", dt, '\n')
param_names = ["u, -r", "a", "D", "phi", "N", "L", "T", "dt"]

@inline ind(i) = mod(i-1, N)+1

# Finite difference coeficcients: https://en.wikipedia.org/wiki/Finite_difference_coefficient

@inline ∇(A, i)  = ( 8*(A[ind(i + 1)] - A[ind(i - 1)]) - (A[ind(i + 2)] - A[ind(i - 2)]) ) / (12*dx)
@inline ∇²(A, i)  = ( 8*( ∇(A,i + 1) - ∇(A,i - 1)) - (∇(A,i + 2) - ∇(A,i - 1)) ) / (12*dx)
@inline t∇²(A, i) = (-5/2*A[i] + 4/3*(A[ind(i+1)] + A[ind(i-1)]) - 1/12*(A[ind(i+2)] + A[ind(i-2)])) / dx^2


function euler!(φ, μ, δφ, ξ, param_r)
    u, α, σ = param_r
    @inbounds for i in 1:N
        @views ruφ² = u * (-1 + (φ[i, 1]^2 + φ[i, 2]^2 ))
        @views μ[i, 1] = ruφ² * φ[i, 1] - t∇²(φ[:, 1], i) + α * φ[i, 2]
        @views μ[i, 2] = ruφ² * φ[i, 2] - t∇²(φ[:, 2], i) - α * φ[i, 1]
    end 
    randn!(ξ)
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

function run_euler(param; init=0, n=1, name_app="")
    u, α, D, bφ = param
    σ = sqrt(2 * D / dt / dx)

    x = LinRange(0, L-dx, N)
    if init==0 φ = zeros(N, 2)
    elseif init==1 φ = [ sin.(2π.*x/L * n)   cos.(2π.*x/L * n) ]
    elseif init==2 φ = [ .05 .* cos.(2 .* 2π.*x/L)   .5 .* cos.(1 .* 2π.*x/L) ]
    elseif init==3 φ = [ 0 .* x     0.2.* cos.(1 .* 2π.*x/L) ] 
    end

    param_r = (u, α, σ)

    φt = zeros(frames, N, 2)
    μ = zeros(N, 2)
    δφ = zeros(N, 2)
    ξ = zeros(N, 2)
    
    φ[:,1] .+= bφ
    φt[1,:,:] .= φ

    loop!(φt, φ, μ, δφ, ξ, param_r)

    param_write = (u, α, D, bφ, N, L, M*dt, dt)
    
    write_file(φt, param_write; name_app=name_app)
end

##############
# Utillities #
##############

function write_file(φt, param; name_app=name_app)
    filename = join(
        param_names[i] * '=' * string(param[i]) * "_"
        for i in range(1, length(param_names))
        )[1:end-1]
    if (name_app!="") filename *= "*" *  name_app end
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

