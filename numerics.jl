using DelimitedFiles
using BenchmarkTools
using Random
using ProgressMeter

seed = 1
rng = Random.Xoshiro(seed)

const dx = L / N
const dt = round(c * (dx)^4; sigdigits=6)
const skip = div(M, frames)
# P rint progress bar
pr = true

print("T    = ", round(M*dt; sigdigits=6), '\n')
print("dT   = ", dt, '\n')
param_names = ["u, -r", "a", "D", "phi1", "phi2", "N", "L", "T", "dt"]

@inline ind(i) = mod(i-1, N)+1

# Finite difference coefficients: https://en.wikipedia.org/wiki/Finite_difference_coefficient

@inline ∇(A, i)  = ( 8*(A[ind(i + 1)] - A[ind(i - 1)]) - (A[ind(i + 2)] - A[ind(i - 2)]) ) / (12*dx)
@inline ∇²(A, i)  = ( 8*( ∇(A,i + 1) - ∇(A,i - 1)) - (∇(A,i + 2) - ∇(A,i - 1)) ) / (12*dx)
@inline t∇²(A, i) = (-5/2*A[i] + 4/3*(A[ind(i+1)] + A[ind(i-1)]) - 1/12*(A[ind(i+2)] + A[ind(i-2)])) / dx^2


function euler_SO2!(φ, μ, δφ, ξ, param_r)
    u, α, σ = param_r
    @inbounds for i in 1:N
        @views ruφ² = u * (-1 + (φ[i, 1]^2 + φ[i, 2]^2 ))
        @views μ[i, 1] = ruφ² * φ[i, 1] - t∇²(φ[:, 1], i) + α * φ[i, 2]
        @views μ[i, 2] = ruφ² * φ[i, 2] - t∇²(φ[:, 2], i) - α * φ[i, 1]
    end 
    randn!(rng, ξ)
    ξ .*= σ
    @inbounds for i in 1:N
        @views δφ[i, 1] = ( ∇²(μ[:, 1], i) - ∇(ξ[:, 1], i) ) * dt
        @views δφ[i, 2] = ( ∇²(μ[:, 2], i) - ∇(ξ[:, 2], i) ) * dt
    end
end


function euler_C4!(φ, μ, δφ, ξ, param_r)
    u, α, σ = param_r
    @inbounds for i in 1:N
        @views ruφ₁² = u * (-1 + 2*φ[i, 1]^2)
        @views ruφ₂² = u * (-1 + 2*φ[i, 2]^2)
        @views μ[i, 1] = ruφ₁² * φ[i, 1] - t∇²(φ[:, 1], i) + α * φ[i, 2]
        @views μ[i, 2] = ruφ₂² * φ[i, 2] - t∇²(φ[:, 2], i) - α * φ[i, 1]
    end 
    randn!(rng, ξ)
    ξ .*= σ
    @inbounds for i in 1:N
        @views δφ[i, 1] = ( ∇²(μ[:, 1], i) - ∇(ξ[:, 1], i) ) * dt
        @views δφ[i, 2] = ( ∇²(μ[:, 2], i) - ∇(ξ[:, 2], i) ) * dt
    end
end


function loop!(φt,  φ, μ, δφ, ξ, param_r, step)
    if step=="SO2" euler! = euler_SO2!
    elseif step=="C4" euler! = euler_C4!
    else throw(ErrorException("Non-valid step string given"))
    end

    @showprogress for i in axes(φt, 1)[2:end]
        for _ in 1:skip
            euler!(φ, μ, δφ, ξ, param_r)
            φ .+= δφ
        end
        check(φ, i)
        φt[i,:,:] .= φ
    end
    print('\n')
end


function set_init!(φ0, init)
    x = LinRange(0, L-dx, N)
    if init==1 φ0 .= [ sin.(2π.*x/L)   cos.(2π.*x/L) ]
    elseif init==2 φ0 .= [ .05 .* cos.(2 .* 2π.*x/L)   .5 .* cos.(1 .* 2π.*x/L) ]
    elseif init==3 φ0 .= [ 0 .* x     0.2.* cos.(1 .* 2π.*x/L) ] 
    elseif init==4 φ0 .= [ zeros(N) cat( ones(N÷2), -ones(N÷2), dims=1) / √2 ]
    else φ0 .= zeros(N, 2)
    end
end


function run_euler(param; init=0, name_app="", step="SO2", φ0=fill(NaN, N, 2), write_folder=write_folder)
    u, α, D, bφ1, bφ2 = param
    σ = sqrt(2 * D / dt / dx)
    param_r = (u, α, σ)

    φt = zeros(frames, N, 2)
    φ = zeros(N, 2)
    μ = zeros(N, 2)
    δφ = zeros(N, 2)
    ξ = zeros(N, 2)
    
    if any(isnan.(φ0)) set_init!(φ0, init) end
    # averag value of initial condition should be zero, as we add bφ
    @assert all(abs.(sum(φ0, dims=1)).<1e-10)

    φ = φ0
    φ[:,1] .+= bφ1
    φ[:,2] .+= bφ2
    φt[1,:,:] .= φ

    loop!(φt, φ, μ, δφ, ξ, param_r, step)

    param_write = (u, α, D, bφ1, bφ2, N, L, M*dt, dt)
    
    write_file(φt, param_write; name_app=name_app, write_folder=write_folder)
end


##############
# Utillities #
##############

function write_file(φt, param; name_app=name_app, write_folder=write_folder)
    if ! isdir(write_folder) mkpath(write_folder) end
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
    if any(isnan, φ) throw(ErrorException("Error, NaN detected" )) end
    if (div(i,n)) - div(i-1,n) == 1 && pr print("\r"*"|"^div(i,n)) end
end

