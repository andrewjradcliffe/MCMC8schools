#
# Date created: 2022-02-21
# Author: aradclif
#
#
############################################################################################
#### Example from Handbook of MCMC

function autoreg(ρ, X₁, τ², N::Int)
    x = Vector{Float64}(undef, N)
    τ = √τ²
    x[1] = ρ * X₁ + randn() * τ
    for n = 2:N
        x[n] = ρ * x[n - 1] + randn() * τ
    end
    x
end

τ² = 1.0
ρ = 0.99
n = 10^4

x = autoreg(ρ, 0.0, τ², n);
p = plot(x, label=nothing);
savefig(p, joinpath(pwd(), "autoreg_trace.pdf"))

# Variance of invariant distribution
v² = τ² / (1 - ρ^2)
# Variance in the CLT
σ² = v² * (1 + ρ) / (1 - ρ)
# Approximate error relative to invariant
σ² / v² == (1 + ρ) / (1 - ρ)
# MCMC sample is about as useful as iid sample from invariant marginal with sample size
n / (σ² / v²)

#
μ̂ₙ = mean(x)
σ̂²ₙ = var(x, corrected=false)
# mcse = √σ̂²ₙ / √n
mcse = √σ² / √n
relmcse = mcse / √v²
# To bring relative MCSE to 1%
σ² / v² / (0.01^2)

function acov(x::Vector{T}, k::Int) where {T<:Real}
    μ̂ₙ = mean(x)
    s = zero(promote_type(T, Float64))
    n = length(x)
    for i = 1:(n - k)
        s += (x[i] - μ̂ₙ) * (x[i + k] - μ̂ₙ)
    end
    s / n
end
function acor(x::Vector{T}, k::Int) where {T<:Real}
    μ̂ₙ = mean(x)
    s = zero(promote_type(T, Float64))
    n = length(x)
    for i = 1:(n - k)
        s += (x[i] - μ̂ₙ) * (x[i + k] - μ̂ₙ)
    end
    s / n / varm(x, μ̂ₙ, corrected=false)
end

acf = map(0:500) do k
    acor(x, k)
end;

acf_true = [ρ^k for k = 0:500];

p = bar(acf, label="empirical");
plot!(p, acf_true, label="ρᵏ");
savefig(p, joinpath(pwd(), "autocorrelation.pdf"))

function batches(n::Int, b::Int)
    m, r = divrem(n, b)
    r == 0 || throw(DomainError(b, "b must divide n equally, without remainder"))
    I = Vector{UnitRange{Int}}(undef, m)
    f = 1
    l = b
    for k = 1:m
        I[k] = f:l
        f += b
        l += b
    end
    I
end

function batchmeans(x::Vector{T}, b::Int) where {T<:Real}
    n = length(x)
    m, r = divrem(n, b)
    r == 0 || throw(DomainError(b, "b must divide n equally, without remainder"))
    μ̂ₖ = Vector{Float64}(undef, m)
    for k = 1:m
        s = zero(T)
        f = b * (k - 1) + 1
        l = b * k
        for i = f:l
            s += x[i]
        end
        μ̂ₖ[k] = s / b
    end
    μ̂ₖ
end

# b should be as large as possible, but maintain a reasonable number of batches (≈20-30)
b = 500 # n ÷ 20
μ̂ₖ = batchmeans(x, b);
acf_b = map(0:length(μ̂ₖ)) do k
    acor(μ̂ₖ, k)
end;

p₁ = bar(acf_b, label=nothing, xlabel="lag", ylabel="ACF");
p₂ = scatter(μ̂ₖ, label=nothing, xlabel="Index", ylabel="Batch mean");
p = plot(p₂, p₁);
savefig(p, joinpath(pwd(), "batchmeans.pdf"))
