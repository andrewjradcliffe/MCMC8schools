#
# Date created: 2022-02-08
# Author: aradclif
#
#
############################################################################################
# # Centered parameterization:
# # yⱼ ~ N(αⱼ, σⱼ²)
# # αⱼ ~ N(μ, τ²)
# # μ ~ U(-∞, ∞)
# # τ ~ U(0, ∞)
# # τ transformed to logτ to place on unconstrained space

# function gradlogπ(α, μ, logτ, y, σ²)
#     τ² = abs2(exp(logτ))
#     # ∇l = @. - (α - y) / σ² - (α - μ) / τ²
#     ∇l = .- (α .- y) ./ σ² .- (α .- μ) ./ τ²
#     dπdμ = 0.
#     for j ∈ axes(α, 1)
#         dπdμ -= (μ - α[j]) / τ² # -
#     end
#     dπdlogτ = - length(α) + 1.
#     for j ∈ axes(α, 1)
#         dπdlogτ += abs2(α[j] - μ) / τ²
#     end
#     [∇l; dπdμ; dπdlogτ]
# end
# gradlogπ(q, y, σ²) = gradlogπ(q[1:8], q[9], q[10], y, σ²)

# function logπq(α, μ, logτ, y, σ²)
#     τ² = abs2(exp(logτ))
#     # s = 0
#     # s += logτ
#     s = logτ
#     for j ∈ axes(α, 1)
#         s += -log(√(2π * σ²[j])) - 0.5abs2(y[j] - α[j]) / σ²[j] - 0.5abs2(α[j] - μ) / τ²
#         # s += - 0.5log(2π * σ²[j]) - 0.5abs2(y[j] - α[j]) / σ²[j] - 0.5abs2(α[j] - μ) / τ²
#     end
#     s -= length(α) * log(√(2π * τ²))
#     # s -= length(α) * 0.5log(2π * τ²)
#     s
#     # # unnormalized density
#     # τ² = abs2(exp(logτ))
#     # s = logτ
#     # for j ∈ axes(α, 1)
#     #     s += -0.5abs2(y[j] - α[j]) / σ²[j] - 0.5abs2(α[j] - μ) / τ²
#     # end
#     # s += -0.5log(2π * τ²)
#     # s
# end
# logπq(q, y, σ²) = logπq(q[1:8], q[9], q[10], y, σ²)

# function inversediagdet(M⁻¹::Matrix{T}) where {T<:Real}
#     s = one(T)
#     for j ∈ axes(M⁻¹, 1)
#         s *= one(T) / M⁻¹[j, j]
#     end
#     s
# end

# function weightedss(p, M⁻¹::Matrix{T}) where {T<:Real}
#     s = zero(T)
#     for j ∈ axes(M⁻¹, 1)
#         s += M⁻¹[j, j] * abs2(p[j])
#         # or:
#         # pⱼ = p[j]
#         # s += pⱼ * pⱼ * M⁻¹[j, j]
#     end
#     s
# end

# function logπp(p, M⁻¹)
#     # s = length(p) * log(2π)
#     # s += log(inversediagdet(M⁻¹))
#     # s += weightedss(p, M⁻¹)
#     # s *= -0.5
#     # s
#     # # Or, simpler:
#     -0.5 * (length(p) * log(2π) + log(inversediagdet(M⁻¹)) + weightedss(p, M⁻¹))
#     # # unnormalized density
#     # -0.5weightedss(p, M⁻¹) -0.5log(inversediagdet(M⁻¹))
#     # -0.5weightedss(p, M⁻¹)
# end

# function symplecticintegrate(q₀, p₀, M⁻¹, ϵ, L, y, σ²)
#     # # Option 1
#     # q′, p′ = q₀, p₀
#     # p′ = p′ + 0.5ϵ * gradlogπ(q₀, y, σ²) # -
#     # for l = 2:L
#     #     q′ = q′ + ϵ * M⁻¹ * p′
#     #     p′ = p′ + ϵ * gradlogπ(q′, y, σ²) # -
#     # end
#     # q′ = q′ + ϵ * M⁻¹ * p′
#     # p′ = p′ + 0.5ϵ * gradlogπ(q′, y, σ²) # -
#     # # Option 2
#     p′ = p₀ + 0.5ϵ * gradlogπ(q₀, y, σ²) # -
#     q′ = q₀ + ϵ * M⁻¹ * p′
#     for l = 2:L
#         p′ = p′ + ϵ * gradlogπ(q′, y, σ²) # -
#         q′ = q′ + ϵ * M⁻¹ * p′
#     end
#     p′ = p′ + 0.5ϵ * gradlogπ(q′, y, σ²) # -
#     q′, -p′
# end

# function logratio(q′, p′, q₀, p₀, M⁻¹, y, σ²)
#     logπq(q′, y, σ²) + logπp(p′, M⁻¹) - logπq(q₀, y, σ²) - logπp(p₀, M⁻¹)
# end

# function metropolis(q′, p′, q₀, p₀, M⁻¹, y, σ²)
#     r = exp(logratio(q′, p′, q₀, p₀, M⁻¹, y, σ²))
#     rand() < r ? (q′, p′) : (q₀, p₀)
# end

# # Using the mass matrix directly
# function momentumrand_M(M)
#     p = randn(size(M, 1))
#     @inbounds for j ∈ axes(M, 1)
#         p[j] *= √(M[j, j])
#     end
#     p
# end

# # Using the lower triangular of the cholesky, AAᵀ = M
# function momentumrand_A(A)
#     p = randn(size(A, 1))
#     @inbounds for j ∈ axes(A, 1)
#         p[j] *= A[j, j]
#     end
#     p
# end

# # Using the inverse mass matrix directly
# function momentumrand_M⁻¹(M⁻¹)
#     p = randn(size(M⁻¹, 1))
#     @inbounds for j ∈ axes(M⁻¹, 1)
#         p[j] *= inv(√(M⁻¹[j, j]))
#     end
#     p
# end

# function hmc(ϵ, L, M, M⁻¹, α₀, μ₀, logτ₀, y, σ²)
#     q₀ = [α₀; μ₀; logτ₀]
#     p₀ = momentumrand_M(M)
#     q′, p′ = symplecticintegrate(q₀, p₀, M⁻¹, ϵ, L, y, σ²)
#     metropolis(q′, p′, q₀, p₀, M⁻¹, y, σ²)
# end

# function hmc(N, ϵ, L, M, M⁻¹, α₀, μ₀, logτ₀, y, σ²)
#     q₀ = [α₀; μ₀; logτ₀]
#     # θ = Vector{Vector{Float64}}()
#     θ = Matrix{Float64}(undef, 10, N)
#     n, t = 0, 0
#     while n < N
#         p₀ = momentumrand_M(M)
#         q′, p′ = symplecticintegrate(q₀, p₀, M⁻¹, ϵ, L, y, σ²)
#         r = exp(logπq(q′, y, σ²) + logπp(p′, M⁻¹) - logπq(q₀, y, σ²) - logπp(p₀, M⁻¹))
#         if rand() < r
#             # push!(θ, q′)
#             q₀ = q′
#             n += 1
#             t += 1
#             θ[:, n] = q′
#         else
#             t += 1
#         end
#     end
#     θ, n, t
# end

# function hmc(N, ϵ, L, M, M⁻¹, y, σ²)
#     J = length(y)
#     q₀, t = initial(ϵ, L, M, M⁻¹, y, σ²)
#     α₀, μ₀, logτ₀ = q₀[1:J], q₀[J + 1], q₀[J + 2]
#     hmc(N, ϵ, L, M, M⁻¹, α₀, μ₀, logτ₀, y, σ²)
# end

# function initial(ϵ, L, M, M⁻¹, y, σ²)
#     J = length(y)
#     𝓈 = sqrt.(diag(M⁻¹))
#     t = 0
#     while true
#         α₀ = randn(J) .* 𝓈[1:J]
#         μ₀ = randn() * 𝓈[J + 1]
#         logτ₀ = randn() * 𝓈[J + 2]
#         q₀ = [α₀; μ₀; logτ₀]
#         q′, p′ = hmc(ϵ, L, M, M⁻¹, α₀, μ₀, logτ₀, y, σ²)
#         t += 1
#         q′ != q₀ && break
#     end
#     return q′, t
# end

# ################ Data
# y = Float64[28, 8, -3, 7, -1, 1, 18, 12];
# σ = Float64[15, 10, 16, 11, 9, 11, 10, 18];
# σ² = abs2.(σ);
# ################ Putting it all together
# # Initialize
# J = length(y)
# # # Stan's initialization -- far too restrictive for this model
# # α₀ = rand(J) .* 4 .- 2
# # μ₀ = rand() * 4 - 2
# # logτ₀ = log(rand() * 15)
# # # Default initialization
# scale = 15;
# # α₀ = randn(J) .* scale;
# # μ₀ = randn() * scale;
# # logτ₀ = randn();
# # Specify ϵ, L, N
# ϵ = 0.1
# L = 10
# N = 1000
# #
# # M = diagm(scale*(ones(10))); M[10, 10] = 1;
# # M⁻¹ = diagm(1 ./ diag(M));
# # Approximate M: inverse of covariance matrix of posterior
# M = diagm(fill(1 / scale^2, 10)); M[10, 10] = 1;
# M⁻¹ = diagm(1 ./ diag(M));

# # Initialization tests -- guarantee a decent starting point
# q₀, t = initial(ϵ, L, M, M⁻¹, y, σ²)
# α₀, μ₀, logτ₀ = q₀[1:J], q₀[J + 1], q₀[J + 2]
# ####

# # Single step, fixed initialization
# @timev q, p = hmc(ϵ, L, M, M⁻¹, α₀, μ₀, logτ₀, y, σ²)
# # Fixed initialization
# @timev θ, n, t = hmc(N, ϵ, L, M, M⁻¹, α₀, μ₀, logτ₀, y, σ²);

# # Auto-initialization
# @timev θ, n, t = hmc(N, ϵ, L, M, M⁻¹, y, σ²);
# acprob = n / t
# mean(θ, dims=2)
# mean(θ[:, (N >> 1):N], dims=2)

# # Multiple chains
# @timev chains = [hmc(N, ϵ, L, M, M⁻¹, y, σ²) for _ = 1:4];
# acprob = getindex.(chains, 2) ./ getindex.(chains, 3)
# θ = cat(first.(chains)..., dims=3);
# mean(θ, dims=(2,3))

# # Straightforward tests
# gradlogπ(α₀, μ₀, logτ₀, y, σ²)
# @benchmark momentumrand_M(M)
# A, Aᵀ = cholesky(M)
# momentumrand_M⁻¹(M⁻¹)
# @benchmark momentumrand_A(A)
# A * randn(10)
# inversediagdet(M⁻¹)
# @benchmark weightedss(p, M⁻¹)
# @benchmark weightedss2(p, M⁻¹)

# # Trajectory test
# q₀, p₀ = [α₀; μ₀; logτ₀], momentumrand_M(M);
# q′, p′ = symplecticintegrate(q₀, p₀, M⁻¹, ϵ, L, y, σ²)
# logπq(q′, y, σ²)
# logπp(p′, M⁻¹)
# logratio(q′, p′, q₀, p₀, M⁻¹, y, σ²)


# # normpdf(x) = inv(√(2π)) * exp(-0.5 * x * x)
# # lognormpdf(x) = -0.5log(2π) - 0.5 * x * x
# # dlognormpdf(x) = -x

# ############################################################################################
# #### 2022-02-14: non-centered parameterization
# # yⱼ ~ N(μ + τ * ηⱼ, σⱼ²)
# # ηⱼ ~ N(0, 1)
# # μ ~ U(-∞, ∞)
# # τ ~ U(0, ∞)
# # τ transformed to logτ to place on unconstrained space

# # Clearly not the most efficient approach -- the proper way is to compute a vector, the
# # elements of which are (y[j] - μ - τ * η[j]). Then, compute dπdμ and dπdlogτ, then
# # mutate in place to complete the computation of ∇l.
# function gradlogπ_nc(η, μ, logτ, y, σ²)
#     τ = exp(logτ)
#     ∇l = Vector{Float64}(undef, length(η))
#     for j ∈ eachindex(η, y, σ²)
#         ∇l[j] = (y[j] - μ - τ * η[j]) * τ / σ²[j] - η[j]
#     end
#     dπdμ = 0.0
#     for j ∈ eachindex(η, y, σ²)
#         dπdμ += (y[j] - μ - τ * η[j]) / σ²[j]
#     end
#     dπdlogτ = 1.0
#     for j ∈ eachindex(η, y, σ²)
#         dπdlogτ += τ * (y[j] - μ - τ * η[j]) * η[j] / σ²[j]
#     end
#     [∇l; dπdμ; dπdlogτ]
# end
# gradlogπ_nc(q, y, σ²) = gradlogπ_nc(q[1:8], q[9], q[10], y, σ²)

# function logπq_nc(η, μ, logτ, y, σ²)
#     τ = exp(logτ)
#     s = logτ
#     for j ∈ eachindex(η, y, σ²)
#         s += -0.5(log(2π * σ²[j]) + (1.0 / σ²[j]) * abs2(y[j] - μ - τ * η[j]) + log(2π) + abs2(η[j]))
#     end
#     s
# end
# logπq_nc(q, y, σ²) = logπq_nc(q[1:8], q[9], q[10], y, σ²)

# function symplecticintegrate(f_gradlogπ, q₀, p₀, M⁻¹, ϵ, L, y, σ²)
#     p′ = p₀ + 0.5ϵ * f_gradlogπ(q₀, y, σ²)
#     q′ = q₀ + ϵ * M⁻¹ * p′
#     for l = 2:L
#         p′ = p′ + ϵ * f_gradlogπ(q′, y, σ²)
#         q′ = q′ + ϵ * M⁻¹ * p′
#     end
#     p′ = p′ + 0.5ϵ * f_gradlogπ(q′, y, σ²)
#     q′, -p′
# end

# function targetrand(M⁻¹::Matrix{T}) where {T<:Real}
#     K = size(M⁻¹, 1)
#     q₀ = randn(K)
#     for k ∈ eachindex(q₀)
#         q₀[k] *= √(M⁻¹[k, k])
#     end
#     q₀
# end

# function hmc_nc_transition(f_gradlogπ, f_logπq, ϵ, L, M, M⁻¹, q₀, y, σ²)
#     p₀ = momentumrand_M(M)
#     q′, p′ = symplecticintegrate(f_gradlogπ, q₀, p₀, M⁻¹, ϵ, L, y, σ²)
#     r = exp(f_logπq(q′, y, σ²) + logπp(p′, M⁻¹) - f_logπq(q₀, y, σ²) - logπp(p₀, M⁻¹))
#     u = rand()
#     u < r ? (q′, p′) : (q₀, p₀)
# end
# function hmc_nc_transition(f_gradlogπ, f_logπq, ϵ, L, M, M⁻¹, y, σ²)
#     q₀ = targetrand(M⁻¹)
#     hmc_nc_transition(f_gradlogπ, f_logπq, ϵ, L, M, M⁻¹, q₀, y, σ²)
# end

# function hmc_nc(f_gradlogπ, f_logπq, N, ϵ, L, M, M⁻¹, q₀, y, σ²)
#     θ = Matrix{Float64}(undef, length(q₀), N)
#     n, t = 0, 0
#     while n < N
#         p₀ = momentumrand_M(M)
#         q′, p′ = symplecticintegrate(f_gradlogπ, q₀, p₀, M⁻¹, ϵ, L, y, σ²)
#         r = exp(f_logπq(q′, y, σ²) + logπp(p′, M⁻¹) - f_logπq(q₀, y, σ²) - logπp(p₀, M⁻¹))
#         if rand() < r
#             q₀ = q′
#             n += 1
#             t += 1
#             θ[:, n] = q′
#         else
#             t += 1
#         end
#     end
#     θ, n, t
# end

# function hmc_nc(f_gradlogπ, f_logπq, N, ϵ, L, M, M⁻¹, y, σ²)
#     q₀ = targetrand(M⁻¹)
#     hmc_nc(f_gradlogπ, f_logπq, N, ϵ, L, M, M⁻¹, q₀, y, σ²)
# end

# # #### Unlike the centered parameterization, the non-centered parameterization
# # will not get stuck at low values of τ, hence, no need to guarantee a "decent" starting point.
# # function initial_nc(f_gradlogπ, f_logπq, ϵ, L, M, M⁻¹, y, σ²)
# #     K = size(M⁻¹, 1)
# #     𝓈 = sqrt.(diag(M⁻¹))
# #     t = 0
# #     while true
# #         q₀ = randn(K) .* 𝓈
# #         q′, p′ = hmc_nc_1iter(f_gradlogπ, f_logπq, ϵ, L, M, M⁻¹, q₀, y, σ²)
# #         t += 1
# #         q′ != q₀ && return q′, t
# #     end
# #     return q′, t
# # end
# # # Initialization tests -- guarantee a decent starting point
# # q₀, t = initial_nc(gradlogπ_nc, logπq_nc, ϵ, L, M, M⁻¹, y, σ²)

# # Data
# y = Float64[28, 8, -3, 7, -1, 1, 18, 12];
# σ = Float64[15, 10, 16, 11, 9, 11, 10, 18];
# σ² = abs2.(σ);

# # Integration time and stepsize. Interestingly, after re-parameterization,
# # ϵ = 0.2, L = 5 gives a reasonable acceptance rate. ϵ = 0.1, L = 10 results in acceptance
# # rate of ≥ 0.99, which indicates the integration could take larger steps.
# ϵ = 0.2
# L = 5
# N = 1000

# # Euclidean-Gaussian metric
# scale = 15;
# M⁻¹ = diagm(ones(10)); M⁻¹[9, 9] = scale^2;
# M = diagm(1 ./ diag(M⁻¹));

# # Single step, random initialization
# q₀, p₀ = hmc_nc_transition(gradlogπ_nc, logπq_nc, ϵ, L, M, M⁻¹, y, σ²)
# # Single step, fixed initialization
# @timev q, p = hmc_nc_transition(gradlogπ_nc, logπq_nc, ϵ, L, M, M⁻¹, q₀, y, σ²)
# # Fixed initialization
# @timev θ, n, t = hmc_nc(gradlogπ_nc, logπq_nc, N, ϵ, L, M, M⁻¹, q₀, y, σ²);

# # Auto-initialization
# @timev θ, n, t = hmc_nc(gradlogπ_nc, logπq_nc, N, ϵ, L, M, M⁻¹, y, σ²);
# acprob = n / t
# mean(θ, dims=2)

# ϑ = originalparam(θ);
# mean(ϑ, dims=2)

# # Multiple chains
# @timev chains = [hmc_nc(gradlogπ_nc, logπq_nc, N, ϵ, L, M, M⁻¹, y, σ²) for _ = 1:4];
# acprob = getindex.(chains, 2) ./ getindex.(chains, 3)
# θ = cat(first.(chains)..., dims=3);
# mean(θ, dims=(2,3))

# ϑ = originalparam(θ);
# mean(ϑ, dims=(2, 3))
# # Convergence check
# neff(ϑ)
# Rhat(ϑ)

# # Trajectory tests
# p₀ = momentumrand_M(M)
# q, p = symplecticintegrate(gradlogπ_nc, q₀, p₀, M⁻¹, ϵ, L, y, σ²)
# #

# # # Experiment with NamedTuple -- no change in allocations or speed
# # 𝒹 = (; :y => y, :σ² => σ²)
# # logπq_nc(q, 𝒹::NamedTuple) = logπq_nc(q[1:8], q[9], q[10], 𝒹.y, 𝒹.σ²)
# # @timev logπq_nc(q₀, y, σ²)
# # @timev logπq_nc(q₀, 𝒹)

# function originalparam(θ::Vector{T}) where {T<:Real}
#     K = length(θ)
#     τ = exp(θ[K])
#     μ = θ[K - 1]
#     ϑ = Vector{promote_type(T, Float64)}(undef, K)
#     for k = 1:(K - 2)
#         ϑ[k] = μ + τ * θ[k]
#     end
#     ϑ[K - 1] = μ
#     ϑ[K] = τ
#     ϑ
# end

# function originalparam(θ::Matrix{T}) where {T<:Real}
#     K, S = size(θ)
#     ϑ = Matrix{promote_type(T, Float64)}(undef, K, S)
#     for s ∈ axes(θ, 2)
#         τ = exp(θ[K, s])
#         μ = θ[K - 1, s]
#         for k = 1:(K - 2)
#             ϑ[k, s] = μ + τ * θ[k, s]
#         end
#         ϑ[K - 1, s] = μ
#         ϑ[K, s] = τ
#     end
#     ϑ
# end
# originalparam(θ::Array{T, 3}) where {T<:Real} = mapslices(originalparam, θ, dims=(1, 2))


# ################################################################
# # Functions to support formal version, other than those provided above
# function symplecticintegrate(f_gradlogπ, q₀, p₀, M⁻¹, ϵ, L)
#     p′ = p₀ + 0.5ϵ * f_gradlogπ(q₀)
#     q′ = q₀ + ϵ * M⁻¹ * p′
#     for l = 2:L
#         p′ = p′ + ϵ * f_gradlogπ(q′)
#         q′ = q′ + ϵ * M⁻¹ * p′
#     end
#     p′ = p′ + 0.5ϵ * f_gradlogπ(q′)
#     q′, -p′
# end

# ################ Thinking in terms of the Hamiltonian
# H(f_logπq, q, p, M⁻¹) = -f_logπq(q) - logπp(p, M⁻¹)
# E_pq(f_logπq, q, E) = E + f_logπq(q)
# logπp(f_logπq, q, E) = E_pq(f_logπq, q, E)

# function transition(f_logπq, f_gradlogπ, q₀, p₀, M, M⁻¹, ϵ, L)
#     q′, p′ = symplecticintegrate(f_gradlogπ, q₀, p₀, M⁻¹, ϵ, L)
#     r = exp(-H(f_logπq, q′, p′, M⁻¹) + H(f_logπq, q₀, p₀, M⁻¹))
#     u = rand()
#     u < r ? (q′, p′) : (q₀, p₀)
# end
# transition(f_logπq, f_gradlogπ, q₀, M, M⁻¹, ϵ, L) =
#     transition(f_logπq, f_gradlogπ, q₀, momentumrand_M(M), M, M⁻¹, ϵ, L)
# transition(f_logπq, f_gradlogπ, M, M⁻¹, ϵ, L) =
#     transition(f_logπq, f_gradlogπ, targetrand(M⁻¹), momentumrand_M(M), M, M⁻¹, ϵ, L)

# function hmc(f_logπq, f_gradlogπ, q₀, M, M⁻¹, ϵ, L, N)
#     θ = Matrix{Float64}(undef, size(M, 1), N)
#     E = Vector{Float64}(undef, N)
#     n, t = 0, 0
#     while n < N
#         p₀ = momentumrand_M(M)
#         q′, p′ = symplecticintegrate(f_gradlogπ, q₀, p₀, M⁻¹, ϵ, L)
#         E₀ = H(f_logπq, q₀, p₀, M⁻¹)
#         E′ = H(f_logπq, q′, p′, M⁻¹)
#         r = exp(-E′ + E₀)
#         if rand() < r
#             n += 1
#             θ[:, n] = q′
#             q₀ = q′
#             E[n] = E′
#         end
#         t += 1
#     end
#     θ, E, n, t
# end
# hmc(f_logπq, f_gradlogπ, M, M⁻¹, ϵ, L, N) =
#     hmc(f_logπq, f_gradlogπ, targetrand(M⁻¹), M, M⁻¹, ϵ, L, N)

# function ebfmi(E::Vector{T}) where {T<:Real}
#     N = length(E)
#     s = zero(T)
#     # As written, this has _clear_ potential for numerical stability issues
#     for n = 2:N
#         s += abs2(E[n] - E[n - 1])
#     end
#     s / (N * var(E, corrected=false))
# end

# function microcanonicalenergy(f_logπq, q, E)
#     Eq = -f_logπq(q)
#     Epq = E - Epq
#     Epq, Eq
# end

# function microcanonicalenergy(f_logπq, Q::Matrix{T}, E::Vector{T}) where {T<:Real}
#     N = length(E)
#     Eq = Vector{Float64}(undef, N)
#     Epq = Vector{Float64}(undef, N)
#     for (n, q) ∈ enumerate(eachcol(Q))
#         Eq[n] = -f_logπq(q)
#         Epq[n] = E[n] - Eq[n]
#     end
#     Epq, Eq
# end

# # Guaranteed initialization for models in which small τ causes sticking
# # (i.e. for centered parameterizations)
# function initialize(f_logπq, f_gradlogπ, M, M⁻¹, ϵ, L)
#     t = 0
#     while true
#         q₀, p₀ = targetrand(M⁻¹), momentumrand_M(M)
#         q′, p′ = transition(f_logπq, f_gradlogπ, q₀, p₀, M, M⁻¹, ϵ, L)
#         t += 1
#         q′ != q₀ && return q′, p′, t
#     end
#     return q′, p′, t
# end

# # Data
# y = Float64[28, 8, -3, 7, -1, 1, 18, 12];
# σ = Float64[15, 10, 16, 11, 9, 11, 10, 18];
# σ² = abs2.(σ);

# # Integration time and stepsize.
# ϵ = 0.2
# L = 5
# N = 1000

# # Euclidean-Gaussian metric
# scale = 15;
# M⁻¹ = diagm(ones(10)); M⁻¹[9, 9] = scale^2;
# M = diagm(1 ./ diag(M⁻¹));

# f_logπ = let y = y, σ² = σ²
#     q -> logπq_nc(q, y, σ²)
# end;
# f_gradlogπ = let y = y, σ² = σ²
#     q -> gradlogπ_nc(q, y, σ²)
# end;

# @benchmark transition(f_logπ, f_gradlogπ, M, M⁻¹, ϵ, L)
# @timev θ, E, n, t = hmc(f_logπ, f_gradlogπ, M, M⁻¹, ϵ, L, N);
# @benchmark hmc(f_logπ, f_gradlogπ, M, M⁻¹, ϵ, L, N)

# ebfmi(E)

# Epq, Eq = microcanonicalenergy(f_logπ, θ, E);
# pₕ = histogram([Epq Eq] .- [mean(Epq) mean(Eq)], labels=["E_pq ≡ π(E|q)" "E_q ≡ π(E)"], alpha=0.5,
#                annotations=((0.1, 1.0), text("E-BFMI = $(round(ebfmi(E), digits=3))", 8)),
#                xlabel="E - <E>");
# savefig(pₕ, joinpath(pwd(), "energy.pdf"))

# # Multiple chains
# @timev chains = [hmc(f_logπ, f_gradlogπ, M, M⁻¹, ϵ, L, N) for _ = 1:4];
# acprob = getindex.(chains, 3) ./ getindex.(chains, 4)
# θ = cat(first.(chains)..., dims=3);
# mean(θ, dims=(2,3))
# E = hcat(getindex.(chains, 2)...);

# ϑ = originalparam(θ);
# mean(ϑ, dims=(2, 3))
# # Convergence check
# neff(ϑ)
# Rhat(ϑ)

# ps = map(1:4) do k
#     Epq, Eq = microcanonicalenergy(f_logπ, θ[:, :, k], E[:, k])
#     pₕ = histogram([Epq Eq] .- [mean(Epq) mean(Eq)], labels=["E_pq ≡ π(E|q)" "E_q ≡ π(E)"],
#                    alpha=0.5,
#                    annotations=((0.1, 1.0), text("E-BFMI = $(round(ebfmi(E[:, k]), digits=3))", 8)),
#                    xlabel="E - <E>")
#     pₕ
# end;
# pₕ = plot(ps...);
# savefig(pₕ, joinpath(pwd(), "energy_$(N_chains)chain.pdf"))

# #### Centered parameterization -- might get stuck due to small τ
# ϵ = 0.1
# L = 10
# N = 1000
# N_chains = 4

# scale = 15;
# M = diagm(fill(1 / scale^2, 10)); M[10, 10] = 1;
# M⁻¹ = diagm(1 ./ diag(M));

# f_logπ_c = let y = y, σ² = σ²
#     q -> logπq(q, y, σ²)
# end;
# f_gradlogπ_c = let y = y, σ² = σ²
#     q -> gradlogπ(q, y, σ²)
# end;


# qq, pp = transition(f_logπ_c, f_gradlogπ_c, M, M⁻¹, ϵ, L)
# # Typically, centered parameterization necessitates a non-small τ starting point
# qᵢ, pᵢ, tᵢ = initialize(f_logπ_c, f_gradlogπ_c, M, M⁻¹, ϵ, L)
# @timev θ, E, n, t = hmc(f_logπ_c, f_gradlogπ_c, qᵢ, M, M⁻¹, ϵ, L, N);


# Epq, Eq = microcanonicalenergy(f_logπ_c, θ, E);
# pₕ = histogram([Epq Eq] .- [mean(Epq) mean(Eq)], labels=["E_pq ≡ π(E|q)" "E_q ≡ π(E)"], alpha=0.5,
#                annotations=((0.1, 1.0), text("E-BFMI = $(round(ebfmi(E), digits=3))", 8)),
#                xlabel="E - <E>");
# savefig(pₕ, joinpath(pwd(), "energy_c.pdf"))

# # Multiple chains
# @timev chains = [hmc(f_logπ_c, f_gradlogπ_c,
#                      first(initialize(f_logπ_c, f_gradlogπ_c, M, M⁻¹, ϵ, L)), # qᵢ's
#                      M, M⁻¹, ϵ, L, N) for _ = 1:N_chains];
# acprob = getindex.(chains, 3) ./ getindex.(chains, N_chains)
# θ = cat(first.(chains)..., dims=3);
# mean(θ, dims=(2,3))
# E = hcat(getindex.(chains, 2)...);

# # Convergence check
# neff(θ)
# Rhat(θ)

# ps = map(1:4) do k
#     Epq, Eq = microcanonicalenergy(f_logπ_c, θ[:, :, k], E[:, k])
#     pₕ = histogram([Epq Eq] .- [mean(Epq) mean(Eq)], labels=["E_pq ≡ π(E|q)" "E_q ≡ π(E)"],
#                    alpha=0.5,
#                    annotations=((0.1, 1.0), text("E-BFMI = $(round(ebfmi(E[:, k]), digits=3))", 8)),
#                    xlabel="E - <E>")
#     pₕ
# end;
# pₕ = plot(ps...);
# savefig(pₕ, joinpath(pwd(), "energy_c_$(N_chains)chain.pdf"))
# ############################################################################################
# #### MCMC diagnostics: potential scale reduction and effective sample size
# # Simulations arranged on dimension 2, chains on dimension 3
# # function withinchainmean(chain::AbstractMatrix)
# #     ψ̄ = mean(chain, dims=2)
# # end
# function betweenchain(ϑ::Array{T, 3}) where {T<:Real}
#     B = size(ϑ, 2) .* var(dropdims(mean(ϑ, dims=2), dims=2), dims=2)
# end

# function withinchain(ϑ::Array{T, 3}) where {T<:Real}
#     W = mean(dropdims(var(ϑ, dims=2), dims=2), dims=2)
# end

# function vhat(ψ::AbstractMatrix{T}) where {T<:Real}
#     n = size(ψ, 1)
#     B = n * var(mean(ψ, dims=1))
#     W = mean(var(ψ, dims=1))
#     v̂⁺ = ((n - 1) / n) * W + (1 / n) * B
# end


# function vhat(ϑ::Array{T, 3}) where {T<:Real}
#     n = size(ϑ, 2)
#     B = n .* var(dropdims(mean(ϑ, dims=2), dims=2), dims=2)
#     W = mean(dropdims(var(ϑ, dims=2), dims=2), dims=2)
#     v̂⁺ = ((n - 1) / n) .* W .+ (1 / n) .* B
# end

# function Rhat(ψ::AbstractMatrix{T}) where {T<:Real}
#     n = size(ψ, 1)
#     B = n * var(mean(ψ, dims=1))
#     W = mean(var(ψ, dims=1))
#     v̂⁺ = ((n - 1) / n) * W + (1 / n) * B
#     R̂ = √(v̂⁺ / W)
# end

# function Rhat(ϑ::Array{T, 3}) where {T<:Real}
#     n = size(ϑ, 2)
#     B = n .* var(dropdims(mean(ϑ, dims=2), dims=2), dims=2)
#     W = mean(dropdims(var(ϑ, dims=2), dims=2), dims=2)
#     v̂⁺ = ((n - 1) / n) .* W .+ (1 / n) .* B
#     R̂ = sqrt.(v̂⁺ ./ W)
# end


# ϑ = θ[:, (N >> 1 + 1):N, :]; #view(θ, :, (N >> 1 + 1):N, :);
# L = size(ϑ, 2)
# n = L >> 1
# ϑ₂ = [view(ϑ, :, 1:n, :);;; view(ϑ, :, (n + 1):L, :)];

# ψ = ϑ₂[1, :, :];

# B = betweenchain(ϑ)
# W = withinchain(ϑ)
# v̂⁺ = ((L - 1) / L) .* W .+ (1 / L) .* B
# R̂ = sqrt.(v̂⁺ ./ W)

# Rhat(ϑ)
# Rhat(ϑ₂)

# function variogram(ψ::AbstractMatrix{T}, t::Int) where {T<:Real}
#     n, m = size(ψ)
#     s = zero(T)
#     for j ∈ axes(ψ, 2)
#         for i = (t + 1):n
#             s += abs2(ψ[i, j] - ψ[i - t, j])
#         end
#     end
#     return s / (m * (n - t))
# end

# autocor(ψ::AbstractMatrix{T}, t::Int) where {T<:Real} = one(T) - variogram(ψ, t) / 2vhat(ψ)

# function findautolag(ψ::AbstractMatrix{T}) where {T<:Real}
#     n = size(ψ, 1)
#     for t = 1:2:(n - 3)
#         ρ̂ₜ₊₁ = autocor(ψ, t)
#         ρ̂ₜ₊₂ = autocor(ψ, t + 1)
#         ρ̂ₜ₊₁ + ρ̂ₜ₊₂ < 0 && return t
#     end
#     return n - 3
# end

# function neff(ψ::AbstractMatrix{T}) where {T<:Real}
#     n, m = size(ψ)
#     # t = findautolag(ψ)
#     Σρ = zero(T)
#     for t = 1:2:(n - 3)
#         ρ̂ₜ₊₁ = autocor(ψ, t)
#         ρ̂ₜ₊₂ = autocor(ψ, t + 1)
#         s = ρ̂ₜ₊₁ + ρ̂ₜ₊₂
#         s < 0 && break
#         Σρ += s
#         # s < 0 ? break : (Σρ += s)
#     end
#     n̂ₑ = n * m / (one(T) + 2Σρ)
# end

# neff(ϑ::Array{T, 3}) where {T<:Real} = mapslices(neff, ϑ, dims=(2, 3))
# neff(ϑ)
# Rhat(ϑ)

# variogram(ψ, 1)
# findautolag(ψ)
# @benchmark neff(ψ)
# mapslices(neff, c̃, dims=(2, 3))
# mapslices(Rhat, c̃, dims=(2, 3))

################################################################
# Import Turing and DynamicHMC.
using DynamicHMC, Turing

# Model definition.
@model function gdemo(x, y)
    s² ~ InverseGamma(2, 3)
    m ~ Normal(0, sqrt(s²))
    x ~ Normal(m, sqrt(s²))
    y ~ Normal(m, sqrt(s²))
end

# Pull 2,000 samples using DynamicNUTS.
chn = sample(gdemo(1.5, 2.0), DynamicNUTS(), 2000)

#### The classic 8-schools
@model function sch8(y, σ)
    τ ~ Uniform(0, 225)
    μ ~ Uniform(-225, 225)
    k = length(y)
    α ~ filldist(Normal(μ, τ), k)
    y .~ Normal.(α, σ)
end
@timev chn2 = sample(sch8(y, σ), DynamicNUTS(), MCMCThreads(), 1000, 4)
θ2 = chn2.value.data
mean(θ2, dims=1)
chn2[:α]
chn3 = sample(sch8(y, σ), HMC(ϵ, L), 1000)

############################################################################################
using LinearAlgebra, Statistics
############################################################################################
#General utilities
function diagdet(M::Matrix{T}) where {T<:Real}
    s = one(T)
    for j ∈ axes(M, 1)
        s *= M[j, j]
    end
    s
end

function inversediagdet(M⁻¹::Matrix{T}) where {T<:Real}
    s = one(T)
    for j ∈ axes(M⁻¹, 1)
        s *= one(T) / M⁻¹[j, j]
    end
    s
end

function weightedss(p, M⁻¹::Matrix{T}) where {T<:Real}
    s = zero(T)
    for j ∈ axes(M⁻¹, 1)
        s += M⁻¹[j, j] * abs2(p[j])
        # or:
        # pⱼ = p[j]
        # s += pⱼ * pⱼ * M⁻¹[j, j]
    end
    s
end

function diagweightedss(x, μ, M::Matrix{T}) where {T<:Real}
    s = zero(T)
    for j ∈ axes(M, 1)
        s += inv(M[j, j]) * abs2(x[j] - μ[j])
    end
    s
end

function logπp(p, M⁻¹)
    # s = length(p) * log(2π)
    # s += log(inversediagdet(M⁻¹))
    # s += weightedss(p, M⁻¹)
    # s *= -0.5
    # s
    # # Or, simpler:
    -0.5 * (length(p) * log(2π) + log(inversediagdet(M⁻¹)) + weightedss(p, M⁻¹))
    # # unnormalized density
    # -0.5weightedss(p, M⁻¹) -0.5log(inversediagdet(M⁻¹))
    # -0.5weightedss(p, M⁻¹)
end

# Using the mass matrix directly
function momentumrand_M(M)
    p = randn(size(M, 1))
    @inbounds for j ∈ axes(M, 1)
        p[j] *= √(M[j, j])
    end
    p
end

# Using the lower triangular of the cholesky, AAᵀ = M
function momentumrand_A(A)
    p = randn(size(A, 1))
    @inbounds for j ∈ axes(A, 1)
        p[j] *= A[j, j]
    end
    p
end

# Using the inverse mass matrix directly
function momentumrand_M⁻¹(M⁻¹)
    p = randn(size(M⁻¹, 1))
    @inbounds for j ∈ axes(M⁻¹, 1)
        p[j] *= inv(√(M⁻¹[j, j]))
    end
    p
end

function targetrand(M⁻¹::Matrix{T}) where {T<:Real}
    K = size(M⁻¹, 1)
    q₀ = randn(K)
    for k ∈ eachindex(q₀)
        q₀[k] *= √(M⁻¹[k, k])
    end
    q₀
end
############################################################################################
# Static HMC implementation:
# This is a high-level implementation, not as efficient as I would normally write.
# Much, so much in-place mutation can be used to minimize memory usage...
################################################################
# Functions to support formal version, other than those provided above
function symplecticintegrate(f_gradlogπ, q₀, p₀, M⁻¹, ϵ, L)
    p′ = p₀ + 0.5ϵ * f_gradlogπ(q₀)
    q′ = q₀ + ϵ * M⁻¹ * p′
    for l = 2:L
        p′ = p′ + ϵ * f_gradlogπ(q′)
        q′ = q′ + ϵ * M⁻¹ * p′
    end
    p′ = p′ + 0.5ϵ * f_gradlogπ(q′)
    q′, -p′
end

################ Thinking in terms of the Hamiltonian
H(f_logπq, q, p, M⁻¹) = -f_logπq(q) - logπp(p, M⁻¹)
E_pq(f_logπq, q, E) = E + f_logπq(q)
logπp(f_logπq, q, E) = E_pq(f_logπq, q, E)

function transition(f_logπq, f_gradlogπ, q₀, p₀, M, M⁻¹, ϵ, L)
    q′, p′ = symplecticintegrate(f_gradlogπ, q₀, p₀, M⁻¹, ϵ, L)
    r = exp(-H(f_logπq, q′, p′, M⁻¹) + H(f_logπq, q₀, p₀, M⁻¹))
    u = rand()
    u < r ? (q′, p′) : (q₀, p₀)
end
transition(f_logπq, f_gradlogπ, q₀, M, M⁻¹, ϵ, L) =
    transition(f_logπq, f_gradlogπ, q₀, momentumrand_M(M), M, M⁻¹, ϵ, L)
transition(f_logπq, f_gradlogπ, M, M⁻¹, ϵ, L) =
    transition(f_logπq, f_gradlogπ, targetrand(M⁻¹), momentumrand_M(M), M, M⁻¹, ϵ, L)

function hmc(f_logπq, f_gradlogπ, q₀, M, M⁻¹, ϵ, L, N)
    θ = Matrix{Float64}(undef, size(M, 1), N)
    E = Vector{Float64}(undef, N)
    n, t = 0, 0
    while n < N
        p₀ = momentumrand_M(M)
        q′, p′ = symplecticintegrate(f_gradlogπ, q₀, p₀, M⁻¹, ϵ, L)
        E₀ = H(f_logπq, q₀, p₀, M⁻¹)
        E′ = H(f_logπq, q′, p′, M⁻¹)
        r = exp(-E′ + E₀)
        if rand() < r
            n += 1
            θ[:, n] = q′
            q₀ = q′
            E[n] = E′
        end
        t += 1
    end
    θ, E, n, t
end
hmc(f_logπq, f_gradlogπ, M, M⁻¹, ϵ, L, N) =
    hmc(f_logπq, f_gradlogπ, targetrand(M⁻¹), M, M⁻¹, ϵ, L, N)

function ebfmi(E::Vector{T}) where {T<:Real}
    N = length(E)
    s = zero(T)
    # As written, this has _clear_ potential for numerical stability issues
    for n = 2:N
        s += abs2(E[n] - E[n - 1])
    end
    s / (N * var(E, corrected=false))
end

function microcanonicalenergy(f_logπq, q, E)
    Eq = -f_logπq(q)
    Epq = E - Epq
    Epq, Eq
end

function microcanonicalenergy(f_logπq, Q::Matrix{T}, E::Vector{T}) where {T<:Real}
    N = length(E)
    Eq = Vector{Float64}(undef, N)
    Epq = Vector{Float64}(undef, N)
    for (n, q) ∈ enumerate(eachcol(Q))
        Eq[n] = -f_logπq(q)
        Epq[n] = E[n] - Eq[n]
    end
    Epq, Eq
end

# Guaranteed initialization for models in which small τ causes sticking
# (i.e. for centered parameterizations)
function initialize(f_logπq, f_gradlogπ, M, M⁻¹, ϵ, L)
    t = 0
    while true
        q₀, p₀ = targetrand(M⁻¹), momentumrand_M(M)
        q′, p′ = transition(f_logπq, f_gradlogπ, q₀, p₀, M, M⁻¹, ϵ, L)
        t += 1
        q′ != q₀ && return q′, p′, t
    end
    return q′, p′, t
end

############################################################################################
# Centered parameterization:
# yⱼ ~ N(αⱼ, σⱼ²)
# αⱼ ~ N(μ, τ²)
# μ ~ U(-∞, ∞)
# τ ~ U(0, ∞)
# τ transformed to logτ to place on unconstrained space

function gradlogπ(α, μ, logτ, y, σ²)
    τ² = abs2(exp(logτ))
    # ∇l = @. - (α - y) / σ² - (α - μ) / τ²
    ∇l = .- (α .- y) ./ σ² .- (α .- μ) ./ τ²
    dπdμ = 0.
        for j ∈ axes(α, 1)
            dπdμ -= (μ - α[j]) / τ² # -
        end
    dπdlogτ = - length(α) + 1.
        for j ∈ axes(α, 1)
            dπdlogτ += abs2(α[j] - μ) / τ²
        end
    [∇l; dπdμ; dπdlogτ]
end
gradlogπ(q, y, σ²) = gradlogπ(q[1:8], q[9], q[10], y, σ²)

function logπq(α, μ, logτ, y, σ²)
    τ² = abs2(exp(logτ))
    # s = 0
    # s += logτ
    s = logτ
    for j ∈ axes(α, 1)
        s += -log(√(2π * σ²[j])) - 0.5abs2(y[j] - α[j]) / σ²[j] - 0.5abs2(α[j] - μ) / τ²
        # s += - 0.5log(2π * σ²[j]) - 0.5abs2(y[j] - α[j]) / σ²[j] - 0.5abs2(α[j] - μ) / τ²
    end
    s -= length(α) * log(√(2π * τ²))
    # s -= length(α) * 0.5log(2π * τ²)
    s
    # # unnormalized density
    # τ² = abs2(exp(logτ))
    # s = logτ
    # for j ∈ axes(α, 1)
    #     s += -0.5abs2(y[j] - α[j]) / σ²[j] - 0.5abs2(α[j] - μ) / τ²
    # end
    # s += -0.5log(2π * τ²)
    # s
end
logπq(q, y, σ²) = logπq(q[1:8], q[9], q[10], y, σ²)

# A faster non-allocating version to replace the simple convenience above:
function logϖq(q, y, σ²)
    @inbounds μ = q[9]
    @inbounds logτ = q[10]
    τ² = abs2(exp(logτ))
    s = logτ
    @inbounds @simd for j ∈ eachindex(y, σ²)
        # s += -0.5log(2π * σ²[j]) - 0.5abs2(y[j] - q[j]) / σ²[j] - 0.5abs2(q[j] - μ) / τ²
        s += -0.5(log(2π * σ²[j]) + abs2(y[j] - q[j]) / σ²[j] + abs2(q[j] - μ) / τ²)
    end
    s += -0.5 * length(y) * log(2π * τ²)
    s
end


############################################################################################
#### 2022-02-14: non-centered parameterization
# yⱼ ~ N(μ + τ * ηⱼ, σⱼ²)
# ηⱼ ~ N(0, 1)
# μ ~ U(-∞, ∞)
# τ ~ U(0, ∞)
# τ transformed to logτ to place on unconstrained space

# Clearly not the most efficient approach -- the proper way is to compute a vector, the
# elements of which are (y[j] - μ - τ * η[j]). Then, compute dπdμ and dπdlogτ, then
# mutate in place to complete the computation of ∇l.
function gradlogπ_nc(η, μ, logτ, y, σ²)
    τ = exp(logτ)
    ∇l = Vector{Float64}(undef, length(η))
    for j ∈ eachindex(η, y, σ²)
        ∇l[j] = (y[j] - μ - τ * η[j]) * τ / σ²[j] - η[j]
    end
    dπdμ = 0.0
    for j ∈ eachindex(η, y, σ²)
        dπdμ += (y[j] - μ - τ * η[j]) / σ²[j]
    end
    dπdlogτ = 1.0
    for j ∈ eachindex(η, y, σ²)
        dπdlogτ += τ * (y[j] - μ - τ * η[j]) * η[j] / σ²[j]
    end
    [∇l; dπdμ; dπdlogτ]
end
gradlogπ_nc(q, y, σ²) = gradlogπ_nc(q[1:8], q[9], q[10], y, σ²)

function logπq_nc(η, μ, logτ, y, σ²)
    τ = exp(logτ)
    s = logτ
    for j ∈ eachindex(η, y, σ²)
        s += -0.5(log(2π * σ²[j]) + (1.0 / σ²[j]) * abs2(y[j] - μ - τ * η[j]) + log(2π) + abs2(η[j]))
    end
    s
end
logπq_nc(q, y, σ²) = logπq_nc(q[1:8], q[9], q[10], y, σ²)

# A faster non-allocating version to replace the simple convenience above:
function gradlogϖ_nc(q::Vector{T}, y, σ²) where {T<:Real}
    @inbounds μ = q[9]
    @inbounds logτ = q[10]
    τ = exp(logτ)
    ∇l = Vector{T}(undef, length(q))
    # for j ∈ eachindex(y, σ²)
    #     ∇l[j] = (y[j] - μ - τ * q[j]) * τ / σ²[j] - q[j]
    # end
    # dπdμ = 0.0
    # for j ∈ eachindex(y, σ²)
    #     dπdμ += (y[j] - μ - τ * q[j]) / σ²[j]
    # end
    # ∇l[9] = dπdμ
    # dπdlogτ = 1.0
    # for j ∈ eachindex(y, σ²)
    #     dπdlogτ += τ * (y[j] - μ - τ * q[j]) * q[j] / σ²[j]
    # end
    # ∇l[10] = dπdlogτ
    # ∇l
    # or, even better
    dπdμ = zero(T)
    dπdlogτ = one(T)
    @inbounds @simd for j ∈ eachindex(y, σ²)
        # yⱼ = y[j]
        # ηⱼ = q[j]
        # σⱼ² = σ²[j]
        # Δ = (yⱼ - μ - τ * ηⱼ)
        # ∇l[j] = Δ * τ / σⱼ² - ηⱼ
        # dπdμ += Δ / σⱼ²
        # dπdlogτ += τ * Δ * ηⱼ / σⱼ²
        yⱼ = y[j]
        ηⱼ = q[j]
        iσⱼ² = inv(σ²[j])
        Δ = yⱼ - μ - τ * ηⱼ
        Δiσⱼ² = Δ * iσⱼ²
        ∇l[j] = τ * Δiσⱼ² - ηⱼ
        dπdμ += Δiσⱼ²
        dπdlogτ += τ * Δiσⱼ² * ηⱼ
    end
    ∇l[9] = dπdμ
    ∇l[10] = dπdlogτ
    ∇l
end

function logϖq_nc(q, y, σ²)
    @inbounds μ = q[9]
    @inbounds logτ = q[10]
    τ = exp(logτ)
    s = logτ
    # @turbo for j ∈ eachindex(y, σ²)
    #     s += -0.5(log(2π * σ²[j]) + (1.0 / σ²[j]) * abs2(y[j] - μ - τ * q[j]) + log(2π) + abs2(q[j]))
    # end
    @inbounds @simd for j ∈ eachindex(y, σ²)
        ηⱼ = q[j]
        σⱼ² = σ²[j]
        s += -0.5(log(2π * σⱼ²) + abs2(y[j] - μ - τ * ηⱼ) / σⱼ² + log(2π) + abs2(ηⱼ))
    end
    s
end

#### Utilities for conversion to original parameterization
function originalparam(θ::Vector{T}) where {T<:Real}
    K = length(θ)
    τ = exp(θ[K])
    μ = θ[K - 1]
    ϑ = Vector{promote_type(T, Float64)}(undef, K)
    for k = 1:(K - 2)
        ϑ[k] = μ + τ * θ[k]
    end
    ϑ[K - 1] = μ
    ϑ[K] = τ
    ϑ
end

function originalparam(θ::Matrix{T}) where {T<:Real}
    K, S = size(θ)
    ϑ = Matrix{promote_type(T, Float64)}(undef, K, S)
    for s ∈ axes(θ, 2)
        τ = exp(θ[K, s])
        μ = θ[K - 1, s]
        for k = 1:(K - 2)
            ϑ[k, s] = μ + τ * θ[k, s]
        end
        ϑ[K - 1, s] = μ
        ϑ[K, s] = τ
    end
    ϑ
end
originalparam(θ::Array{T, 3}) where {T<:Real} = mapslices(originalparam, θ, dims=(1, 2))

############################################################################################
#### MCMC diagnostics: potential scale reduction and effective sample size
# Simulations arranged on dimension 2, chains on dimension 3
# function withinchainmean(chain::AbstractMatrix)
#     ψ̄ = mean(chain, dims=2)
# end
function betweenchain(ϑ::Array{T, 3}) where {T<:Real}
    B = size(ϑ, 2) .* var(dropdims(mean(ϑ, dims=2), dims=2), dims=2)
end

function withinchain(ϑ::Array{T, 3}) where {T<:Real}
    W = mean(dropdims(var(ϑ, dims=2), dims=2), dims=2)
end

function vhat(ψ::AbstractMatrix{T}) where {T<:Real}
    n = size(ψ, 1)
    B = n * var(mean(ψ, dims=1))
    W = mean(var(ψ, dims=1))
    v̂⁺ = ((n - 1) / n) * W + (1 / n) * B
end


function vhat(ϑ::Array{T, 3}) where {T<:Real}
    n = size(ϑ, 2)
    B = n .* var(dropdims(mean(ϑ, dims=2), dims=2), dims=2)
    W = mean(dropdims(var(ϑ, dims=2), dims=2), dims=2)
    v̂⁺ = ((n - 1) / n) .* W .+ (1 / n) .* B
end

function Rhat(ψ::AbstractMatrix{T}) where {T<:Real}
    n = size(ψ, 1)
    B = n * var(mean(ψ, dims=1))
    W = mean(var(ψ, dims=1))
    v̂⁺ = ((n - 1) / n) * W + (1 / n) * B
    R̂ = √(v̂⁺ / W)
end

function Rhat(ϑ::Array{T, 3}) where {T<:Real}
    n = size(ϑ, 2)
    B = n .* var(dropdims(mean(ϑ, dims=2), dims=2), dims=2)
    W = mean(dropdims(var(ϑ, dims=2), dims=2), dims=2)
    v̂⁺ = ((n - 1) / n) .* W .+ (1 / n) .* B
    R̂ = sqrt.(v̂⁺ ./ W)
end

function variogram(ψ::AbstractMatrix{T}, t::Int) where {T<:Real}
    n, m = size(ψ)
    s = zero(T)
    for j ∈ axes(ψ, 2)
        for i = (t + 1):n
            s += abs2(ψ[i, j] - ψ[i - t, j])
        end
    end
    return s / (m * (n - t))
end

autocor(ψ::AbstractMatrix{T}, t::Int) where {T<:Real} = one(T) - variogram(ψ, t) / 2vhat(ψ)

function findautolag(ψ::AbstractMatrix{T}) where {T<:Real}
    n = size(ψ, 1)
    for t = 1:2:(n - 3)
        ρ̂ₜ₊₁ = autocor(ψ, t)
        ρ̂ₜ₊₂ = autocor(ψ, t + 1)
        ρ̂ₜ₊₁ + ρ̂ₜ₊₂ < 0 && return t
    end
    return n - 3
end

function neff(ψ::AbstractMatrix{T}) where {T<:Real}
    n, m = size(ψ)
    # t = findautolag(ψ)
    Σρ = zero(T)
    for t = 1:2:(n - 3)
        ρ̂ₜ₊₁ = autocor(ψ, t)
        ρ̂ₜ₊₂ = autocor(ψ, t + 1)
        s = ρ̂ₜ₊₁ + ρ̂ₜ₊₂
        s < 0 && break
        Σρ += s
        # s < 0 ? break : (Σρ += s)
    end
    n̂ₑ = n * m / (one(T) + 2Σρ)
end

neff(ϑ::Array{T, 3}) where {T<:Real} = mapslices(neff, ϑ, dims=(2, 3))

############################################################################################
using Plots
gr(size=(1200,800))
# Example for non-centered parameterization
# Data
y = Float64[28, 8, -3, 7, -1, 1, 18, 12];
σ = Float64[15, 10, 16, 11, 9, 11, 10, 18];
σ² = abs2.(σ);

# Integration time and stepsize. These are suggested values -- try some others, keeping ϵL=1
ϵ = 0.2
L = 5
N = 1000
N_chains = 8

# Euclidean-Gaussian metric
scale = 15;
M⁻¹ = diagm(ones(10)); M⁻¹[9, 9] = scale^2;
M = diagm(1 ./ diag(M⁻¹));

f_logπ = let y = y, σ² = σ²
    q -> logπq_nc(q, y, σ²)
end;
f_gradlogπ = let y = y, σ² = σ²
    q -> gradlogπ_nc(q, y, σ²)
end;
f_logϖ = let y = y, σ² = σ²
    q -> logϖq_nc(q, y, σ²)
end;
f_gradlogϖ = let y = y, σ² = σ²
    q -> gradlogϖ_nc(q, y, σ²)
end;

@timev θ, E, n, t = hmc(f_logπ, f_gradlogπ, M, M⁻¹, ϵ, L, N);
@timev θ, E, n, t = hmc(f_logϖ, f_gradlogϖ, M, M⁻¹, ϵ, L, N);
ϑ = originalparam(θ);
mean(ϑ, dims=2)

ebfmi(E)

Epq, Eq = microcanonicalenergy(f_logπ, θ, E);
pₕ = histogram([Epq Eq] .- [mean(Epq) mean(Eq)], labels=["E_pq ≡ π(E|q)" "E_q ≡ π(E)"], alpha=0.5,
               annotations=((0.1, 1.0), text("E-BFMI = $(round(ebfmi(E), digits=3))", 8)),
               xlabel="E - <E>");
savefig(pₕ, joinpath(pwd(), "energy.pdf"))

# Multiple chains
@timev chains = [hmc(f_logπ, f_gradlogπ, M, M⁻¹, ϵ, L, N) for _ = 1:N_chains];
acprob = getindex.(chains, 3) ./ getindex.(chains, 4)
θ = cat(first.(chains)..., dims=3);
mean(θ, dims=(2,3))
E = hcat(getindex.(chains, 2)...);

ϑ₀ = originalparam(θ);
ϑ = originalparam(θ[:, (N >> 1 + 1):N, :]);
mean(ϑ, dims=(2, 3))
# Convergence check
n₂ = size(ϑ, 2)
n = n₂ >> 1
ϑ₂ = [view(ϑ, :, 1:n, :);;; view(ϑ, :, (n + 1):n₂, :)];
neff(ϑ₂)
Rhat(ϑ₂)

ps = map(1:N_chains) do k
    Epq, Eq = microcanonicalenergy(f_logπ, θ[:, :, k], E[:, k])
    pₕ = histogram([Epq Eq] .- [mean(Epq) mean(Eq)], labels=["E_pq ≡ π(E|q)" "E_q ≡ π(E)"],
                   alpha=0.5,
                   annotations=((0.1, 1.0), text("E-BFMI = $(round(ebfmi(E[:, k]), digits=3))", 8)),
                   xlabel="E - <E>")
    pₕ
end;
pₕ = plot(ps...);
savefig(pₕ, joinpath(pwd(), "energy_$(N_chains)chain.pdf"))

# Neat trace plots
ps = map(1:N_chains) do k
    plot(permutedims(ϑ₀[:, :, k]), labels=["θ$(tosubscript(i))" for _ = 1:1, i = 1:10])
end;
p = plot(ps..., link=:all);
savefig(p, joinpath(pwd(), "hmc_$(N_chains)chains_trace.pdf"))

# autocorrelation plots
ϑⱼ = ϑ[1, :];
acf = map(0:500) do k
    acor(ϑⱼ, k)
end;
p = bar(acf, label="empirical")
savefig(p, joinpath(pwd(), "hmc_autocorrelation.pdf"))
#
ps = map(1:size(ϑ, 1)) do j
    ϑⱼ = ϑ[j, :]
    acf = map(0:500) do k
        acor(ϑⱼ, k)
    end;
    p = bar(acf, label="θ$(j)")
    p
end;
p = plot(ps...);
savefig(p, joinpath(pwd(), "hmc_autocorrelation.pdf"))

#################### Centered parameterization -- might get stuck due to small τ
ϵ = 0.1
L = 10
N = 1000
N_chains = 8

scale = 15;
M = diagm(fill(1 / scale^2, 10)); M[10, 10] = 1;
M⁻¹ = diagm(1 ./ diag(M));

f_logπ_c = let y = y, σ² = σ²
    q -> logπq(q, y, σ²)
end;
f_gradlogπ_c = let y = y, σ² = σ²
    q -> gradlogπ(q, y, σ²)
end;


# Typically, centered parameterization necessitates a non-small τ starting point, so guarantee it
qᵢ, pᵢ, tᵢ = initialize(f_logπ_c, f_gradlogπ_c, M, M⁻¹, ϵ, L)
@timev θ, E, n, t = hmc(f_logπ_c, f_gradlogπ_c, qᵢ, M, M⁻¹, ϵ, L, N);

ebfmi(E)

Epq, Eq = microcanonicalenergy(f_logπ_c, θ, E);
pₕ = histogram([Epq Eq] .- [mean(Epq) mean(Eq)], labels=["E_pq ≡ π(E|q)" "E_q ≡ π(E)"], alpha=0.5,
               annotations=((0.1, 1.0), text("E-BFMI = $(round(ebfmi(E), digits=3))", 8)),
               xlabel="E - <E>");
savefig(pₕ, joinpath(pwd(), "energy_c.pdf"))

# Multiple chains
@timev chains = [hmc(f_logπ_c, f_gradlogπ_c,
                     first(initialize(f_logπ_c, f_gradlogπ_c, M, M⁻¹, ϵ, L)), # qᵢ's
                     M, M⁻¹, ϵ, L, N) for _ = 1:N_chains];
acprob = getindex.(chains, 3) ./ getindex.(chains, 4)
θ = cat(first.(chains)..., dims=3);
mean(θ, dims=(2,3))
E = hcat(getindex.(chains, 2)...);

# Convergence check
neff(θ)
Rhat(θ)

ϑ = θ[:, (N >> 1 + 1):N, :];
mean(ϑ, dims=(2, 3))
# Convergence check
n₂ = size(ϑ, 2)
n = n₂ >> 1
ϑ₂ = [view(ϑ, :, 1:n, :);;; view(ϑ, :, (n + 1):n₂, :)];
neff(ϑ₂)
Rhat(ϑ₂)

ps = map(1:N_chains) do k
    Epq, Eq = microcanonicalenergy(f_logπ_c, θ[:, :, k], E[:, k])
    pₕ = histogram([Epq Eq] .- [mean(Epq) mean(Eq)], labels=["E_pq ≡ π(E|q)" "E_q ≡ π(E)"],
                   alpha=0.5,
                   annotations=((0.1, 1.0), text("E-BFMI = $(round(ebfmi(E[:, k]), digits=3))", 8)),
                   xlabel="E - <E>")
    pₕ
end;
pₕ = plot(ps...);
savefig(pₕ, joinpath(pwd(), "energy_c_$(N_chains)chain.pdf"))


# Neat trace plots
ps = map(1:N_chains) do k
    plot(permutedims(θ[:, :, k]), labels=["θ$(tosubscript(i))" for _ = 1:1, i = 1:10])
end;
p = plot(ps..., link=:all);
savefig(p, joinpath(pwd(), "hmc_c_$(N_chains)chains_trace.pdf"))

# autocorrelation plots
θⱼ = θ[1, :];
acf = map(0:500) do k
    acor(θⱼ, k)
end;
p = bar(acf, label="empirical")
savefig(p, joinpath(pwd(), "hmc_c_autocorrelation.pdf"))
#
ps = map(1:size(θ, 1)) do j
    θⱼ = θ[j, :]
    acf = map(0:500) do k
        acor(θⱼ, k)
    end;
    p = bar(acf, label="θ$(j)")
    p
end;
p = plot(ps...);
savefig(p, joinpath(pwd(), "hmc_c_autocorrelation.pdf"))

############################################################################################
#### 2022-02-16: Gibbs sampler, centered parameterization
using Distributions

function sample_α(μ, τ², y, σ²)
    α = randn(length(y))
    iτ² = inv(τ²)
    μiτ² = μ / τ²
    for j ∈ eachindex(y, σ²)
        iσⱼ² = inv(σ²[j])
        d = iσⱼ² + iτ²
        α̂ⱼ = (y[j] * iσⱼ² + μiτ²) / d
        Vⱼ = inv(d)
        α[j] = α̂ⱼ + α[j] * √Vⱼ
    end
    α
end

sample_μ(α, τ²) = mean(α) + randn() * √(τ² / length(α))

function sample_τ²(α, μ)
    τ̂² = 0.0
    for j ∈ eachindex(α)
        τ̂² += abs2(α[j] - μ)
    end
    ν = length(α) - 1
    # # Unnecessary since it is immediately canceled, but worthwhile to show:
    # τ̂² /= ν # compltes τ̂² computation
    # τ² = ν * τ̂² / rand(Chisq(ν))
    τ̂² / rand(Chisq(ν)) # τ̂² / rand(Gamma(ν / 2, 2.0)) is slightly faster
end

function gibbs(α₀, μ₀, y, σ², N)
    θ = Matrix{Float64}(undef, 10, N)
    α = α₀
    μ = μ₀
    for n = 1:N
        τ² = sample_τ²(α, μ)
        μ = sample_μ(α, τ²)
        α = sample_α(μ, τ², y, σ²)
        θ[1:8, n] = α
        θ[9, n] = μ
        θ[10, n] = τ²
    end
    θ
end

# Experiment with different initial points
α₀ = randn(8) * 15; # α₀ ~ N(0, 15²I₈)
# α₀ = randn(8);
# α₀ = y;
μ₀ = mean(α₀) + randn();

N = 1000
N_chains = 4

@timev θ = gibbs(α₀, μ₀, y, σ², N);
mean(θ[:, (N >> 1 + 1):N], dims=2)
mean(sqrt.(θ[10, :]))

p = plot(permutedims(θ), labels=["θ$(tosubscript(i))" for _ = 1:1, i = 1:10]);
savefig(p, joinpath(pwd(), "trace.pdf"))
p = plot(θ[9, :], θ[10, :], 1:N, zcolor=(1:N), m = (10, 0.8, :blues, Plots.stroke(0)),
         leg = false, cbar = true);
savefig(p, joinpath(pwd(), "3d.pdf"))

# Multiple chains
@timev chains = [gibbs(α₀, μ₀, y, σ², N) for _ = 1:N_chains];
θ = cat(chains..., dims=3);
mean(θ, dims=(2,3))
mean(sqrt.(θ[10, :, :]))

# Convergence check
neff(θ)
Rhat(θ)

ϑ = θ[:, (N >> 1 + 1):N, :];
mean(ϑ, dims=(2, 3))
# Convergence check
n₂ = size(ϑ, 2)
n = n₂ >> 1
ϑ₂ = [view(ϑ, :, 1:n, :);;; view(ϑ, :, (n + 1):n₂, :)];
neff(ϑ₂)
Rhat(ϑ₂)

# Neat trace plots
ps = map(1:N_chains) do k
    plot(permutedims(θ[:, :, k]), labels=["θ$(tosubscript(i))" for _ = 1:1, i = 1:10])
end;
p = plot(ps..., link=:all);
savefig(p, joinpath(pwd(), "gibbs_$(N_chains)chains_trace.pdf"))

# Complimentary to the above, this plot enables one to visualize when the Gibbs sampler
# (under the centered parameterization) get stuck at small values of τ. The
# effect is quite dramatic!
ps = map(1:N_chains) do k
    θₖ = θ[:, :, k]
    θₖ[10, :] .= sqrt.(θₖ[10, :])
    plot(permutedims(θₖ), labels=["θ$(tosubscript(i))" for _ = 1:1, i = 1:10])
end;
p = plot(ps..., link=:all);
savefig(p, joinpath(pwd(), "gibbs_$(N_chains)chains_trace_sqrttau.pdf"))

# autocorrelation plots
θⱼ = θ[1, :];
acf = map(0:500) do k
    acor(θⱼ, k)
end;
p = bar(acf, label="empirical")
savefig(p, joinpath(pwd(), "gibbs_autocorrelation.pdf"))
#
ps = map(1:size(θ, 1)) do j
    θⱼ = θ[j, :]
    acf = map(0:500) do k
        acor(θⱼ, k)
    end;
    p = bar(acf, label="θ$(j)")
    p
end;
p = plot(ps...);
savefig(p, joinpath(pwd(), "gibbs_autocorrelation.pdf"))

############################################################################################
#### 2022-02-16: Metropolis-Hastings, centered parameterization, d-dimensional vector-jumping

function randmvndiag(μ, Σ)
    z = randn(length(μ))
    for j ∈ eachindex(μ, z)
        z[j] = μ[j] + z[j] * √(Σ[j, j])
    end
    z
end
function randmvndiag(Σ)
    z = randn(size(Σ, 1))
    @inbounds for j ∈ axes(Σ, 1)
        z[j] *= √(Σ[j, j])
    end
    z
end

function logp_mvndiag(x, μ, Σ)
    -0.5(length(x) * log(2π) + log(diagdet(Σ)) + diagweightedss(x, μ, Σ))
end

function metropolishastings(f_logπq, f_logQ, Σ, N)
    θ = Matrix{Float64}(undef, size(Σ, 1), N)
    n, t = 1, 0
    q₀ = randmvndiag(Σ)
    θ[:, n] = q₀
    while n < N
        q′ = randmvndiag(q₀, Σ)
        r = exp(f_logπq(q′) + f_logQ(q₀, q′, Σ) - f_logπq(q₀) - f_logQ(q′, q₀, Σ))
        if rand() < r
            n += 1
            q₀ = q′
            θ[:, n] = q′
        end
        t += 1
    end
    θ, n, t
end

f_logπq = let y = y, σ² = σ²
    q -> logπq(q, y, σ²)
end;

N = 5000
N_chains = 4

# It is instructive to set the scale to values 0.5, 1.0, 5.0, 15.0, etc.
# and check convergence.
scale = .5;
c = 2.4 / √(10)
Σ = c^2 * diagm(fill(scale^2, 10)); Σ[10, 10] = 1;

@timev θ, n, t = metropolishastings(f_logπq, logp_mvndiag, Σ, N);
mean(θ[:, (N >> 1 + 1):N], dims=2)

p = plot(permutedims(θ), labels=["θ$(tosubscript(i))" for _ = 1:1, i = 1:10]);
savefig(p, joinpath(pwd(), "trace.pdf"))
p = plot(θ[9, :], θ[10, :], 1:N, zcolor=(1:N), m = (10, 0.8, :blues, Plots.stroke(0)),
         leg = false, cbar = true);
savefig(p, joinpath(pwd(), "3d.pdf"))

# Multiple chains
@timev chains = [metropolishastings(f_logπq, logp_mvndiag, Σ, N) for _ = 1:N_chains];
θ = cat(first.(chains)..., dims=3);
acprob = getindex.(chains, 2) ./ getindex.(chains, 3)
mean(θ, dims=(2,3))
mean(exp.(θ[10, :, :]))

# Convergence check
neff(θ)
Rhat(θ)

ϑ = θ[:, (N >> 1 + 1):N, :];
mean(ϑ, dims=(2, 3))
# Convergence check
n₂ = size(ϑ, 2)
n = n₂ >> 1
ϑ₂ = [view(ϑ, :, 1:n, :);;; view(ϑ, :, (n + 1):n₂, :)];
neff(ϑ₂)
Rhat(ϑ₂)

# Neat trace plots
ps = map(1:N_chains) do k
    plot(permutedims(θ[:, :, k]), labels=["θ$(tosubscript(i))" for _ = 1:1, i = 1:10])
end;
p = plot(ps..., link=:all);
savefig(p, joinpath(pwd(), "mh_$(N_chains)chains_trace.pdf"))

# autocorrelation plots
θⱼ = θ[1, :];
acf = map(0:500) do k
    acor(θⱼ, k)
end;
p = bar(acf, label="empirical")
savefig(p, joinpath(pwd(), "mh_autocorrelation.pdf"))
#
ps = map(1:size(θ, 1)) do j
    θⱼ = θ[j, :]
    acf = map(0:500) do k
        acor(θⱼ, k)
    end;
    p = bar(acf, label="θ$(j)")
    p
end;
p = plot(ps...);
savefig(p, joinpath(pwd(), "mh_autocorrelation.pdf"))

############################################################################################
#### 2022-02-17: Efficient static HMC
function etransition(f_logπq, f_gradlogπ, q₀, p₀, M, M⁻¹, w₀, ϵ, T, L)
    if ϵ * L ≥ T
        return q₀, p₀
    else
        ϵ = rand(Bool) ? -ϵ : ϵ
        w′ = 0.0
        p′ = p₀ + 0.5ϵ * f_gradlogπ(q₀)
        q′ = q₀ + ϵ * M⁻¹ * p′
        for l = 2:L
            p′ = p′ + ϵ * f_gradlogπ(q′)
            w′ += exp(-H(f_logπq, q′, p′, M⁻¹))
            q′ = q′ + ϵ * M⁻¹ * p′
        end
        p′ = p′ + 0.5ϵ * f_gradlogπ(q′)
        w′ += exp(-H(f_logπq, q′, p′, M⁻¹))
        𝓅₀ = w₀ / (w₀ + w′)
        if rand() ≤ 𝓅₀
            # return etransition(f_logπq, f_gradlogπ, q₀, rand(Bool) ? -p₀ : p₀,
            #                    M, M⁻¹, w₀ + w′, ϵ, T, 2L)
            return etransition(f_logπq, f_gradlogπ, q₀, p₀, M, M⁻¹, w₀ + w′, abs(ϵ), T, 2L)
            # return etransition(f_logπq, f_gradlogπ, q₀, p₀, M, M⁻¹, w₀ + w′, T, 2L
            #                    rand(Bool) ? -ϵ : ϵ, T, 2L)
        else
            # return etransition(f_logπq, f_gradlogπ, q′, rand(Bool) ? -p′ : p′,
            #                    M, M⁻¹, w₀ + w′, ϵ, T, 2L)
            return etransition(f_logπq, f_gradlogπ, q′, p′, M, M⁻¹, w₀ + w′, abs(ϵ), T, 2L)
            # return etransition(f_logπq, f_gradlogπ, q′, p′, M, M⁻¹, w₀ + w′,
            #                    rand(Bool) ? -ϵ : ϵ, T, 2L)
        end
    end
end
function etransition(f_logπq, f_gradlogπ, q₀, M, M⁻¹, ϵ, T)
    p₀ = momentumrand_M(M)
    w₀ = exp(-H(f_logπq, q₀, p₀, M⁻¹))
    etransition(f_logπq, f_gradlogπ, q₀, p₀, M, M⁻¹, w₀, ϵ, T, 1)
end
function etransition(f_logπq, f_gradlogπ, M, M⁻¹, ϵ, T)
    q₀ = targetrand(M⁻¹)
    etransition(f_logπq, f_gradlogπ, q₀, M, M⁻¹, ϵ, T)
end


# Version which explicitly constructs the trajectories, then samples them
function etransition2!(f_logπq, f_gradlogπ, t₀, t′, q₀, p₀, M, M⁻¹, w₀, ϵ, T, L)
    if ϵ * L ≥ T
        return q₀, p₀
    else
        ϵ = rand(Bool) ? -ϵ : ϵ
        resize!(t′, L)
        w′ = 0.0
        p′ = p₀ + 0.5ϵ * f_gradlogπ(q₀)
        q′ = q₀ + ϵ * M⁻¹ * p′
        for l = 2:L
            p′ = p′ + ϵ * f_gradlogπ(q′)
            t′[l - 1] = (q′, p′)
            w′ += exp(-H(f_logπq, q′, p′, M⁻¹))
            q′ = q′ + ϵ * M⁻¹ * p′
        end
        p′ = p′ + 0.5ϵ * f_gradlogπ(q′)
        t′[L] = (q′, p′)
        w′ += exp(-H(f_logπq, q′, p′, M⁻¹))
        𝓅₀ = w₀ / (w₀ + w′)
        if rand() ≤ 𝓅₀
            q, p = rand(t₀)
            signbit(ϵ) ? prepend!(t₀, t′) : append!(t₀, t′)
            return etransition2!(f_logπq, f_gradlogπ, t₀, t′, q, p, M, M⁻¹, w₀ + w′, abs(ϵ), T, 2L)
        else
            q, p = rand(t′)
            signbit(ϵ) ? prepend!(t₀, t′) : append!(t₀, t′)
            return etransition2!(f_logπq, f_gradlogπ, t₀, t′, q, p, M, M⁻¹, w₀ + w′, abs(ϵ), T, 2L)
        end
    end
end
function etransition2(f_logπq, f_gradlogπ, q₀, M, M⁻¹, ϵ, T)
    p₀ = momentumrand_M(M)
    w₀ = exp(-H(f_logπq, q₀, p₀, M⁻¹))
    t₀ = Vector{Tuple{Vector{Float64}, Vector{Float64}}}(undef, 0)
    t′ = Vector{Tuple{Vector{Float64}, Vector{Float64}}}(undef, 0)
    push!(t₀, (q₀, p₀))
    etransition2!(f_logπq, f_gradlogπ, t₀, t′, q₀, p₀, M, M⁻¹, w₀, ϵ, T, 1)
end
function etransition2(f_logπq, f_gradlogπ, M, M⁻¹, ϵ, T)
    q₀ = targetrand(M⁻¹)
    etransition2(f_logπq, f_gradlogπ, q₀, M, M⁻¹, ϵ, T)
end


function ehmc(f_logπq, f_gradlogπ, q₀, M, M⁻¹, ϵ, T, N)
    θ = Matrix{Float64}(undef, size(M, 1), N)
    E = Vector{Float64}(undef, N)
    n, t = 0, 0
    while n < N
        q′, p′ = etransition(f_logπ, f_gradlogπ, q₀, M, M⁻¹, ϵ, T)
        if q₀ != q′
            q₀ = q′
            n += 1
            θ[:, n] = q′
            E[n] = H(f_logπq, q′, p′, M⁻¹)
        end
        t += 1
    end
    θ, E, n, t
end
ehmc(f_logπq, f_gradlogπ, M, M⁻¹, ϵ, T, N) =
    ehmc(f_logπq, f_gradlogπ, targetrand(M⁻¹), M, M⁻¹, ϵ, T, N)

ϵ = 0.125
L = 8
N = 1000
N_chains = 4

@timev etransition(f_logπ, f_gradlogπ, M, M⁻¹, ϵ, ϵ * L)
@timev etransition2(f_logπ, f_gradlogπ, M, M⁻¹, ϵ, ϵ * L)

@timev θ, E, n, t = ehmc(f_logπ, f_gradlogπ, M, M⁻¹, ϵ, ϵ * L, N)
ϑ = originalparam(θ);
mean(ϑ, dims=2)

ebfmi(E)

Epq, Eq = microcanonicalenergy(f_logπ, θ, E);
pₕ = histogram([Epq Eq] .- [mean(Epq) mean(Eq)], labels=["E_pq ≡ π(E|q)" "E_q ≡ π(E)"], alpha=0.5,
               annotations=((0.1, 1.0), text("E-BFMI = $(round(ebfmi(E), digits=3))", 8)),
               xlabel="E - <E>");
savefig(pₕ, joinpath(pwd(), "energy_ehmc.pdf"))

@timev chains = [ehmc(f_logπ, f_gradlogπ, M, M⁻¹, ϵ, ϵ * L, N) for _ = 1:N_chains];
acprob = getindex.(chains, 3) ./ getindex.(chains, 4)
θ = cat(first.(chains)..., dims=3);
mean(θ, dims=(2,3))
E = hcat(getindex.(chains, 2)...);

ϑ₀ = originalparam(θ);
ϑ = originalparam(θ[:, (N >> 1 + 1):N, :]);
mean(ϑ, dims=(2, 3))
# Convergence check
n₂ = size(ϑ, 2)
n = n₂ >> 1
ϑ₂ = [view(ϑ, :, 1:n, :);;; view(ϑ, :, (n + 1):n₂, :)];
neff(ϑ₂)
Rhat(ϑ₂)

############################################################################################
#### 2022-02-17: dynamic HMC
function dotdiff(p′::Vector{T}, q′::Vector{T}, q₀::Vector{T}) where {T<:Real}
    s = zero(T)
    @inbounds @simd for i ∈ eachindex(p′, q′, q₀)
        Δq = q′[i] - q₀[i]
        s += p′[i] * Δq
    end
    s
end
function nutscriterion(q′, p′, q₀, p₀)
    # Δq = q′ - q₀
    # p′ ⋅ Δq < 0.0 && p₀ ⋅ Δq > 0.0
    # p′ ⋅ (q′ - q₀) < 0.0 && p₀ ⋅ (q₀ - q′) < 0.0
    dotdiff(p′, q′, q₀) < 0.0 && dotdiff(p₀, q₀, q′) < 0.0
end

function momentumstep!(p′::Vector{T}, ϵ::T, gradlogπq::Vector{T}) where {T<:Real}
    @inbounds @simd for i ∈ eachindex(p′, gradlogπq)
        p′[i] += ϵ * gradlogπq[i]
    end
    p′
end

function targetstepdiag!(q′::Vector{T}, ϵ::T, M⁻¹::Matrix{T}, p′::Vector{T}) where {T<:Real}
    @inbounds @simd for i ∈ eachindex(q′, p′)
        q′[i] += ϵ * M⁻¹[i, i] * p′[i]
    end
    q′
end

function dtransition(f_logπq, f_gradlogπ, q₀, p₀, M, M⁻¹, w₀, ϵ, L, Lₘ)
    L ≥ Lₘ && return q₀, p₀, L
    ϵ = rand(Bool) ? -ϵ : ϵ
    w′ = 0.0
    p′ = p₀ + 0.5ϵ * f_gradlogπ(q₀)
    q′ = q₀ + ϵ * M⁻¹ * p′
    #### Version for clarity
    # for l = 2:L
    #     p′ = p′ + ϵ * f_gradlogπ(q′)
    #     # nutscriterion(q′, p′, q₀, p₀) && return dtransition(f_logπq, f_gradlogπ, q₀, p₀, M, M⁻¹, w₀, ϵ, 2L)
    #     w′ += exp(-H(f_logπq, q′, p′, M⁻¹))
    #     # nutscriterion(q′, p′, q₀, p₀) && return rand() ≤ w₀ / (w₀ + w′) ? (q₀, p₀) : (q′, p′)
    #     q′ = q′ + ϵ * M⁻¹ * p′
    # end
    #### mutating version
    for l = 2:L
        momentumstep!(p′, ϵ, f_gradlogπ(q′))
        w′ += exp(-H(f_logπq, q′, p′, M⁻¹))
        targetstepdiag!(q′, ϵ, M⁻¹, p′)
    end
    ####
    # p′ = p′ + 0.5ϵ * f_gradlogπ(q′)
    momentumstep!(p′, 0.5ϵ, f_gradlogπ(q′))
    w′ += exp(-H(f_logπq, q′, p′, M⁻¹))
    # 𝓅₀ = w₀ / (w₀ + w′)
    # turned = nutscriterion(q′, p′, q₀, p₀)
    # if rand() ≤ 𝓅₀
    #     # return dtransition(f_logπq, f_gradlogπ, q₀, rand(Bool) ? -p₀ : p₀,
    #     #                    M, M⁻¹, w₀ + w′, ϵ, T, 2L)
    #     return turned ? (q₀, p₀) : dtransition(f_logπq, f_gradlogπ, q₀, p₀, M, M⁻¹, w₀ + w′, ϵ, 2L)
    # else
    #     # return dtransition(f_logπq, f_gradlogπ, q′, rand(Bool) ? -p′ : p′,
    #     #                    M, M⁻¹, w₀ + w′, ϵ, T, 2L)
    #     return turned ? (q′, p′) : dtransition(f_logπq, f_gradlogπ, q′, p′, M, M⁻¹, w₀ + w′, ϵ, 2L)
    # end
    if nutscriterion(q′, p′, q₀, p₀)
        return q′, p′, L
    else
        𝓅₀ = w₀ / (w₀ + w′)
        return rand() ≤ 𝓅₀ ? dtransition(f_logπq, f_gradlogπ, q₀, p₀, M, M⁻¹, w₀ + w′, ϵ, 2L, Lₘ) :
            dtransition(f_logπq, f_gradlogπ, q′, p′, M, M⁻¹, w₀ + w′, ϵ, 2L, Lₘ)
    end
end
function dtransition(f_logπq, f_gradlogπ, q₀, M, M⁻¹, ϵ, Lₘ)
    p₀ = momentumrand_M(M)
    w₀ = exp(-H(f_logπq, q₀, p₀, M⁻¹))
    dtransition(f_logπq, f_gradlogπ, q₀, p₀, M, M⁻¹, w₀, ϵ, 1, Lₘ)
end
function dtransition(f_logπq, f_gradlogπ, M, M⁻¹, ϵ, Lₘ)
    q₀ = targetrand(M⁻¹)
    p₀ = momentumrand_M(M)
    w₀ = exp(-H(f_logπq, q₀, p₀, M⁻¹))
    dtransition(f_logπq, f_gradlogπ, q₀, p₀, M, M⁻¹, w₀, ϵ, 1, Lₘ)
end

function dhmc(f_logπq, f_gradlogπ, q₀, M, M⁻¹, ϵ, Lₘ, N)
    θ = Matrix{Float64}(undef, size(M, 1), N)
    E = Vector{Float64}(undef, N)
    depth = Vector{Int}(undef, N)
    n, t = 0, 0
    while n < N
        # q′, p′, L = dtransition(f_logπ, f_gradlogπ, q₀, M, M⁻¹, ϵ, Lₘ)
        # if q₀ != q′
        #     q₀ = q′
        #     n += 1
        #     θ[:, n] = q′
        #     E[n] = H(f_logπq, q′, p′, M⁻¹)
        #     depth[n] = L
        # end
        # t += 1
        p₀ = momentumrand_M(M)
        w₀ = exp(-H(f_logπq, q₀, p₀, M⁻¹))
        q′, p′, L = dtransition(f_logπ, f_gradlogπ, q₀, p₀, M, M⁻¹, w₀, ϵ, 1, Lₘ)
        E₀ = H(f_logπq, q₀, p₀, M⁻¹)
        E′ = H(f_logπq, q′, p′, M⁻¹)
        r = exp(-E₀ + E′)
        if rand() < r
            n += 1
            θ[:, n] = q′
            q₀ = q′
            E[n] = E′
            depth[n] = L
        end
        t += 1
    end
    θ, E, n, t, depth
end
dhmc(f_logπq, f_gradlogπ, M, M⁻¹, ϵ, Lₘ, N) =
    dhmc(f_logπq, f_gradlogπ, targetrand(M⁻¹), M, M⁻¹, ϵ, Lₘ, N)

ϵ = 0.01
Lₘ = 2^10
N = 1000
N_chains = 4

@timev dtransition(f_logπ, f_gradlogπ, M, M⁻¹, ϵ, Lₘ)

@timev θ, E, n, t, depth = dhmc(f_logπ, f_gradlogπ, M, M⁻¹, ϵ, Lₘ, N);
@timev θ, E, n, t, depth = dhmc(f_logϖ, f_gradlogϖ, M, M⁻¹, ϵ, Lₘ, N);
ϑ = originalparam(θ);
mean(ϑ, dims=2)
mean(log2, depth)
maximum(log2, depth)

ebfmi(E)

Epq, Eq = microcanonicalenergy(f_logπ, θ, E);
pₕ = histogram([Epq Eq] .- [mean(Epq) mean(Eq)], labels=["E_pq ≡ π(E|q)" "E_q ≡ π(E)"], alpha=0.5,
               annotations=((0.1, 1.0), text("E-BFMI = $(round(ebfmi(E), digits=3))", 8)),
               xlabel="E - <E>");
savefig(pₕ, joinpath(pwd(), "energy_dhmc.pdf"))

@timev chains = [dhmc(f_logπ, f_gradlogπ, M, M⁻¹, ϵ, Lₘ, N) for _ = 1:N_chains];
acprob = getindex.(chains, 3) ./ getindex.(chains, 4)
θ = cat(first.(chains)..., dims=3);
mean(θ, dims=(2,3))
E = hcat(getindex.(chains, 2)...);

ϑ₀ = originalparam(θ);
ϑ = originalparam(θ[:, (N >> 1 + 1):N, :]);
mean(ϑ, dims=(2, 3))
# Convergence check
n₂ = size(ϑ, 2)
n = n₂ >> 1
ϑ₂ = [view(ϑ, :, 1:n, :);;; view(ϑ, :, (n + 1):n₂, :)];
neff(ϑ₂)
Rhat(ϑ₂)

ps = map(1:N_chains) do k
    Epq, Eq = microcanonicalenergy(f_logπ, θ[:, :, k], E[:, k])
    pₕ = histogram([Epq Eq] .- [mean(Epq) mean(Eq)], labels=["E_pq ≡ π(E|q)" "E_q ≡ π(E)"],
                   alpha=0.5,
                   annotations=((0.1, 1.0), text("E-BFMI = $(round(ebfmi(E[:, k]), digits=3))", 8)),
                   xlabel="E - <E>")
    pₕ
end;
pₕ = plot(ps...);
savefig(pₕ, joinpath(pwd(), "dhmc_energy_$(N_chains)chain.pdf"))

# Neat trace plots
ps = map(1:N_chains) do k
    plot(permutedims(ϑ₀[:, :, k]), labels=["θ$(tosubscript(i))" for _ = 1:1, i = 1:10])
end;
p = plot(ps..., link=:all);
savefig(p, joinpath(pwd(), "dhmc_$(N_chains)chains_trace.pdf"))

# autocorrelation plots
ϑⱼ = ϑ[1, :];
acf = map(0:500) do k
    acor(ϑⱼ, k)
end;
p = bar(acf, label="empirical")
savefig(p, joinpath(pwd(), "dhmc_autocorrelation.pdf"))
#
ps = map(1:size(ϑ, 1)) do j
    ϑⱼ = ϑ[j, :]
    acf = map(0:500) do k
        acor(ϑⱼ, k)
    end;
    p = bar(acf, label="θ$(j)")
    p
end;
p = plot(ps...);
savefig(p, joinpath(pwd(), "dhmc_autocorrelation.pdf"))
