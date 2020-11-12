using Distributions, Roots
# function MvInv_slow(u, a, kappa, gama, M)
#     f(x) =  a/ Float64(real(gamma(GibbsTypePriors.CC(1- gama)))) * (x.^(-(1+gama)).*exp.( - (u+kappa)*x))

#     ξ = cumsum(rand(Exponential(1), M))
#     J = Array{Float64}(undef, M)
#     higherbound = 10
#     for i in eachindex(J)
#         obj(x) = f(x) - ξ[i]
#         J[i] = find_zero(obj, (0, higherbound))
#         startpos = J[i]
#     end
#     return J
# end
Γ(x) = SpecialFunctions.gamma(x)
using RCall
R"library(expint)"

function gamma_inc_r(a::Float64, x::Float64)
    res::Float64 = R"expint::gammainc($a,$x)"
    return res
end
# Γ(s,x) = SpecialFunctions.gamma_inc_r(s, x, 0)[2] * Γ(s)
Γ(s,x) = gamma_inc_r(s, x)

function MvInv_slow(u, a, κ, γ, M)
    N(x) =  a*κ^γ/SpecialFunctions.gamma(1-γ)*Γ(-γ,κ*x)
    ξ = cumsum(rand(Exponential(1), M))
    J = Array{Float64}(undef, M)
    higherbound = 10^3
    for i in eachindex(J)
        obj(x) = N(x) - ξ[i]
        # println((higherbound, obj.((8eps(), higherbound))))
        J[i] = find_zero(obj, (8eps(), higherbound))
        higherbound = J[i]
    end
    return J
end

function MvInv_slow_arb(u, a, κ, γ, M)
    N(x) =  a*κ^γ*Float64(real(1/gamma(CC(1-γ))*gamma(CC(-γ), CC(κ*x))))
    ξ = cumsum(rand(Exponential(1), M))
    J = Array{Float64}(undef, M)
    higherbound = 10
    for i in eachindex(J)
        obj(x) = N(x) - ξ[i]
        J[i] = find_zero(obj, (8eps(), higherbound))
        higherbound = J[i]
    end
    return J
end


function MvInv_simple(u, a, kappa, gama, N, M)
    J = Array{Float64}(undef, M)
    x = -log.(range(exp(-1e-05), exp(-10), length = N))
    f =  a/ Float64(real(gamma(CC(1- gama)))) * (x.^(-(1+gama)).*exp.( - (u+kappa)*x))
    dx = diff(x)
    h = (f[2:end] + f[1:(end-1)])/2
    Mv = zeros(N)
    ξ = cumsum(rand(Exponential(1), M))
    for i in (N-1):-1:1
        Mv[i] = Mv[i+1] + dx[i]*h[i]
    end
    for i in eachindex(J)
        J[i] = x[argmin(Mv .> ξ[i])]
    end
    return J
end

function fill_v2(M, Mv, W,N, x)
    v = Array{Float64}(undef, M)
    iMv = N
    for i in 1:M
        while iMv > 0 && Mv[iMv] < W[i]
            iMv = iMv -1
        end
        v[i] = x[iMv + 1]
    end
    return(v)
end

function MvInv_D(u, a, kappa, gama, N, M)
    x = -log.(range(exp(-1e-05),exp(-10), length=N))
    f =  a/ Float64(real(gamma(GibbsTypePriors.CC(1- gama)))) * (x.^(-(1+gama)).*exp.( - (u+kappa)*x))
    dx = diff(x)
    h = (f[2:end] + f[1:(end-1)]) / 2
    Mv = append!(reverse(cumsum(reverse(dx.* h))),[0])
    dist = Exponential(1)
    W = rand(dist, M)
    W = cumsum(W)
    return fill_v2(M,Mv,W,N,x)
end
