function Pkn_NGG_arb(k, n,  β, σ)
    if k>n || k==0
        return RR(0)
    else
        σ_arb = RR(σ)
        Vnk_NGG(n, k, β, σ) // σ_arb^k * Cnk(n, k, σ)
    end
end

function Pkn_NGG_robust(k, n,  β, σ; verbose = false)
    res = Pkn_NGG_arb(k, n,  β, σ)
    if has_reasonable_precision(res)
        return res
    else
        start_Pkn_val = res
        i = 1
        start_k_ind = undef
        while !has_reasonable_precision(start_Pkn_val)
            if 2^i*k>n
                error("There seem to be a numerical problem with computing Pkn even for large values of k")
            end
            start_k_ind = 2^i*k
            start_Pkn_val = Pkn_NGG_arb(start_k_ind, n,  β, σ)
            i += 1
        end
        if verbose
            println("k=$k, start_k_ind=$start_k_ind")
        end
        return Pkn_NGG_approx(k, n,  β, σ, start_k_ind, start_Pkn_val)
    end
end

function log_βnk(β, n, k, σ)
    return log(β) + log(n) - 1/σ*log(k)
end
βnk(β, n, k, σ) = exp(log_βnk(β, n, k, σ))

function logxk(n, k, β, σ)
    if n==1
        return log(k*σ ) + log(Cnk(n, k+1, σ)) - log(σ) - log(Cnk(n, k, σ))
    else
       return log(k*σ + βnk(β, n-1, k, σ)) + log(Cnk(n, k+1, σ)) - log(σ) - log(Cnk(n, k, σ))
    end   
end
## change: added logxk1 function
function logxk1(n, k, β, σ)
    return log(k*σ) + log(Cnk(n, k+1, σ)) - log(σ) - log(Cnk(n, k, σ))
end
## change

function Pkn_NGG_approx(k, n, β, σ, start_k_ind, start_Pkn_val)
    if k==start_k_ind
        return start_Pkn_val
    else
        logPkn = log(start_Pkn_val)
        for ki in (start_k_ind-1):-1:k
            # println(logPkn)
            logPkn = logPkn - logxk(n, ki, β, σ)
        end
    #     return P1kn
    #     [logprod_xi(n, j, β, σ) for j in 1:(n-1)]
        return exp(logPkn)
    end
end

Pkn_NGG_approx(k, n,  β, σ) = Pkn_NGG_approx(k, n,  β, σ, Int64(floor(n/2)), Pkn_NGG_arb(Int64(floor(n/2)), n,  β, σ))

function Pkn_NGG_approx(n, β, σ, start_k_ind, start_Pkn_val)
    logP1kn = Array{arb}(undef, start_k_ind)
    logP1kn[start_k_ind] = start_Pkn_val |> log
    for k in (start_k_ind-1):-1:1
        logP1kn[k] = logP1kn[k+1] - logxk(n, k, β, σ)
    end
#     return P1kn
#     [logprod_xi(n, j, β, σ) for j in 1:(n-1)]
    return exp.(logP1kn)
end

Pkn_NGG_approx(n, β, σ) = Pkn_NGG_approx(n, β, σ, n, Pkn_NGG_arb(n, n,  β, σ))

function Pkn_NGG_robust(n,  β, σ; verbose = false)
    P1n = Array{arb}(undef, n)
    int_n_half::Int64 = floor(n/2)
    P1n[int_n_half:n] = Pkn_NGG_arb.(int_n_half:n, n,  β, σ)
    Pkn_val = P1n[int_n_half]
    if !has_reasonable_precision(Pkn_val)
        error("There seem to be a numerical problem with computing Pkn even for k=n")
    else
        k = int_n_half
        while k ≥ 1 && has_reasonable_precision(Pkn_val)
            P1n[k] = Pkn_val
            k -= 1
            Pkn_val = Pkn_NGG_arb(k, n,  β, σ)
        end
        if verbose
            println("k=$k, n=$n")
        end
        if k > 1
            P1n[1:(k+1)] = Pkn_NGG_approx(n, β, σ, k+1, P1n[k+1])
        end
        return P1n
    end
end

Pkn_NGG(k, n,  β, σ) = Pkn_NGG_robust(k, n,  β, σ; verbose = false) |> Float64

## change added Pkn_NGG_approx_full
function Pkn_NGG_approx_full(n, β, σ, f)
        Axnk = Array{arb}(undef, n-1)
        Axnk[1:n-1] = f.(n, 1:n-1, β, σ)
        Sum_xn= f(n, 1, β, σ)
        for i in (3:n)
             Sum_xn= Sum_xn + exp(sum(Axnk[1:i-1]))
        end
        P1n = Array{arb}(undef, n)
        P1n[1] = exp(- log1p(Sum_xn))
        for k in 2:(n)
              P1n[k] = exp(log(P1n[k-1]) + Axnk[k-1])
        end 
        return P1n
end


Pkn_NGG_full_approximation(n, β, σ, f) = convert(Array{Float64,1}, Pkn_NGG_approx_full(n,  β, σ, f))


function Pkn_NGG_fapprox(k,n, β, σ, f)
        Axnk = Array{arb}(undef, n-1)
        Axnk[1:n-1] = f.(n, 1:n-1, β, σ)
        Sum_xn= f(n, 1, β, σ)
        for i in (3:n)
             Sum_xn= Sum_xn + exp(sum(Axnk[1:i-1]))
        end
        Pkn = exp(- log1p(Sum_xn))
        if k==1 
           return exp(- log1p(Sum_xn))
        else 
           for k in 2:(n)
              Pkn = exp(log(Pkn) + Axnk[k-1])
           end
        end
      return Pkn
end


Pkn_NGG_fapproximation(k, n, β, σ, f) =  Pkn_NGG_fapprox(k,n,  β, σ, f)  |> Float64



function Pkn_NGG_approx_partial(k, n, β, σ, start_k_ind, start_Pkn_val, f)
    if k==start_k_ind
        return start_Pkn_val
    else
        logPkn = log(start_Pkn_val)
        for ki in (start_k_ind-1):-1:k
            # println(logPkn)
            logPkn = logPkn - f(n, ki, β, σ)
        end
    #     return P1kn
    #     [logprod_xi(n, j, β, σ) for j in 1:(n-1)]
        return exp(logPkn)
    end
end

Pkn_NGG_approx(k, n,  β, σ, f) = Pkn_NGG_approx_partial(k, n,  β, σ, n, Pkn_NGG_arb(n, n,  β, σ),f) |> Float64

Pkn_NGG_raw(k, n,  β, σ) = Pkn_NGG_arb(k, n,  β, σ) |> Float64

## change
function Pkn_2PD_arb(k::N, n::N, θ::T, σ::T) where {T<:Number, N<:Integer}
    σ_arb = RR(σ)
    Vnk_2PD(n, k, θ, σ) // σ_arb^k * Cnk(n, k, σ)
end

Pkn_2PD(k::N, n::N, θ::T, σ::T) where {T<:Number, N<:Integer} = convert(Float64, Pkn_2PD_arb(k, n, θ, σ))
Pkn_PY_arb = Pkn_2PD_arb
Pkn_PY = Pkn_2PD

function Pkn_Dirichlet_arb(k, n,  θ)
    θ_arb = RR(θ)
    return θ_arb^k // risingfac(θ_arb,n) * unsigned_Stirling1(n,k)
end

Pkn_Dirichlet(k, n,  θ) = convert(Float64, Pkn_Dirichlet_arb(k, n,  θ))

## change

function Pkn_PDM_arb(k::N, n::N,H::N, θ::T, σ::T) where {T<:Number, N<:Integer}
    H_arb = RR(H)
    θ_arb = RR(θ)
    return (fac(H_arb) // (fac(H_arb-RR(k)) *risingfac(1+θ_arb,n-1))) * sum([ (1//(H_arb^i))  * (gamma(i + θ_arb //RR(σ))//gamma(1 + θ_arb //RR(σ))) * unsigned_Stirling2(i,k)* Cnk(n, i, σ)  for i in k:n ])    
end

Pkn_PYM(k, n,H,θ, σ) = convert(Float64, Pkn_PDM_arb(k n,H,θ, σ))


function Pkn_NGGM_arb(k::N, n::N,H::N,  β::T, σ::T) where {T<:Number, N<:Integer}
    H_arb = RR(H)
    return (fac(H_arb) // (H_arb^k * fac(H_arb-RR(k)))) * sum([ (1//(H_arb^i)) * unsigned_Stirling2(i+k,k)* Pkn_NGG_arb(k+i, n,  β, σ)  for i in 0:(n-k)])    
end

Pkn_NGGM(k, n,H,β, σ) = convert(Float64, Pkn_NGGM_arb(k, n, H, β, σ))


## change




