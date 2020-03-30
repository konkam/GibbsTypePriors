function Pkn_NGG_arb(k, n,  β, σ)
    if k>n || k==0
        return RR(0)
    else
        σ_arb = RR(σ)
        Vnk_NGG(n, k, β, σ) // σ_arb^k * Cnk(n, k, σ)
    end
end

function has_reasonable_precision(arb_num)
    return accuracy_bits(arb_num) ≥ 64
end

function Pkn_robust(k, n,  β, σ; verbose = false)
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
        return Pkn_approx(k, n,  β, σ, start_k_ind, start_Pkn_val)
    end
end

function log_βnk(β, n, k, σ)
    return log(β) + log(n) - 1/σ*log(k)
end
βnk(β, n, k, σ) = exp(log_βnk(β, n, k, σ))
function logxk(n, k, β, σ)
    return log(k*σ + βnk(β, n, k, σ)) + log(Cnk(n, k+1, σ)) - log(σ) - log(Cnk(n, k, σ))
end

function Pkn_approx(k, n, β, σ, start_k_ind, start_Pkn_val)
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

Pkn_approx(k, n,  β, σ) = Pkn_approx(k, n,  β, σ, Int64(floor(n/2)), Pkn_NGG_arb(Int64(floor(n/2)), n,  β, σ))

function Pkn_approx(n, β, σ, start_k_ind, start_Pkn_val)
    logP1kn = Array{arb}(undef, start_k_ind)
    logP1kn[start_k_ind] = start_Pkn_val |> log
    for k in (start_k_ind-1):-1:1
        logP1kn[k] = logP1kn[k+1] - logxk(n, k, β, σ)
    end
#     return P1kn
#     [logprod_xi(n, j, β, σ) for j in 1:(n-1)]
    return exp.(logP1kn)
end

Pkn_approx(n, β, σ) = Pkn_approx(n, β, σ, n, Pkn_NGG_arb(n, n,  β, σ))

function Pkn_robust(n,  β, σ; verbose = false)
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
            P1n[1:(k+1)] = Pkn_approx(n, β, σ, k+1, P1n[k+1])
        end
        return P1n
    end
end

Pkn_NGG(k, n,  β, σ) = Pkn_robust(k, n,  β, σ; verbose = false)
