"""
    Pkn_NGG_arb(k, n, β, σ)

Compute the prior probability mass to obtain k clusters through direct computation with arbitrary precision. Not the most stable version.

# Examples
```julia-repl
julia> GibbsTypePriors.Pkn_NGG_arb(50, 100, 0.5, 0.2)
[2.14762363546438648278532704280127812426643249551735595554664760044538213621775301015485328374549413897540160511465187297274194698946624496309582378995095872169826130349407122713909436452654024908060179841056741955863382744781148806296586030164131604285251938875060389310445972984208190921903755337680167605513750214014651799055934001372969787136894687051060898423210469539540879653779017239230826328548585098581761936756315153017814207498994142860020958294285847670933928577444047966837019301853075052805215643829733023817130541787687113481647039349350666819737212357118184059914081469496628351252560981171967614071257194538987858556850738404955654861591001195199137908090799330617531307554849423121028628550128016847465035159250895414945842092857092939501917174224657462412244406519738404915932075136518078258780668230018325226513272503217141555634764038959108448484312201406445495467152388854367276670509646840248971542665576196296601294567533128012684921641261767619308293743726390603225250672867733057200093764412465619355633068601123653038002743467832229736915692389253166213741137134870700536044628527732021542164522846984454807319649387483162094391527262191665810862006754021341268115958393216764683235046467243249773739254152551430288889983085214706056491738375963632435912736357674313108409758462287962329443736179207688980017178194117395394271454788659703559151681296053952789173258126586683260147689572558384333967140709303631030995994161910062213092853666415615122632e-19 +/- 4.26e-1498]
```
"""
function Pkn_NGG_arb(k, n, β, σ)
    if k>n || k==0
        return arb_0
    else
        σ_arb = RR(σ)
        Vnk_NGG(n, k, β, σ) // σ_arb^k * Cnk(n, k, σ)
    end
end

"""
    Pkn_NGG_robust_in(k, n, β, σ)

Compute the prior probability mass to obtain k clusters, either through direct computation or using the recursive computation for the Cnk, with arbitrary precision. Stable for large values of k.

# Examples
```julia-repl
julia> GibbsTypePriors.Pkn_NGG_robust_in(50, 100, 0.5, 0.2)
[2.14762363546438648278532704280127812426643249551735595554664760044538213621775301015485328374549413897540160511465187297274194698946624496309582378995095872169826130349407122713909436452654024908060179841056741955863382744781148806296586030164131604285251938875060389310445972984208190921903755337680167605513750214014651799055934001372969787136894687051060898423210469539540879653779017239230826328548585098581761936756315153017814207498994142860020958294285847670933928577444047966837019301853075052805215643829733023817130541787687113481647039349350666819737212357118184059914081469496628351252560981171967614071257194538987858556850738404955654861591001195199137908090799330617531307554849423121028628550128016847465035159250895414945842092857092939501917174224657462412244406519738404915932075136518078258780668230018325226513272503217141555634764038959108448484312201406445495467152388854367276670509646840248971542665576196296601294567533128012684921641261767619308293743726390603225250672867733057200093764412465619355633068601123653038002743467832229736915692389253166213741137134870700536044628527732021542164522846984454807319649387483162094391527262191665810862006754021341268115958393216764683235046467243249773739254152551430288889983085214706056491738375963632435912736357674313108409758462287962329443736179207688980017178194117395394271454788659703559151681296053952789173258126586683260147689572558384333967140709303631030995994161910062213092853666415615122632e-19 +/- 4.26e-1498]
```
"""
function Pkn_NGG_robust_in(k, n, β, σ; verbose = false)
    if k>n || k==0
        return arb_0
    else
        σ_arb = RR(σ)
        Vnk_NGG_rec(n, k, β, σ) // σ_arb^k * Cnk_robust(n, k, σ; verbose = verbose)
    end
end

"""
    Pkn_NGG_robust(k, n, β, σ)

Compute the prior probability mass to obtain k clusters, either through direct computation or using the recursive computation for the Cnk, with arbitrary precision. Stable for large values of k. If this became unstable for small values of k, the function will print a warning and resort to using an approximation scheme based on the predictive distribution of Gibbs-type processes (Bystrova, 2020)

References:
- Bystrova, D, Kon Kam King, G, Deslandes, F, Arbel, J. Approximating the clusters' prior distribution in Bayesian nonparametric models, 2020 (in preparation).

# Examples
```julia-repl
julia> GibbsTypePriors.Pkn_NGG_robust(1, 100, 1.2, 0.6)
[1.9296391903064605384654157430472418076148726737727875230129641023067313090166540244753747312457129183657025493616396554177981937762780521371044832530798250427147509991595690117491772034391001735237273292732682815399768750800178863279517168277700486277214486303628778846185502064952747248507331917928973821847721861122917899722304242486554537125500815028027048749342590041556047806174326663544410878820833633145820033554401795667161239326297269828805842051466014698678958849980095944598411211901947299487314690595562677257989440632872268784586850284818761100358647158324006198078594045133536472908933706208286326023740694795701448737505642410975566985814946675218483636408975980085854487052384785613476386696535340182246327012856942430959614307701381885382101908730464363930653564743180041556516459003976577914678253451116729865709223986356119211820112193877506211938987401933281878273977621027890553237363870948374184877771847370494105304641393041184948881280139977645736656738738876663695255598855102890090929354876451974185700013265883795024142293167970426325083355312585087826163789033506198771766732204008615074408731001784408725973789012051139760486038536969512748787960733295275525464396278674073554207828665628723052100436461108870416804992391459295951069501611168607219137361112920188489260042042387173028715544129182561207314706608564475318060162678254188901156369626393351965452285428586670794557016665659781242771755478946260451052238e-6 +/- 5.76e-1451]
```
"""
function Pkn_NGG_robust(k, n, β, σ; verbose = false)
    res = Pkn_NGG_robust_in(k, n, β, σ; verbose = verbose)
    if has_reasonable_precision(res)
        return res
    else
        @warn "Direct computation not accurate enough, resorting to approximation"
        start_Pkn_val = res
        i = 1
        start_k_ind = undef
        while !has_reasonable_precision(start_Pkn_val)
            if 2^i*k>n
                error("There seem to be a numerical problem with computing Pkn even for large values of k")
            end
            start_k_ind = 2^i*k
            start_Pkn_val = Pkn_NGG_robust_in(start_k_ind, n, β, σ; verbose = false)#verbose = false bccause this can be very informative
            i += 1
        end
        if verbose
            println("k=$k, start_k_ind=$start_k_ind")
        end
        return Pkn_NGG_approx(k, n, β, σ, start_k_ind, start_Pkn_val)
    end
end

"""
    log_βnk(β, n, k, σ)

log-factor involved in the expansion of the predictive distribution (Bystrova, 2020).

References:
- Bystrova, D, Kon Kam King, G, Deslandes, F, Arbel, J. Approximating the clusters' prior distribution in Bayesian nonparametric models, 2020 (in preparation).
```
"""
function log_βnk(β, n, k, σ)
    return log(β) + log(n) - 1/σ*log(k)
end
"""
    βnk(β, n, k, σ)

    factor involved in the expansion of the predictive distribution (Bystrova, 2020).

    References:
    - Bystrova, D, Kon Kam King, G, Deslandes, F, Arbel, J. Approximating the clusters' prior distribution in Bayesian nonparametric models, 2020 (in preparation).
```
"""
βnk(β, n, k, σ) = exp(log_βnk(β, n, k, σ))

"""
    logxk(n, k, β, σ)

log of the ratio p_{n+1,k+1}/p_{n+1,k} (Bystrova, 2020).

    References:
    - Bystrova, D, Kon Kam King, G, Deslandes, F, Arbel, J. Approximating the clusters' prior distribution in Bayesian nonparametric models, 2020 (in preparation).
```
"""
function logxk(n, k, β, σ)
    if n==1
        return log(k*σ ) + log(Cnk(n, k+1, σ)) - log(σ) - log(Cnk(n, k, σ))
    else
       return log(k*σ + βnk(β, n-1, k, σ)) + log(Cnk(n, k+1, σ)) - log(σ) - log(Cnk(n, k, σ))
    end
end


"""
    Pkn_NGG_approx(k, n, β, σ, start_k_ind, start_Pkn_val)

Compute the prior probability mass to obtain k clusters using an approximation scheme based on the predictive distribution of Gibbs-type processes (Bystrova, 2020). This starts from a hopefully exact value for large k and recursively computes values for smaller k.

References:
- Bystrova, D, Kon Kam King, G, Deslandes, F, Arbel, J. Approximating the clusters' prior distribution in Bayesian nonparametric models, 2020 (in preparation).

# Examples
```julia-repl
julia> GibbsTypePriors.Pkn_NGG_approx(3, 20, 1.2, 0.6, 4, 0.5)
```
"""
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

"""
    Pkn_NGG_approx(k, n, β, σ)

Compute the prior probability mass to obtain k clusters using an approximation scheme based on the predictive distribution of Gibbs-type processes (Bystrova, 2020). This starts from a hopefully exact value for k = n and recursively computes values for smaller k.

References:
- Bystrova, D, Kon Kam King, G, Deslandes, F, Arbel, J. Approximating the clusters' prior distribution in Bayesian nonparametric models, 2020 (in preparation).

# Examples
```julia-repl
julia> GibbsTypePriors.Pkn_NGG_approx(1, 20, 1.2, 0.6)
```
"""
Pkn_NGG_approx(k, n, β, σ) = Pkn_NGG_approx(k, n, β, σ, Int64(floor(n/2)), Pkn_NGG_arb(Int64(floor(n/2)), n, β, σ))

"""
    Pkn_NGG_approx(k, n, β, σ, start_k_ind, start_Pkn_val)

Compute the prior probability mass to obtain k clusters for all k between 1 and start_Pkn_val using an approximation scheme based on the predictive distribution of Gibbs-type processes (Bystrova, 2020). This starts from a hopefully exact value for k = n and recursively computes values for smaller k.

References:
- Bystrova, D, Kon Kam King, G, Deslandes, F, Arbel, J. Approximating the clusters' prior distribution in Bayesian nonparametric models, 2020 (in preparation).

# Examples
```julia-repl
julia> GibbsTypePriors.Pkn_NGG_approx(20, 1.2, 0.6)
```
"""
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

"""
    Pkn_NGG_approx(n, β, σ)

Compute the prior probability mass to obtain k clusters for all k between 1 and n using an approximation scheme based on the predictive distribution of Gibbs-type processes (Bystrova, 2020). This starts from a hopefully exact value for k = n and recursively computes values for smaller k.

References:
- Bystrova, D, Kon Kam King, G, Deslandes, F, Arbel, J. Approximating the clusters' prior distribution in Bayesian nonparametric models, 2020 (in preparation).

# Examples
```julia-repl
julia> GibbsTypePriors.Pkn_NGG_approx(20, 1.2, 0.6)
```
"""
Pkn_NGG_approx(n, β, σ) = Pkn_NGG_approx(n, β, σ, n, Pkn_NGG_arb(n, n, β, σ))


function Pkn_NGG_robust(n, β, σ; verbose = false)
    P1n = Array{arb}(undef, n)
    int_n_half::Int64 = floor(n/2)
    P1n[int_n_half:n] = Pkn_NGG_arb.(int_n_half:n, n, β, σ)
    Pkn_val = P1n[int_n_half]
    if !has_reasonable_precision(Pkn_val)
        error("There seem to be a numerical problem with computing Pkn even for k=n/2")
    else
        k = int_n_half
        while k ≥ 1 && has_reasonable_precision(Pkn_val)
            P1n[k] = Pkn_val
            k -= 1
            Pkn_val = Pkn_NGG_arb(k, n, β, σ)
        end
        if verbose
            println("k=$k, n=$n")
        end
        if k > 1
            @warn "Direct computation not accurate enough, resorting to approximation"
            P1n[1:(k+1)] = Pkn_NGG_approx(n, β, σ, k+1, P1n[k+1])
        end
        return P1n
    end
end

"""
    Pkn_NGG(k, n, β, σ)

Compute the prior probability mass to obtain k clusters, either through direct computation or using the recursive computation for the Cnk, with arbitrary precision. Stable for large values of k, the function may become unstable for small values of k. If this were the case, the function would print a warning and resort to using an approximation scheme based on the predictive distribution of Gibbs-type processes (Bystrova, 2020)

References:
- Bystrova, D, Kon Kam King, G, Deslandes, F, Arbel, J. Approximating the clusters' prior distribution in Bayesian nonparametric models, 2020 (in preparation).

# Examples
```julia-repl
julia> GibbsTypePriors.Pkn_NGG(10, 100, 1.2, 0.6)
```
"""
Pkn_NGG(k, n, β, σ) = Pkn_NGG_robust(k, n, β, σ; verbose = false) |> Float64

"""
    Pkn_NGG(n, β, σ)

Compute the prior probability mass to obtain k clusters over k=0 to k=n, either through direct computation or using the recursive computation for the Cnk, with arbitrary precision. Stable for large values of k. If this became unstable for small values of k, the function will print a warning and resort to using an approximation scheme based on the predictive distribution of Gibbs-type processes (Bystrova, 2020)

References:
- Bystrova, D, Kon Kam King, G, Deslandes, F, Arbel, J. Approximating the clusters' prior distribution in Bayesian nonparametric models, 2020 (in preparation).

# Examples
```julia-repl
julia> GibbsTypePriors.Pkn_NGG(10, 100, 1.2, 0.6)
```
"""
Pkn_NGG(n, β, σ) = Pkn_NGG_robust(k, n, β, σ; verbose = false) |> Float64

function Pkn_2PD_arb(k::N, n::N, θ::T, σ::T) where {T<:Number, N<:Integer}
    σ_arb = RR(σ)
    Vnk_2PD(n, k, θ, σ) // σ_arb^k * Cnk(n, k, σ)
end


"""
    Pkn_2PD(n, β, σ)

Compute the prior probability mass to obtain k clusters for a 2-parameter Poisson-Dirichlet process (also known as Pitman-Yor). 

# Examples
```julia-repl
julia> Pkn_2PD(10, 100, 1.2, 0.6)
0.00787749555803404
```
"""
Pkn_2PD(k::N, n::N, θ::T, σ::T) where {T<:Number, N<:Integer} = convert(Float64, Pkn_2PD_arb(k, n, θ, σ))
Pkn_PY_arb = Pkn_2PD_arb
"""
    Pkn_PY(n, β, σ)

Compute the prior probability mass to obtain k clusters for a 2-parameter Poisson-Dirichlet process (also known as Pitman-Yor). 

# Examples
```julia-repl
julia> Pkn_PY(10, 100, 1.2, 0.6)
0.00787749555803404
```
"""
Pkn_PY = Pkn_2PD


function Pkn_Dirichlet_arb(k, n, θ)
    θ_arb = RR(θ)
    return θ_arb^k // risingfac(θ_arb,n) * unsigned_Stirling1(n,k)
end

"""
    Pkn_Dirichlet(k, n, θ)

Compute the prior probability mass to obtain k clusters for a Dirichlet process. 

# Examples
```julia-repl
julia> Pkn_Dirichlet(7, 10, 1.5)
```
"""
Pkn_Dirichlet(k, n, θ) = convert(Float64, Pkn_Dirichlet_arb(k, n, θ))
