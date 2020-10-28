using Distributions, FreqTables, DataFrames

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

function MvInv(eps, u, alpha, beta, gama, N, M)
    x = -log.(range(exp(-1e-05),exp(-10), length=N))
    f =  alpha/ Float64(real(gamma(CC(1- gama)))) * (x.^(-(1+gama)).*exp.( - (u+beta)*x))
    dx = diff(x)
    h = (f[Not(1)] + f[Not(N)]) / 2
    Mv = append!(reverse(cumsum(reverse(dx.* h))),[0])
    dist = Exponential(1)
    W = rand(dist, M)
    W = cumsum(W)
    return fill_v2(M,Mv,W,N,x)
end


function TruncatedNGG_test(alpha, beta, sigma, Nw,Meps)
    Js = MvInv(Meps, 0,alpha, beta, sigma,50001,Nw)
    sum_Js =0.0
    m2_JS = 0.0
    m3_JS = 0.0
    for i in 1:10000
        Js = MvInv(Meps, 0,alpha, beta, sigma,50001,Nw)
        sum_Js = sum_Js + sum(Js)
        m2_JS = m2_JS + (sum(Js))^2
        m3_JS = m3_JS + (sum(Js))^3
    end
    println([sum_Js/10000, alpha])
    println([m2_JS/10000, alpha^2 + alpha*(1-sigma)])
    println([m3_JS/10000, alpha^3 + 3*alpha^2*(1-sigma) + alpha*risingfac(1-sigma,2)])
    Js = MvInv(Meps, 0,alpha, beta, sigma,50001,Nw)
    Js = Js/(sum(Js))
    return Js
end


function TruncatedNGG(alpha, beta, sigma, Nw,Meps)
    Js = MvInv(Meps, 0,alpha, beta, sigma,50001,Nw)
    Js = Js/(sum(Js))
    return Js
end

function prior_Kn(alpha,beta, sigma, N_s,N_tr ; runs=10^4)
    array_clust_num = Array{Int64}(undef, runs)
    for i in 1:runs
        weights_NGG =  TruncatedNGG(alpha,beta, sigma, N_tr, 0.01)
        c = wsample(1:N_tr, weights_NGG,N_s, replace=true)
        n_c=  length(Set(c))
        array_clust_num[i] = n_c
    end
    p_k_array = freqtable(array_clust_num)
    p_k = p_k_array./sum(p_k_array)
    df_pk= DataFrame(k = names(p_k,1),
                p_k = p_k)
    df_pk_zeros=DataFrame(k = collect(1:N_s)[collect(1:N_s).∉ Ref(df_pk.k)],
                            p_k = zeros(length(collect(1:N_s)[collect(1:N_s).∉ Ref(df_pk.k)])))
    df =  vcat(df_pk, df_pk_zeros)
    E_k = map(*, df.p_k, df.k) |> sum
    V_k = map(*, df.p_k,( df.k .- E_k).^2) |> sum
    return [df, E_k, V_k]
end

#A = prior_Kn(1.0,1.0, 0.25, 100, 100 ; runs=10000)

function exp_K_n_prior_NGG(sigma_ar,beta_ar, Ntr, Ns, it; alpha =1.0)
    #priors_array =  Array{Float64}(undef, length(sigma_ar))
    priors_array =Any[]
    for i in 1:length(sigma_ar)
        NGG_prior =  prior_Kn(alpha, beta_ar[i],sigma_ar[i], Ns, Ntr, runs = it)
        #priors_array[i] =  NGG_prior[1].p_k
        push!(priors_array,NGG_prior[1].p_k)
    end
    #save(priors_list, file=  paste0("FK_NGG_",sigma_array[1],sigma_array[length(sigma_array)],"beta_",beta_array[i],"Ns_",Ns,"N_tr_",Ntr,".Rda"))
    return priors_array
end
