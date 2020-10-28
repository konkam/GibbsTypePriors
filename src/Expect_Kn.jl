"""
    expected_number_of_cluster_2PD(n::N, θ::T, σ::T, ntrunc::N)

Compute the expected number of clusters for a 2-parameter Poisson Dirichlet prior, by direct computation up to ntrunc. If ntrunc = n, this is an exact computation, if it turns out to be too long and that you expect large values of k to have a very small mass, you can set ntrunc < n.

# Examples
```julia-repl
julia> GibbsTypePriors.unsigned_Stirling1(10, 5)
269325
```
"""
expected_number_of_cluster_2PD(n::N, θ::T, σ::T, ntrunc::N) where {T<:Number, N<:Integer} = Pkn_2PD_arb.(1:ntrunc, n, θ, σ) |> ar -> map(*, ar, 1:ntrunc) |> sum

expected_number_of_cluster_2PD(n::N, θ::T, σ::T) where {T<:Number, N<:Integer} = expected_number_of_cluster_2PD(n, θ, σ, n)

expected_number_of_cluster_stable(n::N, σ::T, ntrunc::N) where {T<:Number, N<:Integer} = expected_number_of_cluster_2PD(n, 0., σ, ntrunc)
expected_number_of_cluster_stable(n, σ) = expected_number_of_cluster_stable(n, σ, n)


function expected_number_of_clusters_Dirichlet(n::Int64, theta::Float64, ntrunc::Int64)
    return Pkn_Dirichlet_arb.(1:ntrunc, n, theta) |> ar -> map(*, ar, 1:ntrunc) |> sum
end

expected_number_of_clusters_Dirichlet(n::Int64, theta::Float64) = expected_number_of_clusters_Dirichlet(n, theta, n)


function expected_number_of_clusters_Dirichlet_Multinomial(n::Int64, N::Int64, theta::Float64, ntrunc::Int64)
    return Pkn_Dirichlet_Mult_arb.(1:ntrunc, n,N, theta) |> ar -> map(*, ar, 1:ntrunc) |> sum
end

expected_number_of_clusters_Dirichlet_Multinomial(n::Int64, N::Int64, theta::Float64) = expected_number_of_clusters_Dirichlet_Multinomial(n, N, theta, min(N,n))

function expected_number_of_clusters_NGG(n::Int64, β::Float64, σ::Float64, ntrunc::Int64)
    return Pkn_NGG.(1:ntrunc, n, β, σ) |> ar -> map(*, ar, 1:ntrunc) |> sum
end

expected_number_of_clusters_NGG(n::Int64, β::Float64, σ::Float64) = expected_number_of_clusters_NGG(n, β, σ,n)


function expected_number_of_clusters_NGG_mult(n::Int64, H::Int64, β::Float64, σ::Float64, ntrunc::Int64)
    pk = Pkn_NGG_mult(n, H, β, σ)
    return pk[ 1:ntrunc] |> ar -> map(*, ar, 1:ntrunc) |> sum
end


expected_number_of_clusters_NGG_mult(n::Int64, H::Int64, β::Float64, σ::Float64) = expected_number_of_clusters_NGG_mult(n, H, β, σ,n)

#aliases

"""
    E_2PD(n::N, θ::T, σ::T, ntrunc::N) where {T<:Number, N<:Integer}

Compute the expected number of clusters for a 2-parameter Poisson Dirichlet (also known as Pitman-Yor) prior, by direct computation up to ntrunc (optional argument). If ntrunc = n, this is an exact computation, if it turns out to be too long and that you expect large values of k to have a very small mass, you can set ntrunc < n.

See also: [`E_stable`](@ref), [`E_Dirichlet`](@ref)


# Examples
```julia-repl
julia> E_2PD(10, 0., 0.4)
[2.797180108800000146227991769541811614624157201057103170680445069708107636422615452406497920159297853634382963614055479173158069186418337497999385184360757154429265052398107259805621680470326929725738226729222577234589857779643629658552208848284173508026015152904404059960898279167641852826517393043262627617296658490156505882278711266617066316092443733068941404176241441191328628518109412653612525798394207094225161730723863943059051559749653061842655787927469646092504262924194335937500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 +/- 5.60e-1504]

```
"""
E_2PD = expected_number_of_cluster_2PD
"""
    E_PY(n::N, θ::T, σ::T, ntrunc::N) where {T<:Number, N<:Integer}

Compute the expected number of clusters for a 2-parameter Poisson Dirichlet (also known as Pitman-Yor) prior, by direct computation up to ntrunc (optional argument). If ntrunc = n, this is an exact computation, if it turns out to be too long and that you expect large values of k to have a very small mass, you can set ntrunc < n.

See also: [`E_stable`](@ref), [`E_Dirichlet`](@ref)


# Examples
```julia-repl
julia> E_PY(10, 0., 0.4)
[2.797180108800000146227991769541811614624157201057103170680445069708107636422615452406497920159297853634382963614055479173158069186418337497999385184360757154429265052398107259805621680470326929725738226729222577234589857779643629658552208848284173508026015152904404059960898279167641852826517393043262627617296658490156505882278711266617066316092443733068941404176241441191328628518109412653612525798394207094225161730723863943059051559749653061842655787927469646092504262924194335937500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 +/- 5.60e-1504]
```
"""
E_PY = expected_number_of_cluster_2PD
"""
    E_stable(n::N, σ::T, ntrunc::N) where {T<:Number, N<:Integer}

Compute the expected number of clusters for a stable process prior, by direct computation up to ntrunc (optional argument). If ntrunc = n, this is an exact computation, if it turns out to be too long and that you expect large values of k to have a very small mass, you can set ntrunc < n.

See also: [`E_2PD`](@ref), [`E_PY`](@ref), [`E_Dirichlet`](@ref)


# Examples
```julia-repl
julia> E_stable(10, 0.1)
[1.317283550395312519640862117090075082018602293487207173467704662870901710848586072338941448910368160445217137507034833714840735500227050610150306349069777704715852016445662205119280375270463908555026018111545594765430207866848367807193962524362807035164072216472928669667852836537566353834476540485547750703723020476001586935706322137304949274062823739150683831653442615917431013917124081952539061271081947740823646018757470139485514508444060324299094553588058194998344774262477585580199956893920898437500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 +/- 3.11e-1504]
```
"""
E_stable = expected_number_of_cluster_stable
"""
    E_Dirichlet(n::Int64, theta::Float64, ntrunc::Int64) where {T<:Number, N<:Integer}

Compute the expected number of clusters for a Dirichlet process prior, by direct computation up to ntrunc (optional argument). If ntrunc = n, this is an exact computation, if it turns out to be too long and that you expect large values of k to have a very small mass, you can set ntrunc < n.

See also: [`E_2PD`](@ref), [`E_PY`](@ref), [`E_Dirichlet`](@ref)


# Examples
```julia-repl
julia> E_Dirichlet(10, 0.1)
[1.317283550395312519640862117090075082018602293487207173467704662870901710848586072338941448910368160445217137507034833714840735500227050610150306349069777704715852016445662205119280375270463908555026018111545594765430207866848367807193962524362807035164072216472928669667852836537566353834476540485547750703723020476001586935706322137304949274062823739150683831653442615917431013917124081952539061271081947740823646018757470139485514508444060324299094553588058194998344774262477585580199956893920898437500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 +/- 3.11e-1504]
```
"""
E_Dirichlet = expected_number_of_clusters_Dirichlet

"""
    E_Dirichlet(n::Int64, theta::Float64, ntrunc::Int64) where {T<:Number, N<:Integer}

Compute the expected number of clusters for a Dirichlet process prior, by direct computation up to ntrunc (optional argument). If ntrunc = n, this is an exact computation, if it turns out to be too long and that you expect large values of k to have a very small mass, you can set ntrunc < n.

See also: [`E_2PD`](@ref), [`E_PY`](@ref), [`E_Dirichlet`](@ref)


# Examples
```julia-repl
julia> E_Dirichlet_multinomial(10, 0.1, 100)
[1.317283550395312519640862117090075082018602293487207173467704662870901710848586072338941448910368160445217137507034833714840735500227050610150306349069777704715852016445662205119280375270463908555026018111545594765430207866848367807193962524362807035164072216472928669667852836537566353834476540485547750703723020476001586935706322137304949274062823739150683831653442615917431013917124081952539061271081947740823646018757470139485514508444060324299094553588058194998344774262477585580199956893920898437500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 +/- 3.11e-1504]
```
"""
E_Dirichlet_multinomial  = expected_number_of_clusters_Dirichlet_Multinomial
"""
    E_Dirichlet_multinomial(n::Int64, H::Int64, theta::Float64, ntrunc::Int64)

Compute the expected number of clusters for a Dirichlet multinomial process prior for number of support points H, by direct computation up to ntrunc (optional argument). If ntrunc = min(n,H), this is an exact computation, if it turns out to be too long and that you expect large values of k to have a very small mass, you can set ntrunc < min(n,H).

See also: [`E_2PD`](@ref), [`E_PY`](@ref), [`E_Dirichlet`](@ref)


# Examples
```julia-repl
julia> E_Dirichlet_multinomial(10, 0.1, 100)

```
"""
E_NGG  = expected_number_of_clusters_NGG
"""
    E_NGG(n::Int64, beta::Float64, sigma::Float64, ntrunc::Int64)

Compute the expected number of clusters for a normalized generalized gamma process prior, by direct computation up to ntrunc (optional argument). If ntrunc = n, this is an exact computation, if it turns out to be too long and that you expect large values of k to have a very small mass, you can set ntrunc < n.

See also: [`E_2PD`](@ref), [`E_PY`](@ref), [`E_Dirichlet`](@ref)


# Examples
```julia-repl
julia> E_NGG(10, 0.1, 0.2)

```
"""

E_NGG_multinomial  = expected_number_of_clusters_NGG_mult
"""
    E_NGG_multinomial(n::Int64, H:: Int64, beta::Float64, sigma::Float64, ntrunc::Int64)

Compute the expected number of clusters for a normalized generalized gamma multinomial process prior. If ntrunc = n, this is an exact computation, if it turns out to be too long and that you expect large values of k to have a very small mass, you can set ntrunc < n.

See also: [`E_2PD`](@ref), [`E_PY`](@ref), [`E_Dirichlet`](@ref)


# Examples
```julia-repl
julia> E_NGG_multinomial(10, 100, 0.1, 0.2)
 1.7702496359209865

```
"""




function E_PD_exact(n::N, θ::T, σ::T)::T where {T<:Number, N<:Integer}
    return (θ/σ)*(prod([(θ + σ + i -1)/(θ + i -1) for i in 1:n]) -1)
end

E_PY_exact  = E_PD_exact
"""
    E_PY_exact(n::Int64, theta::Float64, sigma::Float64)

Compute the expected number of clusters for Pitman Yor process prior, by exact computation.

See also: [`E_2PD`](@ref), [`E_PY`](@ref), [`E_Dirichlet`](@ref)


# Examples
```julia-repl
julia> E_PY_exact(10, 0.1, 0.2)

```
"""
function Var_PD_exact(n::N, θ::T, σ::T)::T where {T<:Number, N<:Integer}
           if n==1
               return 0
           else
               E_prev =  E_PD_exact(n-1, θ, σ)
               add_term  =  (E_prev*((n-1)*σ - θ*σ) + (n-1)*θ -  σ*σ*(E_prev^2))/ (n- 1 + θ)^2
               return    Var_PD_exact(n-1, θ, σ)*(n - 1 + θ + 2*σ)/ (n - 1 + θ) + add_term
           end
       end

V_PY_exact  = Var_PD_exact
"""
   V_PY_exact(n::Int64, theta::Float64, sigma::Float64)

Compute the expected number of clusters for Pitman Yor process prior, by exact computation.

See also: [`E_2PD`](@ref), [`E_PY`](@ref), [`E_Dirichlet`](@ref)

# Examples
```julia-repl
julia> V_PY_exact(10, 1.0, 0.2)
2.2826090179905956
```
"""
