expected_number_of_cluster_2PD(n::N, θ::T, σ::T, ntrunc::N) where {T<:Number, N<:Integer} = Pkn_2PD_arb.(1:ntrunc, n, θ, σ) |> ar -> map(*, ar, 1:ntrunc) |> sum

expected_number_of_cluster_stable(n, σ, ntrunc) = expected_number_of_cluster_2PD(n, 0., σ, ntrunc)

function expected_number_of_clusters_Dirichlet(n::Int64, theta::Float64, ntrunc::Int64)
    return Pkn_Dirichlet_arb.(1:ntrunc, n, theta) |> ar -> map(*, ar, 1:ntrunc) |> sum
end
