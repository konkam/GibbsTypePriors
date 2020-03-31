using DataFrames

function expand_grid(iters...)
    vals = Base.Iterators.product(iters...)
    # l = Base.Iterators.take(vals, 1) |> collect
    reduce(vcat, DataFrame([v]) for v in vals)
end
