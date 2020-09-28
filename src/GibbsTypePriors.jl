module GibbsTypePriors

export Pkn_NGG, Pkn_2PD, Pkn_PY, Pkn_Dirichlet, E_2PD, E_PY, E_stable, E_Dirichlet

include("common_functions.jl")
include("Cnk.jl")
include("Vnk.jl")
include("Pkn.jl")
include("Expect_Kn.jl")
end # module
