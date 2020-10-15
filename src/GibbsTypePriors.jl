module GibbsTypePriors

export Pkn_NGG, Pkn_2PD, Pkn_PY, Pkn_Dirichlet, E_2PD, E_PY, E_stable, E_Dirichlet

include("common_functions.jl")
include("Cnk.jl")
include("Cnk_variable_prec.jl")
include("Vnk.jl")
include("Vnk_variable_prec.jl")
include("Pkn.jl")
include("Pkn_variable_prec.jl")
include("Expect_Kn.jl")
end # module
