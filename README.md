![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)<!--
![Lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-retired-orange.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-archived-red.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-dormant-blue.svg) -->
[![Build Status](https://travis-ci.com/konkam/GibbsTypePriors.jl.svg?branch=master)](https://travis-ci.com/konkam/GibbsTypePriors.jl)
[![codecov.io](http://codecov.io/github/konkam/GibbsTypePriors.jl/coverage.svg?branch=master)](http://codecov.io/github/konkam/GibbsTypePriors.jl?branch=master)
[![Coverage Status](https://coveralls.io/repos/github/konkam/GibbsTypePriors/badge.svg?branch=master)](https://coveralls.io/github/konkam/GibbsTypePriors?branch=master)
<!--
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://konkam.github.io/GibbsTypePriors.jl/stable)
[![Documentation](https://img.shields.io/badge/docs-master-blue.svg)](https://konkam.github.io/GibbsTypePriors.jl/dev)
-->

Computing prior number of clusters for Gibbs-type priors.

## Introduction


The following reference gives an overview of Gibbs-type priors and their importance for Bayesian Nonparametrics:

De Blasi, Pierpaolo, Stefano Favaro, Antonio Lijoi, Ramsés H. Mena, Igor Prünster, and Matteo Ruggiero. “Are Gibbs-Type Priors the Most Natural Generalization of the Dirichlet Process?” IEEE Transactions on Pattern Analysis and Machine Intelligence 37, no. 2 (2015): 212–29. https://doi.org/10.1109/TPAMI.2013.217.


## How to install the package

The package is developed for Julia 1.4. Much of its functionality rests on the 'Arb' package[1], via its interface in `Nemo.jl`.

Press `]` in the Julia interpreter to enter the Pkg mode and input:

````julia

pkg> add https://github.com/konkam/GibbsTypePriors
````




# How to use the package

To compute the prior density at clusters of size k=10 for a Normalized Generalised Gamma process of parameters σ=0.8, β = 1.2 and n = 500 data points, use:


````julia
using GibbsTypePriors
Pkn_NGG(10, 500, 1.2, 0.8)
````


````
8.984618037609138e-11
````





The same may be done for the 2-parameter Poisson Dirichlet, also named the Pitman-Yor process:

````julia
Pkn_PY(10, 500, 1.2, 0.8)
````


````
2.562372640654159e-5
````





We also provide the same function for the Dirichlet process:

  ```julia
  Pkn_Dirichlet(10, 500, 1.2)
  ```

# Illustration of the various priors:

````julia
using GibbsTypePriors, DataFrames, DataFramesMeta, RCall
R"library(tidyverse)"

R"p = ggplot(data.frame(x = 1:50,
                        Pkn_NGG = $(Pkn_NGG.(1:50, 50,  48.4185, 0.25)),
                        Pkn_NGG2 = $(Pkn_NGG.(1:50, 50,  1., 0.7353)),
                        Pkn_Dirichlet = $(Pkn_Dirichlet.(1:50, 50,  19.233)),
                        Pkn_2PD = $(Pkn_2PD.(1:50, 50,  12.2157, 0.25)),
                        Pkn_2PD2 = $(Pkn_2PD.(1:50, 50,  1., 0.73001))
                    ) %>%
            gather(Process_type, density, Pkn_NGG:Pkn_2PD2),
    aes(x=x, y = density, colour = Process_type)) +
geom_line() +
ggthemes::scale_colour_ptol() +
theme_minimal()"
R"png('Illustration.png')
plot(p)
dev.off()"
````




 ![](Illustration.png)

References:
[1] Johansson, F. (2017).  Arb:  efficient arbitrary-precision midpoint-radius interval arithmetic.IEEE Transac-tions on Computers, 66:1281–1292.
