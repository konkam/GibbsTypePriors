![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)<!--
![Lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-retired-orange.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-archived-red.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-dormant-blue.svg) -->
[![Build Status](https://travis-ci.com/konkam/GibbsTypePriors.jl.svg?branch=master)](https://travis-ci.com/konkam/GibbsTypePriors.jl)
[![codecov.io](http://codecov.io/github/konkam/GibbsTypePriors.jl/coverage.svg?branch=master)](http://codecov.io/github/konkam/GibbsTypePriors.jl?branch=master)
<!--
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://konkam.github.io/GibbsTypePriors.jl/stable)
[![Documentation](https://img.shields.io/badge/docs-master-blue.svg)](https://konkam.github.io/GibbsTypePriors.jl/dev)
-->

Computing prior number of clusters for Gibbs-type priors.

## Introduction


The following reference gives an overview of Gibbs-type priors and their importance for Bayesian Nonparametrics:

De Blasi, Pierpaolo, Stefano Favaro, Antonio Lijoi, Ramsés H. Mena, Igor Prünster, and Matteo Ruggiero. “Are Gibbs-Type Priors the Most Natural Generalization of the Dirichlet Process?” IEEE Transactions on Pattern Analysis and Machine Intelligence 37, no. 2 (2015): 212–29. https://doi.org/10.1109/TPAMI.2013.217.


## How to install the package

**The package is developed for Julia 1.4.**
Press `]` in the Julia interpreter to enter the Pkg mode and input:

````julia

pkg> add https://github.com/konkam/GibbsTypePriors
````




# How to use the package

To compute the prior density at clusters of size k=10 for a Normalized Generalised Gamma process of parameters σ=0.8, β = 1.2 and n = 500 data points, use:


````julia
using GibbsTypePriors
Pkn(10, 500, 1.2, 0.8)
````


````
[8.984618037609137534209230272655590048020578883093909560975329204663764947
616454323239936826977747686166220699442638102575108090066733474045527410360
450121434830287551953465345194775537869809599211765994321637396240149972171
546510445613750582980320171624333628510225607407143099775585208996543023652
676275861806445460239812676123992225439009997501210839661499955929140859243
745115864620048499707015230200347468860803402195237130937860227937282435173
855624120923262777166460777455949996981610826766461852484489607126252210171
672699398019414520646536125750904460090434261809442516206810344875167684776
472999152201534929920295427711282514061362671348736126531428384753254537233
289157853096242587342910283218682276937685856894635612746655801909217588281
650331471076315910559666588597344197273223630894702054703628015377112870204
902412500493155028411534866762939881293404300772103092549550936499643327178
602002033975387284310325784102934080965885382382461440834753740599339997373
737338566134017795432463249262848651200208718606010830508789242452518351244
608176504311324542498794171034314733779974806974996196548907766596608340725
650997893905413817699024440849869052977198745935404683294056090065267373104
764844213760999093374793690496323646281544305467293457511936906685237578510
33910215691307143433860269144227930197494339073e-11 +/- 3.02e-1331]
````


