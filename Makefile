
### Saves for graphical tests 1000

test/graphical_tests/saves_for_graphical_tests/accuracy_Vnk_1000.jld: test/graphical_tests/accuracy_Vnk_1000_cmp.jl
	julia --project=.  test/graphical_tests/accuracy_Vnk_1000_cmp.jl

test/graphical_tests/saves_for_graphical_tests/accuracy_Vnk_sigma_02_1000.jld: test/graphical_tests/accuracy_Vnk_sigma_02_1000_cmp.jl
	julia --project=.  test/graphical_tests/accuracy_Vnk_sigma_02_1000_cmp.jl

test/graphical_tests/saves_for_graphical_tests/accuracy_Pkn_1000.jld: test/graphical_tests/accuracy_Pkn_1000_cmp.jl
	julia --project=.  test/graphical_tests/accuracy_Pkn_1000_cmp.jl

test/graphical_tests/saves_for_graphical_tests/accuracy_Pkn_sigma_02_1000.jld: test/graphical_tests/accuracy_Pkn_sigma_02_1000_cmp.jl
		julia --project=.  test/graphical_tests/accuracy_Pkn_sigma_02_1000_cmp.jl

test/graphical_tests/saves_for_graphical_tests/accuracy_Cnk_1000.jld: test/graphical_tests/accuracy_Cnk_1000_cmp.jl
	julia --project=.  test/graphical_tests/accuracy_Cnk_1000_cmp.jl

test/graphical_tests/saves_for_graphical_tests/accuracy_Cnk_sigma_02_1000.jld: test/graphical_tests/accuracy_Cnk_sigma_02_1000_cmp.jl
	julia --project=.  test/graphical_tests/accuracy_Cnk_sigma_02_1000_cmp.jl

### Graph for graphical tests 1000

test/graphical_tests/figures_graphical_tests/accuracy_PknCnkVnk_1000.pdf: test/graphical_tests/saves_for_graphical_tests/accuracy_Vnk_1000.jld test/graphical_tests/saves_for_graphical_tests/accuracy_Pkn_1000.jld test/graphical_tests/saves_for_graphical_tests/accuracy_Cnk_1000.jld test/graphical_tests/saves_for_graphical_tests/accuracy_Vnk_sigma_02_1000.jld test/graphical_tests/saves_for_graphical_tests/accuracy_Pkn_sigma_02_1000.jld test/graphical_tests/saves_for_graphical_tests/accuracy_Cnk_sigma_02_1000.jld test/graphical_tests/accuracy_PknCnkVnk_1000_plt.jl
	julia --project=.  test/graphical_tests/accuracy_PknCnkVnk_1000_plt.jl

### Saves for graphical tests 10000

test/graphical_tests/saves_for_graphical_tests/accuracy_Vnk_10000.jld: test/graphical_tests/accuracy_Vnk_10000_cmp.jl
	julia --project=.  test/graphical_tests/accuracy_Vnk_10000_cmp.jl

test/graphical_tests/saves_for_graphical_tests/accuracy_Pkn_10000.jld: test/graphical_tests/accuracy_Pkn_10000_cmp.jl
	julia --project=.  test/graphical_tests/accuracy_Pkn_10000_cmp.jl

test/graphical_tests/saves_for_graphical_tests/accuracy_Cnk_10000.jld: test/graphical_tests/accuracy_Cnk_10000_cmp.jl
	julia --project=.  test/graphical_tests/accuracy_Cnk_10000_cmp.jl

### Graph for graphical tests 10000

test/graphical_tests/figures_graphical_tests/accuracy_PknCnkVnk_10000.pdf: test/graphical_tests/saves_for_graphical_tests/accuracy_Vnk_10000.jld test/graphical_tests/saves_for_graphical_tests/accuracy_Pkn_10000.jld test/graphical_tests/saves_for_graphical_tests/accuracy_Cnk_10000.jld test/graphical_tests/accuracy_PknCnkVnk_10000_plt.jl
	julia --project=.  test/graphical_tests/accuracy_PknCnkVnk_10000_plt.jl

### Save for Cnk graphical tests sigma

test/graphical_tests/saves_for_graphical_tests/accuracy_Cnk_sigma.jld: test/graphical_tests/common_functions_for_tests.jl test/graphical_tests/accuracy_Cnk_sigma_cmp.jl
	julia --project=.  test/graphical_tests/accuracy_Cnk_sigma_cmp.jl

test/graphical_tests/saves_for_graphical_tests/accuracy_Cnk_rec_sigma.jld: test/graphical_tests/accuracy_Cnk_rec_sigma_cmp.jl
	julia --project=.  test/graphical_tests/accuracy_Cnk_rec_sigma_cmp.jl

### Plot for Cnk graphical tests sigma

test/graphical_tests/figures_graphical_tests/accuracy_Cnk_sigma.pdf: test/graphical_tests/saves_for_graphical_tests/accuracy_Cnk_sigma.jld test/graphical_tests/saves_for_graphical_tests/accuracy_Cnk_rec_sigma.jld test/graphical_tests/accuracy_Cnk_sigma_plt.jl
	julia --project=.  test/graphical_tests/accuracy_Cnk_sigma_plt.jl

### Plot dist Pkn using robust approximation

test/graphical_tests/figures_graphical_tests/dist_Pkn_1000.pdf: test/graphical_tests/dist_Pkn_NGG_1000_plt.jl
	julia --project=.  test/graphical_tests/dist_Pkn_NGG_1000_plt.jl


### Plot Pkn approx quality decrease

test/graphical_tests/figures_graphical_tests/Pkn_NGG_approx_quality.pdf: test/graphical_tests/Pkn_NGG_approx_quality.jl
	julia --project=.  test/graphical_tests/Pkn_NGG_approx_quality.jl


### Investigate n limit

test/graphical_tests/saves_for_graphical_tests/accuracy_Vn1.jld: test/graphical_tests/accuracy_Vn1_1000_cmp.jl test/graphical_tests/common_functions_for_tests.jl
	export JULIA_NUM_THREADS=7; julia --project=.  test/graphical_tests/accuracy_Vn1_1000_cmp.jl

test/graphical_tests/figures_graphical_tests/accuracy_Vn1.pdf: test/graphical_tests/saves_for_graphical_tests/accuracy_Vn1.jld test/graphical_tests/accuracy_Vn1_1000_cmp.jl test/graphical_tests/accuracy_Vn1_1000_plt.jl
	julia --project=.  test/graphical_tests/accuracy_Vn1_1000_plt.jl

README.md: README.jmd
	julia --project=.  -e 'using Weave; weave("README.jmd", out_path=:pwd)'


### Run benchmark

benchmark_results.md: run_benchmarks.jl benchmark/benchmarks.jl
	julia --project=.  --project=.  run_benchmarks.jl

### Shortcut

graphical_tests: test/graphical_tests/figures_graphical_tests/accuracy_PknCnkVnk_1000.pdf test/graphical_tests/figures_graphical_tests/accuracy_PknCnkVnk_10000.pdf test/graphical_tests/figures_graphical_tests/accuracy_Cnk_sigma.pdf test/graphical_tests/figures_graphical_tests/dist_Pkn_1000.pdf test/graphical_tests/figures_graphical_tests/Pkn_NGG_approx_quality.pdf
.PHONY: graphical_tests

README: README.md
.PHONY: README

all: graphical_tests README
