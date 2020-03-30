
### Saves for graphical tests 1000

test/graphical_tests/saves_for_graphical_tests/accuracy_Vnk_1000.jld: test/graphical_tests/accuracy_Vnk_1000_cmp.jl
	julia test/graphical_tests/accuracy_Vnk_1000_cmp.jl

test/graphical_tests/saves_for_graphical_tests/accuracy_Pkn_1000.jld: test/graphical_tests/accuracy_Pkn_1000_cmp.jl
	julia test/graphical_tests/accuracy_Pkn_1000_cmp.jl

test/graphical_tests/saves_for_graphical_tests/accuracy_Cnk_1000.jld: test/graphical_tests/accuracy_Cnk_1000_cmp.jl
	julia test/graphical_tests/accuracy_Cnk_1000_cmp.jl

test/graphical_tests/figures_graphical_tests/accuracy_PknCnkVnk_1000.pdf: test/graphical_tests/saves_for_graphical_tests/accuracy_Vnk_1000.jld test/graphical_tests/saves_for_graphical_tests/accuracy_Pkn_1000.jld test/graphical_tests/saves_for_graphical_tests/accuracy_Cnk_1000.jld test/graphical_tests/accuracy_PknCnkVnk_1000_plt.jl
	julia test/graphical_tests/accuracy_PknCnkVnk_1000_plt.jl

test/graphical_tests/saves_for_graphical_tests/accuracy_Vnk_10000.jld: test/graphical_tests/accuracy_Vnk_10000_cmp.jl
	julia test/graphical_tests/accuracy_Vnk_10000_cmp.jl

test/graphical_tests/saves_for_graphical_tests/accuracy_Pkn_10000.jld: test/graphical_tests/accuracy_Pkn_10000_cmp.jl
	julia test/graphical_tests/accuracy_Pkn_10000_cmp.jl

test/graphical_tests/saves_for_graphical_tests/accuracy_Cnk_10000.jld: test/graphical_tests/accuracy_Cnk_10000_cmp.jl
	julia test/graphical_tests/accuracy_Cnk_10000_cmp.jl

test/graphical_tests/figures_graphical_tests/accuracy_PknCnkVnk_10000.pdf: test/graphical_tests/saves_for_graphical_tests/accuracy_Vnk_10000.jld test/graphical_tests/saves_for_graphical_tests/accuracy_Pkn_10000.jld test/graphical_tests/saves_for_graphical_tests/accuracy_Cnk_10000.jld test/graphical_tests/accuracy_PknCnkVnk_10000_plt.jl
	julia test/graphical_tests/accuracy_PknCnkVnk_10000_plt.jl

graphical_tests: test/graphical_tests/figures_graphical_tests/accuracy_PknCnkVnk_1000.pdf test/graphical_tests/figures_graphical_tests/accuracy_PknCnkVnk_10000.pdf
README.md: README.jmd
	julia -e 'using Weave; weave("README.jmd", out_path=:pwd)'
README: README.md