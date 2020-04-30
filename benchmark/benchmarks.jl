using BenchmarkTools
using Random
using GibbsTypePriors

# Define a parent BenchmarkGroup to contain our suite
const SUITE = BenchmarkGroup()

# Add some child groups to our benchmark suite. The most relevant BenchmarkGroup constructor
# for this case is BenchmarkGroup(tags::Vector). These tags are useful for
# filtering benchmarks by topic, which we'll cover in a later section.
SUITE["utf8"] = BenchmarkGroup(["string", "unicode"])

# Add some benchmarks to the "utf8" group
teststr = String(join(rand(MersenneTwister(1), 'a':'d', 10^4)))
SUITE["utf8"]["replace"] = @benchmarkable replace($teststr, "a" => "b")
SUITE["utf8"]["join"] = @benchmarkable join($teststr, $teststr)
SUITE["utf8"]["plots"] = BenchmarkGroup()

# # Add some benchmarks to the "trig" group
# SUITE["trigonometry"] = BenchmarkGroup(["math", "triangles"])
# SUITE["trigonometry"]["circular"] = BenchmarkGroup()
# for f in (sin, cos, tan)
#     for x in (0.0, pi)
#         SUITE["trigonometry"]["circular"][string(f), x] = @benchmarkable ($f)($x)
#     end
# end
#
# SUITE["trigonometry"]["hyperbolic"] = BenchmarkGroup()
# for f in (sin, cos, tan)
#     for x in (0.0, pi)
#         SUITE["trigonometry"]["hyperbolic"][string(f), x] = @benchmarkable ($f)($x)
#     end
# end

SUITE["Cnk"] = BenchmarkGroup()
SUITE["Cnk"]["direct"] = @benchmarkable GibbsTypePriors.Cnk(6, 5, 0.5)
SUITE["Cnk"]["recursive"] = @benchmarkable GibbsTypePriors.Cnk_rec(6, 5, 0.5)
SUITE["Cnk"]["recursive_arb"] = @benchmarkable GibbsTypePriors.Cnk_rec(6, 5, GibbsTypePriors.RR(0.5))
SUITE["Cnk"]["recursive_arb_long"] = @benchmarkable GibbsTypePriors.Cnk_rec(60, 5, GibbsTypePriors.RR(0.5))

SUITE["Vnk"] = BenchmarkGroup()
SUITE["Vnk"]["direct"] = @benchmarkable GibbsTypePriors.Vnk_NGG(950,50, 0.5, 0.2)
