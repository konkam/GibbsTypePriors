using PkgBenchmark
import GibbsTypePriors

results = benchmarkpkg("/home/gkonkamking/Work/GibbsTypePriors", retune = true)
PkgBenchmark.export_markdown("benchmark_results.md", results)
