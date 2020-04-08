using PkgBenchmark
import GibbsTypePriors

results = benchmarkpkg("/home/gkonkamking/Work/GibbsTypePriors", retune = true)
PkgBenchmark.export_markdown("benchmark_results.md", results)

judge_res = judge("GibbsTypePriors",
      "2611aaf8",
      "0115344fe322cdfbe4c83af5da982399d281979b")
