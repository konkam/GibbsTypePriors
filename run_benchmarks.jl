using PkgBenchmark
import GibbsTypePriors

results = benchmarkpkg("/home/gkonkamking/Work/GibbsTypePriors", retune = true)
PkgBenchmark.export_markdown("benchmark_results.md", results)

judge_res = judge("GibbsTypePriors",
      "0510ea6edd45ad53c4bf8c7490239bdfb3fe0f72",
      "0115344fe322cdfbe4c83af5da982399d281979b")

PkgBenchmark.export_markdown("benchmark_cmp.md", judge_res)
