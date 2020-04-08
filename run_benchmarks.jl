using PkgBenchmark
import GibbsTypePriors

results = benchmarkpkg("/home/gkonkamking/Work/GibbsTypePriors", retune = true)
PkgBenchmark.export_markdown("benchmark_results.md", results)

judge_res = judge("GibbsTypePriors",
      "7ec7963270ffe82d8958166ebb70c0facadb1233",
      "0115344fe322cdfbe4c83af5da982399d281979b")

PkgBenchmark.export_markdown("benchmark_cmp.md", judge_res)
