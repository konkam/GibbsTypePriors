# Benchmark Report for */home/gkonkamking/Work/GibbsTypePriors*

## Job Properties
* Time of benchmark: 7 Apr 2020 - 19:38
* Package commit: fc013f
* Julia commit: b8e9a9
* Julia command flags: None
* Environment variables: None

## Results
Below is a table of this job's results, obtained by running the benchmarks.
The values listed in the `ID` column have the structure `[parent_group, child_group, ..., key]`, and can be used to
index into the BaseBenchmarks suite to retrieve the corresponding benchmarks.
The percentages accompanying time and memory values in the below table are noise tolerances. The "true"
time/memory value for a given benchmark is expected to fall within this percentage of the reported value.
An empty cell means that the value was zero.

| ID                                                 | time            | GC time | memory          | allocations |
|----------------------------------------------------|----------------:|--------:|----------------:|------------:|
| `["trigonometry", "circular", "(\"cos\", 0.0)"]`   |  30.000 ns (5%) |         |                 |             |
| `["trigonometry", "circular", "(\"cos\", π)"]`     |  26.000 ns (5%) |         |                 |             |
| `["trigonometry", "circular", "(\"sin\", 0.0)"]`   |  28.000 ns (5%) |         |                 |             |
| `["trigonometry", "circular", "(\"sin\", π)"]`     |  26.000 ns (5%) |         |                 |             |
| `["trigonometry", "circular", "(\"tan\", 0.0)"]`   |  28.000 ns (5%) |         |                 |             |
| `["trigonometry", "circular", "(\"tan\", π)"]`     |  26.000 ns (5%) |         |                 |             |
| `["trigonometry", "hyperbolic", "(\"cos\", 0.0)"]` |  30.000 ns (5%) |         |                 |             |
| `["trigonometry", "hyperbolic", "(\"cos\", π)"]`   |  26.000 ns (5%) |         |                 |             |
| `["trigonometry", "hyperbolic", "(\"sin\", 0.0)"]` |  29.000 ns (5%) |         |                 |             |
| `["trigonometry", "hyperbolic", "(\"sin\", π)"]`   |  26.000 ns (5%) |         |                 |             |
| `["trigonometry", "hyperbolic", "(\"tan\", 0.0)"]` |  31.000 ns (5%) |         |                 |             |
| `["trigonometry", "hyperbolic", "(\"tan\", π)"]`   |  26.000 ns (5%) |         |                 |             |
| `["utf8", "join"]`                                 | 117.034 ms (5%) |         | 156.27 MiB (1%) |          20 |
| `["utf8", "replace"]`                              | 115.183 μs (5%) |         |  12.06 KiB (1%) |           6 |

## Benchmark Group List
Here's a list of all the benchmark groups executed by this job:

- `["trigonometry", "circular"]`
- `["trigonometry", "hyperbolic"]`
- `["utf8"]`

## Julia versioninfo
```
Julia Version 1.4.0
Commit b8e9a9ecc6 (2020-03-21 16:36 UTC)
Platform Info:
  OS: Linux (x86_64-pc-linux-gnu)
      Ubuntu 18.04.4 LTS
  uname: Linux 5.3.0-45-generic #37~18.04.1-Ubuntu SMP Fri Mar 27 15:58:10 UTC 2020 x86_64 x86_64
  CPU: Intel(R) Core(TM) i7-8665U CPU @ 1.90GHz: 
              speed         user         nice          sys         idle          irq
       #1  4100 MHz       7399 s          3 s       1484 s     200083 s          0 s
       #2  3951 MHz       7191 s        135 s       1499 s     199902 s          0 s
       #3  4185 MHz       7240 s          9 s       1451 s     200474 s          0 s
       #4  4196 MHz       7858 s          0 s       1507 s     199648 s          0 s
       #5  4168 MHz       8030 s          2 s       1341 s     199622 s          0 s
       #6  3957 MHz       7096 s          2 s       1387 s     200636 s          0 s
       #7  4077 MHz       9255 s          2 s       1531 s     197954 s          0 s
       #8  3990 MHz       8329 s          0 s       1818 s     198876 s          0 s
       
  Memory: 31.242061614990234 GB (22306.6953125 MB free)
  Uptime: 2097.0 sec
  Load Avg:  1.07763671875  0.634765625  0.416015625
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-8.0.1 (ORCJIT, skylake)
```