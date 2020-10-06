# Benchmark Report for */home/gkonkamking/Work/GibbsTypePriors*

## Job Properties
* Time of benchmark: 6 Oct 2020 - 11:54
* Package commit: dirty
* Julia commit: 539f3c
* Julia command flags: None
* Environment variables: None

## Results
Below is a table of this job's results, obtained by running the benchmarks.
The values listed in the `ID` column have the structure `[parent_group, child_group, ..., key]`, and can be used to
index into the BaseBenchmarks suite to retrieve the corresponding benchmarks.
The percentages accompanying time and memory values in the below table are noise tolerances. The "true"
time/memory value for a given benchmark is expected to fall within this percentage of the reported value.
An empty cell means that the value was zero.

| ID                              | time            | GC time   | memory          | allocations |
|---------------------------------|----------------:|----------:|----------------:|------------:|
| `["Cnk", "direct"]`             |   4.737 μs (5%) |           |   4.03 KiB (1%) |          46 |
| `["Cnk", "recursive"]`          |  24.243 ns (5%) |           |                 |             |
| `["Cnk", "recursive_arb"]`      |  10.616 μs (5%) |           |   6.53 KiB (1%) |         115 |
| `["Cnk", "recursive_arb_long"]` | 353.090 μs (5%) |           | 196.88 KiB (1%) |        3535 |
| `["Vnk", "direct"]`             |    4.698 s (5%) | 17.577 ms |   3.51 GiB (1%) |     6216312 |
| `["utf8", "join"]`              | 122.558 ms (5%) |           | 156.27 MiB (1%) |          20 |
| `["utf8", "replace"]`           | 115.707 μs (5%) |           |  12.00 KiB (1%) |           4 |

## Benchmark Group List
Here's a list of all the benchmark groups executed by this job:

- `["Cnk"]`
- `["Vnk"]`
- `["utf8"]`

## Julia versioninfo
```
Julia Version 1.5.2
Commit 539f3ce943 (2020-09-23 23:17 UTC)
Platform Info:
  OS: Linux (x86_64-pc-linux-gnu)
      Ubuntu 18.04.5 LTS
  uname: Linux 5.4.0-48-generic #52~18.04.1-Ubuntu SMP Thu Sep 10 12:50:22 UTC 2020 x86_64 x86_64
  CPU: Intel(R) Core(TM) i7-8665U CPU @ 1.90GHz: 
              speed         user         nice          sys         idle          irq
       #1  3584 MHz    1226502 s        751 s     328614 s   10153271 s          0 s
       #2  3521 MHz    1109123 s       1584 s     304850 s   10296110 s          0 s
       #3  3561 MHz    1168216 s        966 s     312550 s   10228821 s          0 s
       #4  3528 MHz    1239989 s       1538 s     432388 s    9912799 s          0 s
       #5  3604 MHz    1251008 s        570 s     308406 s   10162030 s          0 s
       #6  3547 MHz    1257546 s       1267 s     344702 s   10101698 s          0 s
       #7  3583 MHz    1238304 s        754 s     335781 s   10115662 s          0 s
       #8  3553 MHz    1283264 s        505 s     327022 s   10107916 s          0 s
       
  Memory: 31.17916488647461 GB (2212.984375 MB free)
  Uptime: 169482.0 sec
  Load Avg:  1.59765625  1.470703125  1.46044921875
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-9.0.1 (ORCJIT, skylake)
```