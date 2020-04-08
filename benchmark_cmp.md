# Benchmark Report for *GibbsTypePriors*

## Job Properties
* Time of benchmarks:
    - Target: 8 Apr 2020 - 17:49
    - Baseline: 8 Apr 2020 - 17:49
* Package commits:
    - Target: 2611aa
    - Baseline: 011534
* Julia commits:
    - Target: b8e9a9
    - Baseline: b8e9a9
* Julia command flags:
    - Target: None
    - Baseline: None
* Environment variables:
    - Target: None
    - Baseline: None

## Results
A ratio greater than `1.0` denotes a possible regression (marked with :x:), while a ratio less
than `1.0` denotes a possible improvement (marked with :white_check_mark:). Only significant results - results
that indicate possible regressions or improvements - are shown below (thus, an empty table means that all
benchmark results remained invariant between builds).

| ID                                                 | time ratio                   | memory ratio |
|----------------------------------------------------|------------------------------|--------------|
| `["Cnk", "direct"]`                                | 0.93 (5%) :white_check_mark: |   0.99 (1%)  |
| `["Cnk", "recursive"]`                             |               69.08 (5%) :x: | Inf (1%) :x: |
| `["Vnk", "direct"]`                                |                1.09 (5%) :x: |   0.99 (1%)  |
| `["trigonometry", "circular", "(\"tan\", 0.0)"]`   |                1.07 (5%) :x: |   1.00 (1%)  |
| `["trigonometry", "hyperbolic", "(\"sin\", π)"]`   |                1.12 (5%) :x: |   1.00 (1%)  |
| `["trigonometry", "hyperbolic", "(\"tan\", 0.0)"]` |                1.07 (5%) :x: |   1.00 (1%)  |
| `["trigonometry", "hyperbolic", "(\"tan\", π)"]`   |                1.07 (5%) :x: |   1.00 (1%)  |

## Benchmark Group List
Here's a list of all the benchmark groups executed by this job:

- `["Cnk"]`
- `["Vnk"]`
- `["trigonometry", "circular"]`
- `["trigonometry", "hyperbolic"]`
- `["utf8"]`

## Julia versioninfo

### Target
```
Julia Version 1.4.0
Commit b8e9a9ecc6 (2020-03-21 16:36 UTC)
Platform Info:
  OS: Linux (x86_64-pc-linux-gnu)
      Ubuntu 18.04.4 LTS
  uname: Linux 5.3.0-45-generic #37~18.04.1-Ubuntu SMP Fri Mar 27 15:58:10 UTC 2020 x86_64 x86_64
  CPU: Intel(R) Core(TM) i7-8665U CPU @ 1.90GHz: 
              speed         user         nice          sys         idle          irq
       #1  2631 MHz     227073 s        160 s      64683 s    3365657 s          0 s
       #2  2617 MHz     249195 s        327 s      61730 s    3350697 s          0 s
       #3  2642 MHz     234262 s        147 s      63629 s    3351572 s          0 s
       #4  2666 MHz     237267 s         45 s      69076 s    3345254 s          0 s
       #5  2647 MHz     232759 s        145 s      62680 s    3364947 s          0 s
       #6  2631 MHz     226981 s        365 s      62537 s    3372626 s          0 s
       #7  2693 MHz     230682 s        181 s      63391 s    3361916 s          0 s
       #8  2677 MHz     230548 s        142 s      64552 s    3366312 s          0 s
       
  Memory: 31.242061614990234 GB (14024.94921875 MB free)
  Uptime: 81927.0 sec
  Load Avg:  1.38818359375  1.16650390625  0.9462890625
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-8.0.1 (ORCJIT, skylake)
```

### Baseline
```
Julia Version 1.4.0
Commit b8e9a9ecc6 (2020-03-21 16:36 UTC)
Platform Info:
  OS: Linux (x86_64-pc-linux-gnu)
      Ubuntu 18.04.4 LTS
  uname: Linux 5.3.0-45-generic #37~18.04.1-Ubuntu SMP Fri Mar 27 15:58:10 UTC 2020 x86_64 x86_64
  CPU: Intel(R) Core(TM) i7-8665U CPU @ 1.90GHz: 
              speed         user         nice          sys         idle          irq
       #1  3936 MHz     227188 s        160 s      64731 s    3369251 s          0 s
       #2  3979 MHz     249453 s        327 s      61782 s    3354147 s          0 s
       #3  3948 MHz     234611 s        147 s      63668 s    3354946 s          0 s
       #4  3927 MHz     237380 s         45 s      69124 s    3348848 s          0 s
       #5  3857 MHz     232857 s        145 s      62734 s    3368545 s          0 s
       #6  3619 MHz     230084 s        365 s      62802 s    3373038 s          0 s
       #7  3707 MHz     230782 s        181 s      63454 s    3365508 s          0 s
       #8  3897 MHz     230656 s        142 s      64598 s    3369916 s          0 s
       
  Memory: 31.242061614990234 GB (14030.8125 MB free)
  Uptime: 81965.0 sec
  Load Avg:  1.19775390625  1.1435546875  0.9501953125
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-8.0.1 (ORCJIT, skylake)
```