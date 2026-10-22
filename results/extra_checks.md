# Extra checks on the random.org main series

38,016,682 experiments.

## 1. Fifteen thresholds

| t | needs abs(M - A) >= | observed | exact | obs/exact | z |
|---|---|---|---|---|---|
| 2.0 | 23 | 1,675,231 | 1,675,146.9 | 1.000 | +0.07 |
| 2.1 | 24 | 1,347,052 | 1,347,943.7 | 0.999 | -0.78 |
| 2.2 | 25 | 1,077,163 | 1,076,866.6 | 1.000 | +0.29 |
| 2.3 | 26 | 854,353 | 854,090.6 | 1.000 | +0.29 |
| 2.4 | 27 | 672,178 | 672,479.6 | 1.000 | -0.37 |
| 2.5 | 28 | 525,349 | 525,617.4 | 0.999 | -0.37 |
| 2.6 | 30 | 314,663 | 314,071.8 | 1.002 | +1.06 |
| 2.7 | 31 | 240,215 | 240,085.2 | 1.001 | +0.27 |
| 2.8 | 32 | 182,396 | 182,159.7 | 1.001 | +0.55 |
| 2.9 | 33 | 137,509 | 137,175.0 | 1.002 | +0.90 |
| 3.0 | 34 | 102,729 | 102,522.5 | 1.002 | +0.65 |
| 3.5 | 40 | 15,542 | 15,222.5 | 1.021 | +2.59 |
| 4.0 | 45 | 2,557 | 2,510.2 | 1.019 | +0.93 |
| 4.5 | 51 | 211 | 222.5 | 0.948 | -0.77 |
| 5.0 | 56 | 22 | 23.7 | 0.929 | -0.35 |

Largest |z| = 2.59. In 20,000 simulated binomial series, a largest |z| at least this big occurred in 9.0% of them.

## 2. Kurtosis and central peak

Standard error of excess kurtosis = sqrt(24/N) = 0.000795.

| score | variance (expected) | excess kurtosis | binomial value | z |
|---|---|---|---|---|
| difference M - A | 124.9633 (125) | -0.00215 | -0.00400 | +2.33 |
| total M + A | 124.9865 (125) | -0.00444 | -0.00400 | -0.55 |
| morning | 62.4765 (62.5) | -0.00804 | -0.00800 | -0.05 |
| afternoon | 62.4984 (62.5) | -0.00750 | -0.00800 | +0.63 |

If the standard error varied from experiment to experiment, with relative variance of the variance v, the excess kurtosis of both scores would rise by about 3v.
From the total score, 3v <= 0.00112 (upper 95% limit), so v <= 0.000372: the variance varies by at most 1.9% and the standard error by at most about 1.0%.

| score | cells | observed | expected | z |
|---|---|---|---|---|
| difference M - A | abs <= 0 | 1,358,765 | 1,355,851 | +2.55 |
| difference M - A | abs <= 1 | 4,062,383 | 4,056,751 | +2.96 |
| difference M - A | abs <= 2 | 6,732,108 | 6,725,497 | +2.81 |
| total M + A | abs <= 0 | 1,355,049 | 1,355,851 | -0.70 |
| total M + A | abs <= 1 | 4,056,475 | 4,056,751 | -0.14 |
| total M + A | abs <= 2 | 6,726,295 | 6,725,497 | +0.34 |

Nominal two-sided p for the excess within one hit of M = A (z = 2.96): 0.0031. This window was chosen after looking at the data, among overlapping comparisons, so the nominal p overstates the evidence.

## 3. Upper and lower tails of the difference

| t | morning ahead (M - A >= L) | afternoon ahead (A - M >= L) | expected each |
|---|---|---|---|
| 3.0 | 51,263 | 51,466 | 51,261 |
| 3.5 | 7,639 | 7,903 | 7,611 |
| 4.0 | 1,244 | 1,313 | 1,255 |

Afternoon-ahead count at t = 3.5: z = 3.34, nominal two-sided p = 0.0008 (chosen after looking at the data).

## 4. random.org digits, +6 displacement

chi-square 78.5 on 44 df, p = 0.00106; Bonferroni threshold for 35 fits = 0.00143.
Mean 26.8057 (expected 26.8), variance 24.1046 (expected 24.1).
First and second halves of the series: p = 0.003 and 0.096.

## 5. Bit layout and runs

The bits of two daily random.org files survive in October 2012 workspaces (8,388,500 used from each). Taken 500 at a time in the order they were read, the first 250 as morning and the next 250 as afternoon, they reproduce the corresponding experiments of the main series (counted from 0). Runs and lag-1 correlation test for dependence between successive bits.

| day | experiments | reproduced exactly | proportion of ones | runs | expected | z | lag-1 correlation (SE) |
|---|---|---|---|---|---|---|---|
| 2012-06-17 | 37,194,609 to 37,211,385 | True | 0.49997 | 4,195,568 | 4,194,250 | +0.91 | -0.00031 (0.00035) |
| 2012-08-04 | 37,999,905 to 38,016,681 | True | 0.50025 | 4,195,136 | 4,194,250 | +0.61 | -0.00021 (0.00035) |

## 6. Without the six days with odd bit-level statistics

Days left out: 2012-01-04, 2011-09-20, 2008-04-19, 2011-06-27, 2011-10-13, 2007-04-29 (100,662 experiments, 0.26%).

| t | z, all days | z, without the six days | change |
|---|---|---|---|
| 2.0 | +0.07 | +0.18 | +0.11 |
| 2.1 | -0.78 | -0.72 | +0.06 |
| 2.2 | +0.29 | +0.33 | +0.04 |
| 2.3 | +0.29 | +0.27 | -0.01 |
| 2.4 | -0.37 | -0.39 | -0.02 |
| 2.5 | -0.37 | -0.41 | -0.04 |
| 2.6 | +1.06 | +1.08 | +0.02 |
| 2.7 | +0.27 | +0.27 | +0.00 |
| 2.8 | +0.55 | +0.54 | -0.01 |
| 2.9 | +0.90 | +0.89 | -0.02 |
| 3.0 | +0.65 | +0.62 | -0.02 |
| 3.5 | +2.59 | +2.54 | -0.05 |
| 4.0 | +0.93 | +0.91 | -0.03 |
| 4.5 | -0.77 | -0.73 | +0.04 |
| 5.0 | -0.35 | -0.33 | +0.01 |

Excess within one hit of M = A: z = 2.96 with all days, 2.92 without the six. Variance of M - A: 124.9633 with all days, 124.9657 without (binomial 125).

