# Sequential checks on the main series

2,266 of 2,266 daily random.org files; 19,008,341,000 bits, exactly those the experiments used. Each file is tested separately and the counts summed.

## Proportions under the null hypothesis of independent fair bits (expected 0.5)

| statistic | pairs or bits | proportion | z | two-sided p |
|---|---|---|---|---|
| ones | 19,008,341,000 | 0.5000001 | +0.03 | 0.973 |
| changes between successive bits (runs) | 19,008,338,734 | 0.4999990 | -0.26 | 0.792 |
| agreement at lag 2 | 19,008,336,468 | 0.4999930 | -1.94 | 0.053 |
| agreement at lag 3 | 19,008,334,202 | 0.4999991 | -0.24 | 0.814 |
| agreement at lag 4 | 19,008,331,936 | 0.4999962 | -1.04 | 0.299 |
| agreement at lag 5 | 19,008,329,670 | 0.5000007 | +0.18 | 0.858 |
| agreement at lag 6 | 19,008,327,404 | 0.4999964 | -0.98 | 0.325 |
| agreement at lag 7 | 19,008,325,138 | 0.4999962 | -1.06 | 0.289 |
| agreement at lag 8 | 19,008,322,872 | 0.5000022 | +0.62 | 0.537 |
| agreement at lag 9 | 19,008,320,606 | 0.4999931 | -1.92 | 0.055 |
| agreement at lag 10 | 19,008,318,340 | 0.4999986 | -0.39 | 0.699 |
| agreement at lag 11 | 19,008,316,074 | 0.5000006 | +0.18 | 0.859 |
| agreement at lag 12 | 19,008,313,808 | 0.5000032 | +0.87 | 0.384 |
| agreement at lag 13 | 19,008,311,542 | 0.5000035 | +0.97 | 0.330 |
| agreement at lag 14 | 19,008,309,276 | 0.4999959 | -1.14 | 0.253 |
| agreement at lag 15 | 19,008,307,010 | 0.5000046 | +1.28 | 0.202 |
| agreement at lag 16 | 19,008,304,744 | 0.4999951 | -1.35 | 0.178 |
| agreement at lag 250 | 19,007,774,500 | 0.5000029 | +0.80 | 0.424 |
| agreement at lag 500 | 19,007,208,000 | 0.4999961 | -1.08 | 0.282 |

19 tests; Bonferroni threshold for p = .05 is 0.0026. Lag 250 compares the same trial in the morning and afternoon halves of an experiment; lag 500 the same trial in consecutive experiments.

## Run lengths

9,504,148,928 complete runs. Chi-square against P(length = k) = 2^-k, lengths 1 to 24 and 25 or more: 22.6 on 24 df, p = 0.541.

| length | observed | expected | observed / expected |
|---|---|---|---|
| 1 | 4,751,999,202 | 4,752,074,464 | 1.0000 |
| 2 | 2,376,097,407 | 2,376,037,232 | 1.0000 |
| 3 | 1,188,018,020 | 1,188,018,616 | 1.0000 |
| 4 | 594,050,497 | 594,009,308 | 1.0001 |
| 5 | 297,005,533 | 297,004,654 | 1.0000 |
| 10 | 9,285,510 | 9,281,395 | 1.0004 |
| 15 | 290,338 | 290,044 | 1.0010 |
| 20 | 9,092 | 9,064 | 1.0031 |
| 25+ | 574 | 566 | 1.0133 |

## Day to day

If the days differed only by chance, each column of per-file z scores would have mean 0, standard deviation 1 and mean square 1.

| per-file z score | mean | SD | mean of z squared | excess as z | chi-square p |
|---|---|---|---|---|---|
| ones | +0.001 | 0.988 | 0.976 | -0.8 | 0.787 |
| runs (lag 1) | -0.006 | 1.179 | 1.389 | +13.1 | < .001 |
| agreement at lag 2 | -0.041 | 1.056 | 1.116 | +3.9 | < .001 |
| agreement at lag 3 | -0.005 | 1.001 | 1.001 | +0.0 | 0.482 |
| agreement at lag 4 | -0.022 | 1.009 | 1.018 | +0.6 | 0.273 |
| agreement at lag 5 | +0.004 | 0.980 | 0.960 | -1.4 | 0.914 |
| agreement at lag 6 | -0.021 | 1.021 | 1.043 | +1.4 | 0.075 |
| agreement at lag 7 | -0.022 | 1.007 | 1.014 | +0.5 | 0.319 |
| agreement at lag 8 | +0.013 | 0.995 | 0.990 | -0.3 | 0.624 |
| agreement at lag 9 | -0.040 | 1.019 | 1.040 | +1.4 | 0.089 |
| agreement at lag 10 | -0.008 | 1.002 | 1.003 | +0.1 | 0.455 |
| agreement at lag 11 | +0.004 | 1.002 | 1.003 | +0.1 | 0.450 |
| agreement at lag 12 | +0.018 | 1.020 | 1.040 | +1.4 | 0.088 |
| agreement at lag 13 | +0.020 | 0.988 | 0.977 | -0.8 | 0.778 |
| agreement at lag 14 | -0.024 | 0.976 | 0.952 | -1.6 | 0.948 |
| agreement at lag 15 | +0.027 | 1.014 | 1.029 | +1.0 | 0.162 |
| agreement at lag 16 | -0.028 | 0.998 | 0.996 | -0.1 | 0.550 |
| agreement at lag 250 | +0.017 | 1.005 | 1.010 | +0.3 | 0.366 |
| agreement at lag 500 | -0.023 | 1.003 | 1.006 | +0.2 | 0.414 |

For the runs, a correlation r between successive bits shifts a file's z score by about -r x sqrt(8,388,499). The variance of the runs z scores beyond 1 (0.389) therefore implies that the correlation between successive bits varies from day to day with a standard deviation of about 0.000215. A correlation r changes the variance of a 250-trial hit count by a factor of about 1 + 2r, so day-to-day variation of this size changes it by about 0.04%.

Correlation between the runs z scores of consecutive days (in date order): -0.005 (standard error about 0.021 if days are independent).

### Runs by year

| year | files | mean z | SD of z | mean of z squared | implied SD of correlation |
|---|---|---|---|---|---|
| 2006 | 226 | -0.066 | 0.976 | 0.957 | 0.000000 |
| 2007 | 362 | -0.047 | 1.030 | 1.063 | 0.000085 |
| 2008 | 366 | -0.050 | 1.041 | 1.086 | 0.000100 |
| 2009 | 365 | -0.069 | 1.021 | 1.046 | 0.000070 |
| 2010 | 365 | +0.056 | 0.996 | 0.995 | 0.000000 |
| 2011 | 365 | +0.015 | 1.201 | 1.444 | 0.000230 |
| 2012 | 217 | +0.171 | 2.004 | 4.043 | 0.000599 |

### The ten days with the most extreme runs z scores

| day | runs z | ones z |
|---|---|---|
| 2012-01-04 | +25.33 | +1.40 |
| 2011-09-20 | +9.35 | +0.89 |
| 2011-06-27 | +4.44 | +0.86 |
| 2011-10-13 | +4.23 | -0.33 |
| 2008-02-07 | +3.63 | +0.87 |
| 2008-01-15 | +3.51 | +0.00 |
| 2007-07-13 | +3.51 | -0.03 |
| 2008-12-19 | +3.49 | -0.74 |
| 2009-01-22 | +3.22 | -1.56 |
| 2010-06-17 | -3.13 | +0.75 |
