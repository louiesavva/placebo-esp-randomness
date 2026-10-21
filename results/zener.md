# random.org Zener-card series

473,594 experiments of 625 trials, hit probability 1/5 (two random symbols out of five match). CR = (X - 125)/10.

## Moments and fit

| statistic | observed | binomial | z |
|---|---|---|---|
| mean | 125.0111 | 125 | +0.76 |
| variance | 99.6330 | 100 | -1.79 |
| skewness | 0.0591 | 0.0600 | -0.24 |
| excess kurtosis | 0.01783 | 0.00040 | +2.45 |

Chi-square goodness of fit to Binomial(625, 1/5): 75.7 on 82 df, p = 0.674 (cells with expected count below 5 pooled into the two tails).

## Experiments with |CR| >= t

| t | needs abs(X - 125) >= | observed | exact | normal | normal, corrected | obs/exact | obs/normal | exact/normal | (obs - exact)/SE |
|---|---|---|---|---|---|---|---|---|---|
| 2.0 | 20 | 24,018 | 24,144.9 | 21,548.7 | 24,236.7 | 0.995 | 1.115 | 1.120 | -0.84 |
| 2.1 | 21 | 18,960 | 19,034.2 | 16,921.0 | 19,116.4 | 0.996 | 1.121 | 1.125 | -0.55 |
| 2.2 | 22 | 14,849 | 14,873.1 | 13,169.2 | 14,944.4 | 0.998 | 1.128 | 1.129 | -0.20 |
| 2.3 | 23 | 11,526 | 11,519.1 | 10,157.7 | 11,578.9 | 1.001 | 1.135 | 1.134 | +0.07 |
| 2.4 | 24 | 8,831 | 8,842.6 | 7,764.6 | 8,891.0 | 0.999 | 1.137 | 1.139 | -0.12 |
| 2.5 | 25 | 6,741 | 6,728.0 | 5,881.7 | 6,765.6 | 1.002 | 1.146 | 1.144 | +0.16 |
| 2.6 | 26 | 5,155 | 5,074.0 | 4,415.0 | 5,101.7 | 1.016 | 1.168 | 1.149 | +1.14 |
| 2.7 | 27 | 3,893 | 3,792.9 | 3,283.9 | 3,812.0 | 1.026 | 1.185 | 1.155 | +1.63 |
| 2.8 | 28 | 2,884 | 2,810.4 | 2,420.2 | 2,822.4 | 1.026 | 1.192 | 1.161 | +1.39 |
| 2.9 | 29 | 2,123 | 2,064.3 | 1,767.3 | 2,070.5 | 1.028 | 1.201 | 1.168 | +1.30 |
| 3.0 | 30 | 1,536 | 1,503.1 | 1,278.6 | 1,505.0 | 1.022 | 1.201 | 1.176 | +0.85 |
| 3.5 | 35 | 283 | 271.0 | 220.3 | 265.5 | 1.044 | 1.284 | 1.230 | +0.73 |
| 4.0 | 40 | 30 | 40.0 | 30.0 | 37.0 | 0.750 | 1.000 | 1.333 | -1.58 |
| 4.5 | 45 | 4 | 4.9 | 3.2 | 4.1 | 0.813 | 1.243 | 1.529 | -0.42 |

The exact count exceeds the normal prediction at every threshold, because each threshold coincides with a possible score.

## Upper and lower tails separately

| t | upper: observed | exact | normal | lower: observed | exact | normal |
|---|---|---|---|---|---|---|
| 2.0 | 12,774 | 12,865.7 | 10,774.3 | 11,244 | 11,279.2 | 10,774.3 |
| 2.5 | 3,835 | 3,834.2 | 2,940.9 | 2,906 | 2,893.8 | 2,940.9 |
| 3.0 | 953 | 938.4 | 639.3 | 583 | 564.7 | 639.3 |
| 3.5 | 201 | 188.6 | 110.2 | 82 | 82.4 | 110.2 |

With p = 1/5 the binomial is skewed to the right. Together with the lattice effect, this makes the upper tail (more hits than expected, the direction an ESP test looks for) much heavier than the normal curve, while from t = 2.5 upwards the lower tail is lighter.

## True size of nominal one-sided z tests, 625 trials

| nominal | smallest qualifying X | exact size | observed proportion (95% interval) |
|---|---|---|---|
| z >= 1.645 (5%) | 142 | 5.11% | 5.08% (5.01% to 5.14%) |
| z >= 1.96 (2.5%) | 145 | 2.72% | 2.70% (2.65% to 2.74%) |
| z >= 2.326 (1%) | 149 | 1.05% | 1.04% (1.01% to 1.07%) |

