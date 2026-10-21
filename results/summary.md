# Results

## 1. random.org main series (38,016,682 experiments of 250 morning + 250 afternoon trials)

Critical ratio = (morning hits - afternoon hits) / sqrt(125). Under chance, morning - afternoon + 250 ~ Binomial(500, 1/2) exactly.

| t | observed | exact binomial | normal | continuity-corrected | (obs-exact)/SE |
|---|---|---|---|---|---|
| 2.0 | 1,675,231 | 1,675,147 | 1,729,769 | 1,679,248 | +0.07 |
| 2.2 | 1,077,163 | 1,076,867 | 1,057,126 | 1,080,690 | +0.29 |
| 2.5 | 525,349 | 525,617 | 472,142 | 528,671 | -0.37 |
| 2.6 | 314,663 | 314,072 | 354,406 | 316,523 | +1.06 |
| 3.0 | 102,729 | 102,523 | 102,637 | 103,880 | +0.65 |
| 4.0 | 2,557 | 2,510 | 2,408 | 2,618 | +0.93 |
| 5.0 | 22 | 24 | 22 | 26 | -0.35 |

Nominal |CR| >= 1.96 test: requires |deviation| >= 22; exact size 0.0544.

| threshold | exact two-sided probability | observed / expected | 95% range if exact |
|---|---|---|---|
| 3.0 | 0.0027 | 1.002 | 0.99-1.01 |
| 4.0 | 6.6e-05 | 1.019 | 0.96-1.04 |
| 5.0 | 6.2e-07 | 0.929 | 0.63-1.44 |

## 2. Gilmore: does the standard error itself vary?

| quantity | observed | expected | test |
|---|---|---|---|
| mean of difference | -0.00113 | 0 | -0.63 |
| variance of difference | 124.9633 | 125 | -1.28 |
| excess kurtosis of difference | -0.00215 | -0.0040 (binomial); 0 (normal) | +2.33 vs binomial; -2.71 vs normal |
| excess kurtosis, morning | -0.00804 | -0.0080 | -0.05 |
| excess kurtosis, afternoon | -0.00750 | -0.0080 | +0.63 |
| experiments scoring exactly 0 | 1,358,765 | 1,355,851 | +2.55 |
| chi-square, whole distribution | 121.5 on 112 df |  | p = 0.254 |
| spread of daily variances (2266 days) | SD 1.364 | SD 1.363 | dispersion 2267 on 2266 df, p = 0.49 |
| spread of daily means | 2290 | 2266 df | p = 0.36 |
| corr(morning, afternoon) | +0.000093 | 0 | +0.57 |
| corr(morning^2, afternoon^2) | +0.000197 | 0 | +1.22 |

## 3. Every series against the exact binomial

| source | design | score | experiments | n | p | variance / expected | chi-square p | |CR|>=3: obs/exact | obs/normal |
|---|---|---|---|---|---|---|---|---|---|
| random.org (Example 1B) | morning-afternoon | difference | 38,016,682 | 500 | 0.5 | 0.9997 | 0.254 | 1.002 | 1.001 |
| random.org (Example 1B) | morning-afternoon | total | 38,016,682 | 500 | 0.5 | 0.9999 | 0.151 | 0.995 | 0.994 |
| random.org (Example 1B) | morning-afternoon | morning | 38,016,682 | 250 | 0.5 | 0.9996 | 0.032 | 1.000 | 1.066 |
| random.org (Example 1B) | morning-afternoon | afternoon | 38,016,682 | 250 | 0.5 | 1.0000 | 0.903 | 0.999 | 1.065 |
| R Mersenne Twister (Example 1A) | morning-afternoon | difference | 1,000,004 | 500 | 0.5 | 1.0007 | 0.699 | 0.972 | 0.971 |
| R Mersenne Twister (Example 1A) | morning-afternoon | total | 1,000,004 | 500 | 0.5 | 0.9998 | 0.948 | 1.022 | 1.020 |
| R Mersenne Twister (Example 1A) | morning-afternoon | morning | 1,000,004 | 250 | 0.5 | 1.0025 | 0.401 | 1.013 | 1.080 |
| R Mersenne Twister (Example 1A) | morning-afternoon | afternoon | 1,000,004 | 250 | 0.5 | 0.9980 | 0.422 | 0.989 | 1.054 |
| Entropy Key (Example 2) | digit matching | straight | 1,000,000 | 500 | 0.1 | 1.0004 | 0.176 | 0.997 | 0.859 |
| Entropy Key (Example 2) | digit matching | +1 displacement | 1,000,000 | 490 | 0.1 | 0.9990 | 0.636 | 0.974 | 1.223 |
| Entropy Key (Example 2) | digit matching | -1 displacement | 1,000,000 | 490 | 0.1 | 1.0017 | 0.006 | 1.021 | 1.282 |
| random.org (Example 3) | digit matching | straight | 1,108,074 | 268 | 0.1 | 0.9989 | 0.766 | 1.013 | 1.172 |
| random.org (Example 3) | digit matching | +1 displacement | 1,108,074 | 268 | 0.1 | 0.9973 | 0.138 | 0.976 | 1.129 |
| random.org (Example 3) | digit matching | -1 displacement | 1,108,074 | 268 | 0.1 | 0.9999 | 0.298 | 1.016 | 1.175 |
| random.org (Example 3) | digit matching | +2 displacement | 1,108,074 | 268 | 0.1 | 1.0002 | 0.465 | 1.029 | 1.190 |
| random.org (Example 3) | digit matching | -2 displacement | 1,108,074 | 268 | 0.1 | 1.0002 | 0.074 | 1.006 | 1.164 |
| random.org (Example 3) | digit matching | +3 displacement | 1,108,074 | 268 | 0.1 | 0.9990 | 0.737 | 1.018 | 1.178 |
| random.org (Example 3) | digit matching | -3 displacement | 1,108,074 | 268 | 0.1 | 1.0022 | 0.545 | 1.005 | 1.163 |
| random.org (Example 3) | digit matching | +4 displacement | 1,108,074 | 268 | 0.1 | 0.9981 | 0.818 | 0.997 | 1.154 |
| random.org (Example 3) | digit matching | -4 displacement | 1,108,074 | 268 | 0.1 | 1.0001 | 0.408 | 1.012 | 1.171 |
| random.org (Example 3) | digit matching | +5 displacement | 1,108,074 | 268 | 0.1 | 1.0006 | 0.699 | 0.984 | 1.139 |
| random.org (Example 3) | digit matching | -5 displacement | 1,108,074 | 268 | 0.1 | 1.0005 | 0.977 | 1.029 | 1.191 |
| random.org (Example 3) | digit matching | +6 displacement | 1,108,074 | 268 | 0.1 | 0.9994 | 0.001 | 1.005 | 1.163 |
| random.org (Example 3) | digit matching | -6 displacement | 1,108,074 | 268 | 0.1 | 1.0002 | 0.464 | 1.017 | 1.176 |
| random.org (Example 3) | digit matching | +7 displacement | 1,108,074 | 268 | 0.1 | 1.0010 | 0.491 | 0.983 | 1.137 |
| random.org (Example 3) | digit matching | -7 displacement | 1,108,074 | 268 | 0.1 | 0.9999 | 0.704 | 0.981 | 1.135 |
| random.org, 125 + 125 | morning-afternoon | total | 2,266,000 | 250 | 0.5 | 0.9992 | 0.956 | 0.998 | 1.064 |
| random.org, 125 + 125 | morning-afternoon | morning | 2,266,000 | 125 | 0.5 | 1.0005 | 0.494 | 1.003 | 0.827 |
| random.org, 125 + 125 | morning-afternoon | afternoon | 2,266,000 | 125 | 0.5 | 0.9975 | 0.090 | 0.999 | 0.824 |
| mcr.RData run (source not recorded) | morning-afternoon | difference | 1,000,000 | 500 | 0.5 | 1.0007 | 0.435 | 1.038 | 1.037 |
| mcr.RData run (source not recorded) | morning-afternoon | morning | 1,000,000 | 250 | 0.5 | 1.0009 | 0.732 | 1.020 | 1.088 |
| mcr.RData run (source not recorded) | morning-afternoon | afternoon | 1,000,000 | 250 | 0.5 | 0.9982 | 0.824 | 1.006 | 1.073 |
| aao1rng run (source not recorded) | morning-afternoon | difference | 255,732 | 500 | 0.5 | 1.0015 | 0.614 | 1.003 | 1.002 |
| aao1rng run (source not recorded) | morning-afternoon | morning | 255,732 | 250 | 0.5 | 0.9968 | 0.957 | 1.070 | 1.141 |
| aao1rng run (source not recorded) | morning-afternoon | afternoon | 255,732 | 250 | 0.5 | 1.0038 | 0.058 | 1.030 | 1.098 |

## 4. Digit matching (p = 1/10): two-sided counts against normal and exact

**Entropy Key, straight, 500 pairs** (1,000,000 experiments)

| t | observed | normal | exact | obs/normal | obs/exact |
|---|---|---|---|---|---|
| 2.0 | 43,689 | 45,500 | 43,686 | 0.96 | 1.000 |
| 2.5 | 13,872 | 12,419 | 13,776 | 1.12 | 1.007 |
| 3.0 | 2,320 | 2,700 | 2,326 | 0.86 | 0.997 |
| 3.5 | 532 | 465 | 529 | 1.14 | 1.005 |
| 4.0 | 111 | 63 | 109 | 1.75 | 1.017 |
| 4.5 | 11 | 7 | 12 | 1.62 | 0.951 |

Upper tail at 3.5 SD: observed 461; exact 455; normal 233; continuity-corrected 230.

**random.org, straight, 268 pairs** (1,108,074 experiments)

| t | observed | normal | exact | obs/normal | obs/exact |
|---|---|---|---|---|---|
| 2.0 | 45,912 | 50,418 | 46,144 | 0.91 | 0.995 |
| 2.5 | 11,564 | 13,762 | 11,675 | 0.84 | 0.990 |
| 3.0 | 3,505 | 2,992 | 3,461 | 1.17 | 1.013 |
| 3.5 | 867 | 516 | 870 | 1.68 | 0.997 |
| 4.0 | 118 | 70 | 123 | 1.68 | 0.958 |
| 4.5 | 26 | 8 | 29 | 3.45 | 0.886 |

Upper tail at 3.5 SD: observed 835; exact 832; normal 258; continuity-corrected 373.

## 5. True size of nominal one-sided tests

| design | n | p | z>=1.645 (5%) | z>=1.96 (2.5%) | z>=2.326 (1%) |
|---|---|---|---|---|---|
| One Zener run | 25 | 0.2 | 4.68% (X>=9) | 4.68% (X>=9) | 1.73% (X>=10) |
| Four Zener runs | 100 | 0.2 | 5.58% (X>=27) | 3.42% (X>=28) | 1.12% (X>=30) |
| Eight Zener runs | 200 | 0.2 | 4.94% (X>=50) | 2.36% (X>=52) | 1.03% (X>=54) |
| Four-choice, 40 trials | 40 | 0.25 | 5.44% (X>=15) | 2.62% (X>=16) | 1.16% (X>=17) |
| Four-choice, 100 trials | 100 | 0.25 | 4.46% (X>=33) | 2.76% (X>=34) | 0.94% (X>=36) |
| Binary, 100 trials | 100 | 0.5 | 4.43% (X>=59) | 2.84% (X>=60) | 1.05% (X>=62) |
| Binary, 500 trials | 500 | 0.5 | 4.89% (X>=269) | 2.72% (X>=272) | 0.88% (X>=277) |
| Binary, 1,000 trials | 1000 | 0.5 | 4.68% (X>=527) | 2.68% (X>=531) | 1.05% (X>=537) |
