# Placebo ESP experiments and the normal approximation

Data and code for two papers by Louie Savva:

- *Anomalous Significance Revisited: Nineteen Billion Placebo Trials and the Normal Approximation* (research paper, in preparation)
- *Too Many Miracles* (non-technical article, under review)

Between 2006 and 2013 the author ran tens of millions of placebo "ESP" experiments, in which random numbers stood in for both the target and the guess, so no psi was possible. An analysis in 2012 found too many extreme critical ratios and took this as evidence that real random sources misbehave in their tails. That conclusion was wrong. The excess is real when the counts are judged against the normal approximation to the binomial, and disappears when they are judged against the exact binomial. This repository holds every experiment's hit counts, the code that recovered them from the original R workspaces, and the code for every table in both papers.

## Contents

```
code/
  rdata_reader.py    reads R workspaces (.RData) in Python, without R
  extract.py         original .RData files  ->  data/*.npz
  analyse.py         data/*.npz  ->  results/ (every table except the Iversen simulation)
  extra_checks.py    follow-up checks on the main series  ->  results/extra_checks.md
  iversen_check.py   simulation of the Iversen et al. (1971) test battery
  figure1.py         Figure 1 (black and white): error of the normal approximations by threshold
  find_rdata.py      utility used to search the author's disks for the original files
data/                hit counts per experiment, unsigned 8-bit, numpy .npz (67 MB)
results/             output of the scripts above
```

## The data

| File | Series | Random source | Design | Experiments | Score(s) | Original workspace | 2012 analysis |
|---|---|---|---|---|---|---|---|
| `random_org_main_morning.npz`, `random_org_main_afternoon.npz` | Main series | random.org (atmospheric noise) | Spencer-Brown's morning-afternoon: 250 + 250 binary trials, p = 1/2 | 38,016,682 | hits in each half (of 250) | `rand1.RData` (checked against `Spencp.RData`) | Example 1B |
| `prng.npz` | Pseudo-random | Mersenne Twister in R | morning-afternoon, 250 + 250 | 1,000,004 | hits in each half | `prng.RData` | Example 1A |
| `entropy_key_digits.npz` | Hardware | Simtec Electronics Entropy Key (USB hardware generator using avalanche noise, read from `/dev/random`) | Spencer-Brown's digit matching: 500 pairs, p = 1/10; +1 and -1 displacements of 490 pairs | 1,000,000 | `straight`, `plus1`, `minus1` | `mcrx.RData` | Example 2 |
| `random_org_digits.npz` | random.org digits | random.org | digit matching, 268 pairs, p = 1/10, displacements +1 to +7 and -1 to -7 | 1,108,074 | `straight`, `plus1`..`plus7`, `minus1`..`minus7` | `randorgorigA.RData` | Example 3 |
| `random_org_125.npz` | Short series | random.org | morning-afternoon, 125 + 125 | 2,266,000 | hits in each half (of 125) | `randmorn.RData` | not used |
| `other_runs.npz` | Unlabelled runs | not recorded | morning-afternoon, 250 + 250 | 1,000,000 and 255,732 | `mcr_*`, `aao1rng_*` | `mcr.RData`, `aao1rng.RData` | not used |

`random_org_files.txt` lists the 2,266 daily random.org files used by the main series, 12 March 2006 to 4 August 2012, in order. Each file supplied 8,388,608 random bits, divided into 16,777 blocks of 500: one bit per trial (250 morning, 250 afternoon), serving as the hit indicator. This is equivalent to matching a random target with a random guess, since the match of two independent fair bits is itself a fair bit. The remaining 108 bits of each file were not used. Experiment *i* comes from file *i* // 16,777. The script that read the files was not saved with the analysis workspaces, so the order in which bits were taken from each byte is not recorded; under the null hypothesis it does not affect the distribution of hit counts.

Load any series with numpy:

```python
import numpy as np
m = np.load('data/random_org_main_morning.npz')['morning']      # hits out of 250
a = np.load('data/random_org_main_afternoon.npz')['afternoon']
x = m.astype(int) - a + 250                                      # exactly Binomial(500, 1/2) under chance
```

### Scores

For the morning-afternoon design, with *M* and *A* the morning and afternoon hits, Spencer-Brown's critical ratio is CR = (*M* - *A*) / sqrt(125). Because 250 - *A* has the same distribution as *A*, *X* = *M* - *A* + 250 is exactly Binomial(500, 1/2), the same as one experiment of 500 trials. The total *M* + *A* has the same distribution and is analysed too.

For digit matching, hits are Binomial(*n*, 1/10) with *n* = 500, 490 or 268, and CR = (*X* - *np*) / sqrt(*np*(1 - *p*)).

### Provenance

- **One copy of each series.** Every file in `data/` holds one complete, checked copy of its series. Nothing was excluded: the main series is exactly 2,266 files × 16,777 experiments = 38,016,682.
- **How the counts were recovered.** The original R workspaces stored critical ratios rather than hit counts. `extract.py` converts them back to whole-number counts and stops if any value does not convert exactly.
- **Where the main series comes from.** `rand1.RData`. An older copy of the same experiments, `Spencp.RData`, has 1,367 damaged values (0.004% of experiments), in three blocks of about a thousand experiments, probably caused by memory or disk errors. It is used only for the list of random.org files and as a cross-check (its undamaged values agree exactly), and it is why counts in the 2012 analysis differ from Table 1 by up to six.
- **What went wrong in 2012.** The original analysis computed p-values as `2*pnorm(-abs(cr))`, the normal approximation. That is the source of the apparent anomaly.
- **Unrecorded source.** The source of the two runs in `other_runs.npz` was not recorded. They are included for completeness and fit the binomial as well as the rest.

## Reproducing the results

Python 3.9 or later.

```
pip install -r requirements.txt
python code/analyse.py                    # a minute or two; writes results/summary.md and CSVs
python code/extra_checks.py               # writes results/extra_checks.md
python code/iversen_check.py 20000 10000  # about 25 minutes; results/iversen_battery.txt
python code/figure1.py                    # results/figure1.png and figure1.pdf (needs matplotlib)
```

To rebuild `data/` from the original workspaces: `python code/extract.py /folder/with/RData/files`.

`iversen_check.py` uses a fixed seed (2026), so `python code/iversen_check.py 20000 10000` reproduces the figures in `results/iversen_battery.txt` exactly (only the progress and timing lines differ).

## Main results

Full tables are in `results/summary.md` and `results/extra_checks.md`.

- **Main series, against the exact binomial.** At the seven thresholds in Table 1 (|CR| >= 2.0 to 5.0), no count differs from its exact expectation by more than 1.1 standard errors. Across fifteen thresholds the largest discrepancy is 2.1% at 3.5 (2.6 SE); a discrepancy that large somewhere among the fifteen occurs in 9% of simulated binomial series.
- **Main series, against the normal approximation.** Errors run in both directions: 10% too few extremes predicted at |CR| >= 2.5, 13% too many at 2.6. Roughly every 0.85 in the threshold the error jumps from negative to positive, because the lattice of possible CRs (spacing 0.0894) drifts against a 0.1 grid of thresholds; between the jumps it drifts back down through zero.
- **Every source and design.** The main series fits the exact binomial closely, and so do nearly all the others: variances within 0.3%, counts of |CR| >= 3 within 3%. Three of 35 chi-square fits give p < .05 (1.75 expected). One exploratory comparison, the +6 displacement of the random.org digits, departs significantly (p = .001, surviving a Bonferroni correction for 35 fits); it is not concentrated in the tails and remains unexplained.
- **Digit matching (p = 1/10).** The normal curve badly undercounts the upper tail. At CR >= 3.5 the Entropy Key gave 461 experiments (exact 455, normal 233) and the random.org digits 835 (exact 832, normal 258).
- **Does the standard error vary?** No. The daily variances are exactly as dispersed as chance predicts (p = .49), and the kurtosis of the total score limits any experiment-to-experiment variation in the standard error to about 1%. The difference score shows a small central excess (z = 3.0 within one hit of *M* = *A*) not shared by the total score; it is reported in `results/extra_checks.md`.
- **Iversen et al. (1971).** Their test battery is well calibrated under ideal randomness. Their dice (16 extreme scores of 267, 12.3 expected, P = .21) and pseudo-random series (5, 3.11 expected, P = .22) are within chance. Only the RAND digits are unusual (10, P = .003). Probabilities come from 20,000 simulated 100-block data sets and 10,000 simulated dice sets, with exact 95% intervals in the output.

## Licence and citation

Code (`code/`): MIT licence, see `LICENSE`. Data and results (`data/`, `results/`): CC BY 4.0, see `LICENSE-DATA`.

Please cite as: Savva, L. (2026). *Placebo ESP experiments: data and code* [Data set]. https://github.com/louiesavva/placebo-esp-randomness

## Use of AI tools

The re-analysis, the code and the drafting of both papers were carried out with the assistance of Claude (Anthropic). The author directed the work and supplied the original data.
