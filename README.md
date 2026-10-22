# Placebo ESP experiments and the normal approximation

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23541088.svg)](https://doi.org/10.5281/zenodo.23541088)

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
  zener.py           the Zener-card series  ->  results/zener.md
  rebuild_from_archive.py  rebuilds main-series experiments from random.org's public archive files
  make_rdata.R       data/*.npz  ->  data/placebo_esp.RData (the whole data set as one R file)
  check_rdata.py     confirms the R file holds exactly the same values as the .npz files
  figure1.py         Figure 1 (black and white): error of the normal approximations by threshold
  find_rdata.py      utility used to search the author's disks for the original files
data/                hit counts per experiment: numpy .npz files (about 70 MB) and the same data as one R file
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
| `random_org_zener.npz` | Zener cards | random.org (the same 2,266 days, binary files) | 625 trials (25 runs of 25), guess and target each one of five symbols, p = 1/5 | 473,594 | `hits` (of 625) | `rhine1.RData` (2013) | not used |
| `random_org_2012-06-17_bits.npz`, `random_org_2012-08-04_bits.npz` | Raw bits | random.org | the bits of two days of the main series, as divided into experiments | 2 × 16,777 | `bits` (packed), `n_bits` | `x.RData`, `rled.RData` | not used |

`random_org_files.txt` lists the 2,266 daily random.org files used by the main series, 12 March 2006 to 4 August 2012, in order. Each file supplied 8,388,608 random bits, divided into 16,777 blocks of 500: one bit per trial (250 morning, 250 afternoon), serving as the hit indicator. This is equivalent to matching a random target with a random guess, since the match of two independent fair bits is itself a fair bit. 108 bits of each file were not used (see below). Experiment *i* comes from file *i* // 16,777. The bits of two of these days (17 June and 4 August 2012) survive in October 2012 workspaces, already divided into experiments: taken 500 at a time, the first 250 as morning and the next 250 as afternoon, they reproduce those days' experiments exactly (`results/extra_checks.md`, section 5). random.org still distributes these daily files at [archive.random.org/text](https://archive.random.org/text) (free via BitTorrent, one torrent per month); each is 8,388,608 characters, each 0 or 1. Comparing the archive copy of 17 June 2012 with that day's surviving bits shows that the main series used each file's characters in order, from the 108th to the 8,388,607th: the first 107 characters and the last one are the 108 bits not used. `code/rebuild_from_archive.py` rebuilds the main series from the archive files and checks it against `data/`. Run on all 2,266 archive files (21 October 2026), it rebuilt every one of the 38,016,682 experiments exactly.

The Zener-card series reread the same 2,266 days from random.org's binary files. Each file was read as 262,144 32-bit integers, each reduced modulo 5 to one of five symbols; trial *j* matched symbol *j* (the guess) against symbol *j* + 131,072 (the target). The first 447 pairs of each file were not used, leaving 209 experiments of 625 trials. `random_org_zener_files.txt` lists the files in the order they were processed: experiment *i* comes from line *i* // 209 + 1. Reducing a 32-bit integer modulo 5 makes one symbol more likely than the others by about one part in 860 million, which is negligible.

Every series is also in one R file, `data/placebo_esp.RData` (69 MB), with one data frame per series, the file lists, the surviving bits and a short description of each object:

```r
load("data/placebo_esp.RData")
about                                      # what each object holds
x <- main$morning - main$afternoon + 250   # exactly Binomial(500, 1/2) under chance
```

It is built from the .npz files by `code/make_rdata.R`, and `code/check_rdata.py` confirms that every value is identical.

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
- **Rebuilt from public data.** Every experiment in the main series has been rebuilt from random.org's public archive of the original daily files and matches `data/` exactly (`code/rebuild_from_archive.py`).
- **Records from 2012.** A workspace saved in October 2012 holds frequency tables of morning hits, afternoon hits and their difference for all 38,016,682 experiments; they agree exactly with `data/` (checked by `extract.py`).
- **Unrecorded source.** The source of the two runs in `other_runs.npz` was not recorded; `aao1rng` used a random number generator on a computer the author had at the time. They are included for completeness and fit the binomial as well as the rest.

## Reproducing the results

Python 3.9 or later.

```
pip install -r requirements.txt
python code/analyse.py                    # a minute or two; writes results/summary.md and CSVs
python code/extra_checks.py               # writes results/extra_checks.md
python code/iversen_check.py 20000 10000  # about 25 minutes; results/iversen_battery.txt and iversen_replicates_20000_10000.npz
python code/figure1.py                    # results/figure1.png and figure1.pdf (needs matplotlib)
python code/zener.py                      # results/zener.md
Rscript code/make_rdata.R                 # data/placebo_esp.RData (needs R; about 8 minutes)
python code/check_rdata.py                # confirms the R file matches the .npz files
```

To rebuild `data/` from the original workspaces: `python code/extract.py /folder/with/RData/files`.

`iversen_check.py` uses a fixed seed (2026), so `python code/iversen_check.py 20000 10000` reproduces the figures in `results/iversen_battery.txt` exactly (only the progress and timing lines differ). The result for each simulated data set (its number of extreme scores, and for the 100-block sets the variance of its covariance scores) is saved in `results/iversen_replicates_20000_10000.npz`, so the probabilities can be checked without rerunning the simulation.

## Main results

Full tables are in `results/summary.md` and `results/extra_checks.md`.

- **Main series, against the exact binomial.** At the seven thresholds in Table 1 (|CR| >= 2.0 to 5.0), no count differs from its exact expectation by more than 1.1 standard errors. Across fifteen thresholds the largest discrepancy is 2.1% at 3.5 (2.6 SE); a discrepancy that large somewhere among the fifteen occurs in 9% of simulated binomial series.
- **Main series, against the normal approximation.** Errors run in both directions: 10% too few extremes predicted at |CR| >= 2.5, 13% too many at 2.6. Roughly every 0.85 in the threshold the error jumps from negative to positive, because the lattice of possible CRs (spacing 0.0894) drifts against a 0.1 grid of thresholds; between the jumps it drifts back down through zero.
- **Every source and design.** The main series fits the exact binomial closely, and so do nearly all the others: variances within 0.3%, counts of |CR| >= 3 within 3%. Three of 35 chi-square fits give p < .05 (1.75 expected). One exploratory comparison, the +6 displacement of the random.org digits, departs significantly (p = .001, surviving a Bonferroni correction for 35 fits); it is not concentrated in the tails and remains unexplained.
- **Digit matching (p = 1/10).** The normal curve badly undercounts the upper tail. At CR >= 3.5 the Entropy Key gave 461 experiments (exact 455, normal 233) and the random.org digits 835 (exact 832, normal 258).
- **Zener cards (p = 1/5, 625 trials).** The fit to the exact binomial is good (chi-square p = .67). Here the standard deviation is exactly 10, so every critical-ratio threshold on a 0.1 grid falls on a possible score and the plain normal curve always predicts too few extremes. The exact binomial expects 12% more than the normal curve at |CR| >= 2, 18% more at 3 and 23% more at 3.5, and the observed counts follow the exact values. In the upper tail, the direction an ESP test looks for, 201 experiments reached CR >= 3.5 (exact 189, normal 110).
- **Does the standard error vary?** No. The daily variances are exactly as dispersed as chance predicts (p = .49), and the kurtosis of the total score limits any experiment-to-experiment variation in the standard error to about 1%. The difference score shows a small central excess (z = 3.0 within one hit of *M* = *A*) not shared by the total score; it is reported in `results/extra_checks.md`.
- **Dependence between successive bits.** All 19,008,341,000 bits the main series used, rebuilt from random.org's archive, were tested file by file (`code/sequential_checks.py`, `results/sequential_checks.md`). Across the series, runs, run lengths and agreement at lags 2 to 16, 250 and 500 are as expected for independent fair bits. Day by day, the runs and the lag-2 agreement vary more than chance allows, because of a handful of days: 4 January 2012 (correlation between successive bits about -0.009, and about -0.003 at lag 2), 20 September 2011 (about -0.003 at lag 1) and 19 April 2008 (about +0.003 at lag 2), plus three marginal days. Without them the rest are close to chance. Correlations of this size change the variance of a 250-trial hit count by about 2% on the worst day and are not visible in the hit counts.
- **Iversen et al. (1971).** Their test battery is well calibrated under ideal randomness. Their dice (16 extreme scores of 267, 12.3 expected, P = .21) and pseudo-random series (5, 3.11 expected, P = .22) are within chance. Only the RAND digits are unusual (10, P = .003). Probabilities come from 20,000 simulated 100-block data sets and 10,000 simulated dice sets, with exact 95% intervals in the output.

## Licence and citation

Code (`code/`): MIT licence, see `LICENSE`. Data and results (`data/`, `results/`): CC BY 4.0, see `LICENSE-DATA`.

Please cite as: Savva, L. (2026). *Placebo ESP experiments: data and code* (Version 1.0.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.23541088

Version 1.0.0 is archived permanently at Zenodo; this repository is where the data and code are maintained.

## Use of AI tools

The re-analysis, the code and the drafting of both papers were carried out with the assistance of Claude (Anthropic). The author directed the work and supplied the original data.
