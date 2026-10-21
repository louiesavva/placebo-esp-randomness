# make_rdata.R - build data/placebo_esp.RData, the whole data set in one R file.
#
#     Rscript code/make_rdata.R
#
# Reads the .npz files in data/ (zip archives of .npy arrays) with base R only and saves one
# object per series. Hit counts are stored as integers. The values are identical to the .npz
# files; code/check_rdata.py confirms this.
#
#   main               morning, afternoon: hits of 250 each, 38,016,682 experiments (random.org)
#   main_files         the 2,266 daily random.org files, in order (16,777 experiments each)
#   short              morning, afternoon: hits of 125 each, 2,266,000 experiments (random.org)
#   prng               morning, afternoon: hits of 250 each, 1,000,004 experiments (Mersenne Twister in R)
#   entropy_key        straight (of 500), plus1, minus1 (of 490): digit matches, p = 1/10, 1,000,000 experiments
#   digits             straight, plus1..plus7, minus1..minus7 (of 268): digit matches, p = 1/10, 1,108,074 experiments
#   other_mcr          morning, afternoon (of 250): 1,000,000 experiments, source not recorded
#   other_aao1rng      morning, afternoon (of 250): 255,732 experiments, a generator on the author's computer
#   zener              hits of 625, p = 1/5, 473,594 experiments (random.org)
#   zener_files        the 2,266 daily files in the order processed (209 experiments each)
#   bits_2012_06_17,   the 8,388,500 bits used from each of two days of the main series, in order:
#   bits_2012_08_04    experiment k uses bits 500(k-1)+1 to 500k, the first 250 as morning
#   about              a short description of each object
#
# Load with:  load("data/placebo_esp.RData")

args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
data_dir <- file.path(here, "..", "data")

read_npz <- function(file) {
  tmp <- tempfile(); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))
  utils::unzip(file.path(data_dir, file), exdir = tmp)
  out <- list()
  for (f in list.files(tmp, pattern = "\\.npy$")) {
    con <- file(file.path(tmp, f), "rb")
    magic <- readBin(con, "raw", 6)
    stopifnot(identical(magic, as.raw(c(0x93, charToRaw("NUMPY")))))
    ver <- readBin(con, "integer", 2, size = 1, signed = FALSE)
    hlen <- if (ver[1] == 1) readBin(con, "integer", 1, size = 2, signed = FALSE, endian = "little")
            else readBin(con, "integer", 1, size = 4, endian = "little")
    header <- rawToChar(readBin(con, "raw", hlen))
    descr <- sub(".*'descr': *'([^']+)'.*", "\\1", header)
    shape <- sub(".*'shape': *\\(([^)]*)\\).*", "\\1", header)
    n <- if (nchar(gsub("[ ,]", "", shape)) == 0) 1 else as.numeric(strsplit(gsub(" ", "", shape), ",")[[1]][1])
    if (descr == "|u1") {
      val <- readBin(con, "integer", n, size = 1, signed = FALSE)
    } else if (descr == "<i8" && n == 1) {                   # a single 64-bit integer
      w <- readBin(con, "integer", 2, size = 4, endian = "little")
      val <- (w[1] %% 2^32) + w[2] * 2^32
    } else stop("unexpected dtype ", descr, " in ", file)
    close(con)
    stopifnot(length(val) == n)
    out[[sub("\\.npy$", "", f)]] <- val
  }
  out
}

frame <- function(file, cols) {
  x <- read_npz(file)
  as.data.frame(lapply(stats::setNames(cols, cols), function(k) as.integer(x[[k]])))
}
unpack <- function(file) {
  x <- read_npz(file)
  bits <- as.integer(rawToBits(as.raw(x$bits)))           # rawToBits gives the low bit first
  bits <- as.vector(matrix(bits, nrow = 8)[8:1, ])          # numpy packbits stores the high bit first
  bits[seq_len(x$n_bits)]
}
lines_of <- function(file) readLines(file.path(data_dir, file))

main <- data.frame(morning = as.integer(read_npz("random_org_main_morning.npz")$morning),
                   afternoon = as.integer(read_npz("random_org_main_afternoon.npz")$afternoon))
main_files <- lines_of("random_org_files.txt")
short <- frame("random_org_125.npz", c("morning", "afternoon"))
prng <- frame("prng.npz", c("morning", "afternoon"))
entropy_key <- frame("entropy_key_digits.npz", c("straight", "plus1", "minus1"))
digits <- frame("random_org_digits.npz", c("straight", paste0("plus", 1:7), paste0("minus", 1:7)))
o <- read_npz("other_runs.npz")
other_mcr <- data.frame(morning = as.integer(o$mcr_morning), afternoon = as.integer(o$mcr_afternoon))
other_aao1rng <- data.frame(morning = as.integer(o$aao1rng_morning), afternoon = as.integer(o$aao1rng_afternoon))
zener <- frame("random_org_zener.npz", "hits")
zener_files <- lines_of("random_org_zener_files.txt")
bits_2012_06_17 <- unpack("random_org_2012-06-17_bits.npz")
bits_2012_08_04 <- unpack("random_org_2012-08-04_bits.npz")

about <- c(
  main = "random.org, morning-afternoon design: hits of 250 in each half. X = morning - afternoon + 250 is Binomial(500, 1/2).",
  main_files = "Daily random.org files for main, in order; experiment i (from 1) comes from file (i - 1) %/% 16777 + 1.",
  short = "random.org, morning-afternoon design with 125 + 125 trials.",
  prng = "Mersenne Twister in R, morning-afternoon design, 250 + 250 trials.",
  entropy_key = "Simtec Entropy Key, digit matching: straight of 500 pairs; plus1 and minus1 of 490 pairs; p = 1/10.",
  digits = "random.org digits, digit matching of 268 pairs, displacements +1 to +7 and -1 to -7; p = 1/10.",
  other_mcr = "Morning-afternoon run, 250 + 250 trials; source not recorded.",
  other_aao1rng = "Morning-afternoon run, 250 + 250 trials; a generator on the author's computer, not otherwise recorded.",
  zener = "random.org Zener-card series: hits of 625 (25 runs of 25), p = 1/5.",
  zener_files = "Daily random.org binary files for zener, in the order processed; experiment i comes from file (i - 1) %/% 209 + 1.",
  bits_2012_06_17 = "The bits used from random.org's file for 17 June 2012, as divided into experiments 500 at a time.",
  bits_2012_08_04 = "The bits used from random.org's file for 4 August 2012, as divided into experiments 500 at a time.")

stopifnot(nrow(main) == 38016682, length(main_files) == 2266, nrow(short) == 2266000, nrow(prng) == 1000004,
          nrow(entropy_key) == 1e6, nrow(digits) == 1108074, nrow(other_mcr) == 1e6, nrow(other_aao1rng) == 255732,
          nrow(zener) == 473594, length(bits_2012_06_17) == 8388500, length(bits_2012_08_04) == 8388500)
k <- (match("2012-06-17.txt", main_files) - 1) * 16777
m <- matrix(bits_2012_06_17, nrow = 500)
stopifnot(all(colSums(m[1:250, ]) == main$morning[k + 1:16777]), all(colSums(m[251:500, ]) == main$afternoon[k + 1:16777]))

out <- file.path(data_dir, "placebo_esp.RData")
save(main, main_files, short, prng, entropy_key, digits, other_mcr, other_aao1rng, zener, zener_files,
     bits_2012_06_17, bits_2012_08_04, about, file = out, compress = "xz")
cat("written", out, round(file.size(out) / 2^20, 1), "MB\n")
