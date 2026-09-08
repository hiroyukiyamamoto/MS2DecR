#!/usr/bin/env Rscript
# Run from data/tims. Preserves measured frame times and sums over IM.
suppressPackageStartupMessages(library(data.table))
out_dir <- Sys.getenv('FIG6_CUBE_OUT_DIR', 'output')
input <- file.path(out_dir, 'both_deconv_raw_peaks.tsv')
if (!file.exists(input)) input <- paste0(input, '.gz')
read_table <- function(path) {
  if (grepl('[.]gz$', path)) {
    con <- gzfile(path, 'rt')
    on.exit(close(con))
    return(as.data.table(read.delim(con, check.names = FALSE)))
  }
  fread(path)
}
p <- read_table(input)
stopifnot(nrow(p) > 0, all(is.finite(p$intensity)), all(p$intensity > 0))
x <- p[, .(intensity = sum(intensity)), by = .(frame, retention_time, mz_bin)]
axis_rt <- unique(x[, .(frame, retention_time)])
setorder(axis_rt, retention_time, frame)
axis_mz <- sort(unique(x$mz_bin))
stopifnot(!anyDuplicated(axis_rt$frame))
Y <- matrix(0, nrow(axis_rt), length(axis_mz))
Y[cbind(match(x$frame, axis_rt$frame), match(x$mz_bin, axis_mz))] <- x$intensity
demo <- list(
  mz = axis_mz, Y = Y, rt = axis_rt$retention_time, frame = axis_rt$frame,
  metadata = list(
    accession = 'PXD033904',
    run = '20220302_tims1_nElute_8cm_DOl_Phospho_60min_rep1_Slot1-95_1_1835',
    precursor_mz = 829.8746, rt_range_min = c(17.6, 18.7),
    im_range = c(0.93, 0.97), fragment_mz_range = c(150, 1500),
    mz_bin_width = 0.01, mz_axis = 'bin lower edge', rt_unit = 'seconds',
    aggregation = 'sum over IM within each measured frame',
    conversion = 'Bruker timsdata.dll bundled with MS-DIAL 5.5.260319'
  )
)
stopifnot(nrow(Y) == length(demo$rt), ncol(Y) == length(demo$mz),
          all(is.finite(Y)), all(Y >= 0), all(diff(demo$rt) > 0),
          isTRUE(all.equal(sum(Y), sum(p$intensity))))
saveRDS(demo, file.path(out_dir, 'fig6_ms2decr_demo_rt_mz.rds'), compress = 'xz')

# Check range and intensity conservation for all three extracted windows.
windows <- fread(file.path(out_dir, 'cube_windows.tsv'))
checks <- lapply(seq_len(nrow(windows)), function(i) {
  w <- windows[i]
  find_table <- function(suffix) {
    f <- file.path(out_dir, paste0(w$label, suffix, '.tsv'))
    if (!file.exists(f)) f <- paste0(f, '.gz')
    read_table(f)
  }
  peaks <- find_table('_raw_peaks')
  cube <- find_table('_cube_binned')
  stopifnot(nrow(peaks) > 0, nrow(cube) > 0,
    all(peaks$rt_min >= w$rt_min_lower & peaks$rt_min <= w$rt_min_upper),
    all(peaks$inv_ion_mobility >= w$im_lower & peaks$inv_ion_mobility <= w$im_upper),
    all(peaks$mz >= 150 & peaks$mz <= 1500),
    all(is.finite(peaks$intensity)), all(peaks$intensity > 0),
    isTRUE(all.equal(sum(peaks$intensity), sum(cube$intensity))))
  data.table(label = w$label, raw_rows = nrow(peaks), cube_rows = nrow(cube),
             total_intensity = sum(peaks$intensity), range_check = TRUE,
             intensity_conserved = TRUE)
})
fwrite(rbindlist(checks), file.path(out_dir, 'validation.tsv'), sep = '\t')
png(file.path(out_dir, 'demo_rt_profile.png'), width = 1000, height = 600)
plot(demo$rt / 60, rowSums(Y), type = 'l', xlab = 'RT (min)',
     ylab = 'Summed MS2 intensity', main = 'Fig.6 raw MS2: IM 0.930-0.970')
dev.off()

# Compact text data with base R; no additional compression package required.
compress_table <- function(path) {
  src <- file(path, 'rb')
  dst <- gzfile(paste0(path, '.gz'), 'wb', compression = 9)
  on.exit({close(src); close(dst)})
  repeat {
    block <- readBin(src, 'raw', n = 1048576L)
    if (!length(block)) break
    writeBin(block, dst)
  }
}
tables <- list.files(out_dir, pattern = '_(raw_peaks|cube_binned)[.]tsv$', full.names = TRUE)
for (path in tables) compress_table(path)
cat('Validated demo:', nrow(Y), 'frames x', ncol(Y), 'm/z bins\n')
cat('Uncompressed TSVs may be removed after checking the .tsv.gz files.\n')
