#!/usr/bin/env Rscript

# Extract RT x IM x fragment m/z cubes from Bruker diaPASEF .d data.
# This is intended as an MS2DecR-3D input candidate: do not preselect
# fragment transitions; keep the centroid peaks inside the precursor DIA
# window, RT range, IM range, and fragment m/z range.
#
# Run:
#   Rscript scripts/extract_fig6_3d_cube_opentimsr.R

options(stringsAsFactors = FALSE)

required <- c("opentimsr", "data.table")
missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing package(s): ", paste(missing_pkgs, collapse = ", "))
}

library(data.table)

d_path <- Sys.getenv(
  "FIG6_D_PATH",
  "raw/fig6_run.d"
)
out_dir <- Sys.getenv("FIG6_CUBE_OUT_DIR", "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

precursor_mz <- as.numeric(Sys.getenv("FIG6_PRECURSOR_MZ", "829.8746"))
fragment_mz_range <- c(
  as.numeric(Sys.getenv("FIG6_FRAGMENT_MZ_LOWER", "150")),
  as.numeric(Sys.getenv("FIG6_FRAGMENT_MZ_UPPER", "1500"))
)

rt_bin_width <- as.numeric(Sys.getenv("FIG6_RT_BIN_WIDTH_MIN", "0.01"))
im_bin_width <- as.numeric(Sys.getenv("FIG6_IM_BIN_WIDTH", "0.001"))
mz_bin_width <- as.numeric(Sys.getenv("FIG6_MZ_BIN_WIDTH", "0.01"))
chunk_size <- as.integer(Sys.getenv("FIG6_FRAME_CHUNK_SIZE", "80"))

windows <- data.table(
  label = c("both_deconv", "S1_box", "S8_box"),
  rt_min_lower = c(17.60, 18.40, 17.75),
  rt_min_upper = c(18.70, 18.50, 17.85),
  im_lower = c(0.930, 0.932, 0.953),
  im_upper = c(0.970, 0.948, 0.968)
)

custom_window <- Sys.getenv("FIG6_CUBE_WINDOW", "")
if (nzchar(custom_window)) {
  parts <- strsplit(custom_window, ",", fixed = TRUE)[[1]]
  if (length(parts) != 5) {
    stop("FIG6_CUBE_WINDOW must be: label,rt_lower,rt_upper,im_lower,im_upper")
  }
  windows <- data.table(
    label = parts[1],
    rt_min_lower = as.numeric(parts[2]),
    rt_min_upper = as.numeric(parts[3]),
    im_lower = as.numeric(parts[4]),
    im_upper = as.numeric(parts[5])
  )
}

table_dt <- function(handle, table_name) {
  x <- opentimsr::table2df(handle, table_name)
  if (is.list(x) && !is.data.frame(x)) {
    x <- x[[table_name]]
  }
  as.data.table(x)
}

bin_floor <- function(x, width) {
  floor(x / width) * width
}

extract_window <- function(D, frames, prec_windows, spec) {
  rt_range_sec <- c(spec$rt_min_lower, spec$rt_min_upper) * 60
  frame_ids <- frames[Time >= rt_range_sec[1] & Time <= rt_range_sec[2], Id]
  if (length(frame_ids) == 0) {
    stop("No frames in RT window for ", spec$label)
  }

  chunks <- split(frame_ids, ceiling(seq_along(frame_ids) / chunk_size))
  all_peaks <- vector("list", length(chunks))
  query_columns <- c("frame", "scan", "tof", "intensity", "mz", "inv_ion_mobility", "retention_time")

  for (i in seq_along(chunks)) {
    message(sprintf("[%s] query chunk %d/%d (%d frames)", spec$label, i, length(chunks), length(chunks[[i]])))
    peaks <- as.data.table(opentimsr::query(D, frames = chunks[[i]], columns = query_columns))

    in_scan_window <- rep(FALSE, nrow(peaks))
    for (j in seq_len(nrow(prec_windows))) {
      in_scan_window <- in_scan_window |
        (peaks$scan >= prec_windows$ScanNumBegin[j] & peaks$scan <= prec_windows$ScanNumEnd[j])
    }
    peaks <- peaks[in_scan_window]

    peaks <- peaks[
      retention_time >= rt_range_sec[1] &
        retention_time <= rt_range_sec[2] &
        inv_ion_mobility >= spec$im_lower &
        inv_ion_mobility <= spec$im_upper &
        mz >= fragment_mz_range[1] &
        mz <= fragment_mz_range[2] &
        intensity > 0
    ]

    if (nrow(peaks) > 0) {
      peaks[, rt_min := retention_time / 60]
      peaks[, `:=`(
        rt_bin_min = bin_floor(rt_min, rt_bin_width),
        im_bin = bin_floor(inv_ion_mobility, im_bin_width),
        mz_bin = bin_floor(mz, mz_bin_width)
      )]
      all_peaks[[i]] <- peaks
    }
  }

  peaks <- rbindlist(all_peaks, use.names = TRUE, fill = TRUE)
  if (nrow(peaks) == 0) {
    warning("No peaks extracted for ", spec$label)
    return(NULL)
  }

  cube <- peaks[
    ,
    .(
      intensity = sum(intensity),
      n_peaks = .N,
      mz_weighted = sum(mz * intensity) / sum(intensity),
      rt_weighted_min = sum(rt_min * intensity) / sum(intensity),
      im_weighted = sum(inv_ion_mobility * intensity) / sum(intensity)
    ),
    by = .(rt_bin_min, im_bin, mz_bin)
  ]
  setorder(cube, rt_bin_min, im_bin, mz_bin)

  summary <- data.table(
    label = spec$label,
    rt_min_lower = spec$rt_min_lower,
    rt_min_upper = spec$rt_min_upper,
    im_lower = spec$im_lower,
    im_upper = spec$im_upper,
    precursor_mz = precursor_mz,
    fragment_mz_lower = fragment_mz_range[1],
    fragment_mz_upper = fragment_mz_range[2],
    raw_peak_rows = nrow(peaks),
    cube_rows = nrow(cube),
    total_intensity = sum(peaks$intensity),
    max_peak_intensity = max(peaks$intensity),
    max_cube_intensity = max(cube$intensity)
  )

  fwrite(peaks, file.path(out_dir, paste0(spec$label, "_raw_peaks.tsv")), sep = "\t")
  fwrite(cube, file.path(out_dir, paste0(spec$label, "_cube_binned.tsv")), sep = "\t")
  fwrite(summary, file.path(out_dir, paste0(spec$label, "_summary.tsv")), sep = "\t")
  summary
}

# Require the Bruker conversion used to validate this dataset.
bruker_dll <- Sys.getenv(
  "BRUKER_TIMSDATA_DLL",
  "lib/Bruker/timsdata.dll"
)
if (file.exists(bruker_dll)) {
  message("Using Bruker timsdata.dll for m/z and ion-mobility conversion: ", bruker_dll)
  opentimsr::setup_bruker_so(bruker_dll)
} else {
  stop("Set BRUKER_TIMSDATA_DLL to the MS-DIAL bundled Bruker timsdata.dll.")

}

message("Opening .d: ", d_path)
D <- opentimsr::OpenTIMS(d_path)

frames <- table_dt(D, "Frames")
dia_info <- table_dt(D, "DiaFrameMsMsInfo")
dia_windows <- table_dt(D, "DiaFrameMsMsWindows")
dia_windows[, mz_lower := IsolationMz - IsolationWidth / 2]
dia_windows[, mz_upper := IsolationMz + IsolationWidth / 2]

prec_windows <- dia_windows[mz_lower <= precursor_mz & mz_upper >= precursor_mz]
if (nrow(prec_windows) == 0) {
  stop("No DIA isolation window contains precursor m/z ", precursor_mz)
}

target_frame_ids <- dia_info[WindowGroup %in% prec_windows$WindowGroup, unique(Frame)]
target_frames <- frames[Id %in% target_frame_ids]

fwrite(windows, file.path(out_dir, "cube_windows.tsv"), sep = "\t")
fwrite(prec_windows, file.path(out_dir, "precursor_containing_windows.tsv"), sep = "\t")

summaries <- vector("list", nrow(windows))
for (i in seq_len(nrow(windows))) {
  summaries[[i]] <- extract_window(D, target_frames, prec_windows, windows[i])
}

summary_all <- rbindlist(summaries, use.names = TRUE, fill = TRUE)
fwrite(summary_all, file.path(out_dir, "cube_summary.tsv"), sep = "\t")

message("Done: ", normalizePath(out_dir, winslash = "/"))
