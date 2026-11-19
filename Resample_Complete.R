#------------------------------------------------------------------------------#
#--------------- RESAMPLE THE AIR TEMPERATURE RASTERS -------------------------#
#------------------------------------------------------------------------------#

library(raster)
library(terra)
library(dplyr)
library(tidyr)
library(readr)
library(tidyverse)
library(plyr)
library(rgdal)

# Define input and output directories
input_dir   <- "E:/PhD/Rasters/MSc Rasters/Air Temperature/RESOLUTION/AT5cm"
output_base <- "E:/PhD/Rasters/MSc Rasters/Air Temperature/RESOLUTION"

setwd(input_dir)
raster_list <- list.files(pattern = "\\.tif$")

# Target resolutions (m)
resolutions <- c(
  seq(0.1, 2, by = 0.10),   # 5 cm to 2 m in 10 cm steps
  seq(3, 100, by = 1)       # 3 m to 100 m in 1 m steps
)

# Assume 5 cm nominal native resolution (m)
nominal_orig_res <- 0.05  

# *** IMPORTANT: set your nodata/background codes here ***
# If your outside-buffer / nodata is 0 or -9999, list them:
nodata_vals <- c(0, -9999)

# Log of original resolutions
res_log <- data.frame(
  file       = character(),
  orig_res_x = numeric(),
  orig_res_y = numeric(),
  stringsAsFactors = FALSE
)

for (res_value in resolutions) {
  
  # Create output folder for the current resolution
  output_folder <- file.path(output_base, paste0("AT", res_value))
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  
  for (raster_file in raster_list) {
    
    # Read original raster
    air_raster <- raster(raster_file)
    
    # --- Handle nodata / background values ---
    # Convert specified nodata codes to NA so they are excluded from the mean
    for (nv in nodata_vals) {
      air_raster[air_raster == nv] <- NA
    }
    # Optional: if you *know* valid temps are e.g. 10–60 °C, you can also clip:
    # air_raster[air_raster < 10 | air_raster > 60] <- NA
    
    # Actual stored resolution (for logging only)
    r_xy       <- res(air_raster)   # c(xres, yres)
    orig_res_x <- r_xy[1]
    orig_res_y <- r_xy[2]
    
    # Log original resolution once per file
    bn <- basename(raster_file)
    if (!bn %in% res_log$file) {
      res_log <- rbind(
        res_log,
        data.frame(
          file       = bn,
          orig_res_x = orig_res_x,
          orig_res_y = orig_res_y,
          stringsAsFactors = FALSE
        )
      )
    }
    
    # Optional: warn if stored res is far from 0.05 m
    if (abs(orig_res_x - nominal_orig_res) > 0.01) {
      warning("File ", bn, " has x resolution ", orig_res_x,
              " (diff from 0.05 = ", orig_res_x - nominal_orig_res, ")")
    }
    
    # Aggregation factor based on nominal 0.05 m base res
    fact <- round(res_value / nominal_orig_res)
    if (fact < 1) {
      stop("Requested resolution ", res_value,
           " is smaller than nominal original resolution ", nominal_orig_res)
    }
    
    # Aggregate using mean of non-NA pixels
    resampled_raster <- aggregate(air_raster, fact = fact,
                                  fun = mean, na.rm = TRUE)
    
    # New filename: replace leading "AT" with "AT<res>"
    old_name <- basename(raster_file)
    new_name <- sub("^AT", paste0("AT", res_value), old_name)
    output_path <- file.path(output_folder, new_name)
    
    # Overwrite existing rasters (fix incorrect previous outputs)
    writeRaster(resampled_raster, output_path, overwrite = TRUE)
  }
}

# Write log of original resolutions
log_path <- file.path(output_base, "resolution_log.csv")
write.csv(res_log, log_path, row.names = FALSE)
message("[ok] Original resolutions logged to: ", log_path)

#------------------------------------------------------------------------------#
#--------------------- CALCULATE TPI RASTERS ----------------------------------#
#------------------------------------------------------------------------------#

library(terra)
library(tools)

# Define species-specific transformation functions
species_list <- list(
  ANTIOPA = function(x) (x + 19.1266157150269) / (34.6577774047852 + 19.1266157150269),
  APHRODITE = function(x) (x + 14.7759920120239) / (32.3702377319336 + 14.7759920120239),
  ARTHEMIS = function(x) (x + 21.3989551544189) / (37.1850135803223 + 21.3989551544189),
  ATALANTA = function(x) (x + 11.4172444343567) / (35.5383041381836 + 11.4172444343567),
  AUGUSTINUS = function(x) (x + 5.15009803771973) / (29.1876739501953 + 5.15009803771973),
  BAPTISIAE = function(x) (x - 3.1226978302002) / (33.3321212768555 - 3.1226978302002),
  BELLONA = function(x) (x + 13.1947817802429) / (31.7832191467285 + 13.1947817802429),
  BIMACULA = function(x) (x - 9.27992362976074) / (32.2608283996582 - 9.27992362976074),
  COCYTA = function(x) (x + 4.80765491724014) / (27.7682361602783 + 4.80765491724014),
  COMYNTAS = function(x) (x + 8.78979110717773) / (35.4972869873047 + 8.78979110717773),
  CYBELE = function(x) (x + 12.4393273353577) / (34.2306396484375 + 12.4393273353577),
  CYMELA = function(x) (x + 0.323590278625488) / (33.9261047363281 + 0.323590278625488),
  DION = function(x) (x - 8.6043044090271) / (34.9274810791016 - 8.6043044090271),
  EGEREMET = function(x) (x - 9.1182834148407) / (31.2812335968018 - 9.1182834148407),
  EURYTHEME = function(x) (x + 17.635590171814) / (39.2337394714355 + 17.635590171814),
  HARRISII = function(x) (x + 0.259680366516113) / (29.3815567016602 + 0.259680366516113),
  HENRICI = function(x) (x - 0.126528882980347) / (27.3422073364258 - 0.126528882980347),
  LINEOLA = function(x) (x - 4.67663316726685) / (30.1312366485596 - 4.67663316726685),
  LOGAN = function(x) (x - 8.97750358581543) / (26.9660892486572 - 8.97750358581543),
  LUCIA = function(x) (x + 12.9523355007172) / (28.2466384887695 + 12.9523355007172),
  LUCILIUS = function(x) (x - 3.72501909732819) / (33.0711490631104 - 3.72501909732819),
  LYGDAMUS = function(x) (x + 19.9541408538818) / (34.6273559570312 + 19.9541408538818),
  ORIGENES = function(x) (x - 1.66916532516479) / (33.538646697998 - 1.66916532516479),
  PEGALA = function(x) (x + 18.0576820373535) / (36.316365814209 + 18.0576820373535),
  PHILODICE = function(x) (x + 18.4147871017456) / (35.6953071594238 + 18.4147871017456),
  PLEXIPPUS = function(x) (x + 17.1790409088135) / (35.0375343322754 + 17.1790409088135),
  POLIOS = function(x) (x + 3.85859483480453) / (26.1612518310547 + 3.85859483480453),
  PROGNE = function(x) (x + 8.35500826835632) / (30.0811725616455 + 8.35500826835632),
  RAPAE = function(x) (x + 15.2569845199585) / (36.1623306274414 + 15.2569845199585),
  SASSACUS = function(x) (x - 5.40159873962402) / (28.5587604522705 - 5.40159873962402),
  SELENE = function(x) (x + 6.61254180669785) / (32.2023822784424 + 6.61254180669785),
  THAROS = function(x) (x + 11.8877655029297) / (35.9251731872559 + 11.8877655029297),
  THEMISTOCLES = function(x) (x + 5.56193317174912) / (36.8051856994629 + 5.56193317174912),
  TULLIA = function(x) (x + 12.6538070678711) / (36.0443229675292 + 12.6538070678711),
  VESTRIS = function(x) (x + 8.60373512879014) / (35.840274810791 + 8.60373512879014)
)

# Directories
base_in_dir  <- "E:/PhD/Rasters/MSc Rasters/Air Temperature/RESOLUTION"
base_out_dir <- "E:/PhD/Rasters/MSc Rasters/TPI/RESOLUTION"

# ---- SELECT AT FOLDERS TO PROCESS ----
# We want:
#  - all non-integer AT resolutions below 2 m (0 < res < 2 and not an integer)
#  - all integer AT resolutions >= 2 m (2, 3, 4, ... no explicit upper cap)

res_folders_all <- list.dirs(base_in_dir, full.names = TRUE, recursive = FALSE)
res_names       <- basename(res_folders_all)

# parse numeric resolution from "AT<res>"
res_vals <- suppressWarnings(as.numeric(sub("^AT", "", res_names)))

is_integer_like <- function(x, tol = 1e-6) is.finite(x) && (abs(x - round(x)) < tol)

keep <- vapply(
  res_vals,
  FUN = function(v) {
    if (!is.finite(v) || v <= 0) return(FALSE)
    (v < 2 && !is_integer_like(v)) || (v >= 2 && is_integer_like(v))
  },
  FUN.VALUE = logical(1)
)

res_folders <- res_folders_all[keep]

message("AT folders selected for TPI generation: ",
        paste(basename(res_folders), collapse = ", "))

# ---- MAIN LOOP: CREATE TPI FOLDERS + RASTERS ----

for (in_folder in res_folders) {
  res_name   <- basename(in_folder)          # e.g., "AT0.1", "AT1.9", "AT2", "AT51"
  res_num    <- sub("^AT", "", res_name)     # "0.1", "1.9", "2", "51"
  tpi_folder <- paste0("TPI", res_num)
  out_res_dir <- file.path(base_out_dir, tpi_folder)
  
  if (dir.exists(out_res_dir)) {
    message("Skipping ", res_name, ": TPI folder already exists (", tpi_folder, ").")
    next
  }
  
  message("Processing ", res_name, " ➝ creating ", tpi_folder)
  
  # Create folders
  dir.create(out_res_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_res_dir, "SITE"), showWarnings = FALSE)
  species_dir <- file.path(out_res_dir, "SPECIES")
  dir.create(species_dir, showWarnings = FALSE)
  for (species in names(species_list)) {
    dir.create(file.path(species_dir, species), showWarnings = FALSE)
  }
  
  # List rasters
  rasters <- list.files(in_folder, pattern = "\\.tif$", full.names = TRUE)
  message(length(rasters), " rasters found in ", res_name)
  
  for (r in rasters) {
    base_name <- tools::file_path_sans_ext(basename(r))
    message("  Processing raster: ", base_name)
    rast_in <- rast(r)
    
    for (species in names(species_list)) {
      message("    Applying species: ", species)
      fun      <- species_list[[species]]
      rast_out <- app(rast_in, fun)
      out_path <- file.path(species_dir, species,
                            paste0(base_name, "_", species, ".tif"))
      writeRaster(rast_out, out_path, overwrite = TRUE)
      message("      Saved: ", out_path)
    }
  }
  
  message("Finished processing ", tpi_folder, "\n")
}

message("✅ TPI generation complete for selected resolutions.")

#------------------------------------------------------------------------------#
#------------------- COPY AND ORGANIZE BY SITE --------------------------------#
#------------------------------------------------------------------------------#

# ---- Copy TPI/SPECIES rasters into TPI/SITE/<sXX>/ by site_month_day ----
# Expects filenames like: AT<res>_<site>_<month>_<day>_<SPECIES>.tif

library(tools)

BASE_TPI_ROOT <- "E:/PhD/Rasters/MSc Rasters/TPI/RESOLUTION"

# Mapping of <site>_<month>_<day> -> sXX   (NO ".tif" suffix)
date_to_site <- c(
  "1_5_18"="s11", "1_7_22"="s12",
  "2_5_18"="s21", "2_7_28"="s22",
  "3_5_17"="s31", "3_8_20"="s32",
  "4_5_17"="s41", "4_7_28"="s42",
  "5_5_19"="s51", "5_7_23"="s52",
  "6_5_24"="s61", "6_8_2" ="s62",
  "7_6_2" ="s71", "7_7_26"="s72",
  "8_6_9" ="s81", "8_8_16"="s82",
  "9_6_9" ="s91", "9_8_20"="s92",
  "10_7_21"="s101","10_8_26"="s102",
  "11_7_21"="s111","11_8_26"="s112",
  "12_6_10"="s121","12_8_16"="s122",
  "13_7_21"="s131","14_7_16"="s141","14_8_24"="s142",
  "15_5_31"="s151","15_8_3" ="s152"
)

# --- Helpers ---

# Extract "site_month_day" as "1_7_16" from "AT60_1_7_16_ANTIOPA.tif"
extract_date_key <- function(fname) {
  bn <- basename(fname)
  m  <- regexec("^AT[^_]*_(\\d+_\\d+_\\d+)_.*\\.tif$", bn, perl = TRUE)
  r  <- regmatches(bn, m)[[1]]
  if (length(r) >= 2) return(r[2])
  # Fallback if pattern slightly varies
  p <- strsplit(bn, "_", fixed = TRUE)[[1]]
  if (length(p) >= 5) return(paste(p[2:4], collapse = "_"))
  NA_character_
}

# Parse meters from TPI folder name
tpi_meters <- function(tpi_name){
  nm <- gsub("^TPI[ _-]*", "", tpi_name, ignore.case = TRUE)
  nm <- gsub(",", ".", nm, fixed = TRUE)
  if (grepl("cm$", nm, ignore.case = TRUE)) {
    return(as.numeric(sub("cm$", "", nm, ignore.case = TRUE)) / 100)
  }
  as.numeric(sub("m$", "", nm, ignore.case = TRUE))
}

# --- Choose which TPI folders to process ---

tpi_dirs_all <- list.dirs(BASE_TPI_ROOT, full.names = TRUE, recursive = FALSE)
tpi_dirs_all <- tpi_dirs_all[grepl("^TPI", basename(tpi_dirs_all), ignore.case = TRUE)]

tpi_info <- data.frame(
  path   = tpi_dirs_all,
  name   = basename(tpi_dirs_all),
  meters = sapply(basename(tpi_dirs_all), tpi_meters),
  stringsAsFactors = FALSE
)

# Process all TPI resolutions with positive meter value
res_filter_tpi <- function(m) is.finite(m) && (m > 0)
tpi_dirs <- tpi_info$path[sapply(tpi_info$meters, res_filter_tpi)]

if (!length(tpi_dirs)) stop("No TPI folders matched the resolution filter (m > 0).")

# --- Main loop: copy SPECIES -> SITE/sXX ---

for (tpi in tpi_dirs) {
  species_root <- file.path(tpi, "SPECIES")
  site_root    <- file.path(tpi, "SITE")
  
  if (!dir.exists(species_root)) { 
    message("[skip] No SPECIES/ in ", basename(tpi)) 
    next 
  }
  dir.create(site_root, showWarnings = FALSE, recursive = TRUE)
  
  tif_paths <- list.files(
    species_root, pattern = "\\.tif$", 
    recursive = TRUE, full.names = TRUE
  )
  if (!length(tif_paths)) { 
    message("[skip] No .tif under ", species_root) 
    next 
  }
  
  message("TPI folder: ", basename(tpi), 
          "  (#rasters: ", length(tif_paths), ")")
  
  n_ok <- 0L; n_warn <- 0L
  
  for (i in seq_along(tif_paths)) {
    src <- tif_paths[i]
    
    date_key <- extract_date_key(src)    # e.g. "1_7_16"
    if (is.na(date_key)) { 
      n_warn <- n_warn + 1L
      message("  [warn] Could not parse date key: ", basename(src))
      next 
    }
    
    site_id <- unname(date_to_site[date_key])
    if (is.na(site_id)) { 
      n_warn <- n_warn + 1L
      message("  [warn] No mapping for ", date_key, " (", basename(src), ")")
      next 
    }
    
    dest_dir <- file.path(site_root, site_id)
    dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)
    
    dest_path <- file.path(dest_dir, basename(src))
    
    if (i <= 3L) {
      message("    -> example dest: ", dest_path)
    }
    
    ok <- tryCatch(
      file.copy(src, dest_path, overwrite = TRUE),
      error = function(e) FALSE
    )
    
    if (isTRUE(ok)) {
      n_ok <- n_ok + 1L
    } else {
      n_warn <- n_warn + 1L
      message("  [warn] Copy failed: ", basename(src))
    }
  }
  
  message("  -> Copied (overwritten if existed): ", n_ok, 
          "  Warnings: ", n_warn)
}

message("✅ Done copying SPECIES rasters into per-resolution SITE/sXX folders.")

################################################################################
################################# TPI ##########################################
################################################################################
#------------------------------------------------------------------------------#
#------------- CREATE MEAN TPI BY SITE AND RESOLUTION CSV ---------------------#
#------------------------------------------------------------------------------#
# --- Libraries (lean) ---
library(dplyr)
library(readr)
library(data.table)
library(raster)
library(stringi)
library(tools)

# -------------------------
# Config: input / output
# -------------------------
BASE_RES_DIR <- "E:/PhD/Rasters/MSc Rasters/TPI/RESOLUTION"   # contains TPI0.1, TPI2, TPI32, ...
BASE_OUT     <- "E:/PhD/Chap 4 - Resolution"                   # <- all CSVs go here now
GLOBAL_MASTER_PATH <- file.path(BASE_OUT, "Resample_TPI_mean.csv")
dir.create(BASE_OUT, showWarnings = FALSE, recursive = TRUE)

# -------------------------
# Helpers
# -------------------------
# Mean TPI of a raster layer (ignore NA)
mean_tpi <- function(x){
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

# Parse token after "TPI" (supports integers, decimals, optional trailing 'm')
# Examples: TPI6 -> "6", TPI0.05 -> "0.05", TPI2m -> "2"
resolution_token <- function(res_folder_name){
  nm <- gsub("^TPI[ _-]*", "", res_folder_name, ignore.case = TRUE)
  nm <- trimws(nm)
  nm <- gsub(",", ".", nm, fixed = TRUE)
  m  <- regexec("^(\\d+(?:\\.\\d+)?)(?:m)?$", nm, perl = TRUE)
  r  <- regmatches(nm, m)[[1]]
  if (length(r) == 0) stop(sprintf("Could not parse resolution from '%s'", res_folder_name))
  r[2]  # numeric string in meters
}

# Convert token to meters label used in column name (strip any trailing 'm')
#   "6" -> "6"; "0.05" -> "0.05"; "2" -> "2"
meters_label <- function(token){
  token <- trimws(token)
  token <- gsub(",", ".", token, fixed = TRUE)
  sub("m$", "", token, ignore.case = TRUE)
}

tpi_colname <- function(m_label) paste0("TPI", m_label)

species_vec <- c(
  "antiopa","aphrodite","archippus","atalanta","augustinus","baptisiae",
  "bellona","bimacula","cocyta","comyntas","cybele","cymela","dion",
  "egeremet","eurytheme","henrici","lineola","lucia","lucilius","lygdamus",
  "origenes","pegala","philodice","plexippus","polios","progne","rapae",
  "sassacus","selene","tharos","themistocles","tullia","vestris"
)

# -------------------------
# Process one resolution dir
# -------------------------
process_resolution_dir <- function(res_dir){
  res_name  <- basename(res_dir)
  token     <- resolution_token(res_name)      # "6", "0.05", "2", ...
  m_label   <- meters_label(token)             # same as token here
  col_TPI   <- tpi_colname(m_label)            # "TPI6", "TPI0.05", ...
  
  site_root <- file.path(res_dir, "SITE")
  if (!dir.exists(site_root)){ 
    message(sprintf("[skip] No SITE/ in %s", res_dir)); 
    return(invisible(NULL)) 
  }
  
  # Put per-site & per-resolution master under BASE_OUT/TPI/<col_TPI>/
  out_dir_res <- file.path(BASE_OUT, "TPI", col_TPI)
  dir.create(out_dir_res, showWarnings = FALSE, recursive = TRUE)
  
  site_folders <- list.dirs(site_root, full.names = TRUE, recursive = FALSE)
  if (!length(site_folders)){ 
    message(sprintf("[warn] No site folders in %s", site_root)); 
    return(invisible(NULL)) 
  }
  
  # helper to map a *filename* (no path, no extension) to a species
  # expects strings like:
  #   AT0.1_1_5_18_ANTIOPA
  #   TPI10_2_7_18_SELENE
  map_species <- function(fname_noext){
    # just in case we ever get a full path or extension
    x <- tools::file_path_sans_ext(basename(fname_noext))
    # handle odd duplicates like ".1", ".2" if raster adds them
    x <- sub("\\.\\d+$", "", x)
    # special lucilius case (if you ever had "LUC_" prefix)
    x <- gsub("^LUC_", "LUCILIUS_", x, ignore.case = TRUE)
    
    # assume species is last underscore-delimited token
    parts <- strsplit(x, "_", fixed = TRUE)[[1]]
    if (length(parts) == 0) return(NA_character_)
    cand  <- parts[length(parts)]         # e.g. "ANTIOPA"
    cand  <- gsub("^LUC$", "LUCILIUS", cand, ignore.case = TRUE)
    
    sv    <- species_vec
    sv_l  <- tolower(sv)
    cand_l <- tolower(cand)
    
    hit <- match(cand_l, sv_l)
    if (!is.na(hit)) return(sv[hit])
    
    # fallback: maybe the whole string is species name
    x_l <- tolower(x)
    hit2 <- match(x_l, sv_l)
    if (!is.na(hit2)) return(sv[hit2])
    
    # fallback: species name at start of full string
    hit3 <- which(startsWith(x_l, sv_l))
    if (length(hit3) == 1) return(sv[hit3])
    
    NA_character_
  }
  
  site_csvs <- character(0)
  
  for (fp in site_folders){
    tif_list <- list.files(fp, pattern = "\\.tif$", full.names = TRUE)
    if (!length(tif_list)){ 
      message(sprintf("[warn] No .tif in %s", fp)); 
      next 
    }
    
    stk <- try(stack(tif_list), silent = TRUE)
    if (inherits(stk, "try-error")){ 
      message(sprintf("[error] Could not stack rasters in %s", fp)); 
      next 
    }
    
    dts <- as.data.frame(stk)
    
    # --- compute mean per raster (column) ---
    means <- vapply(dts, mean_tpi, numeric(1))
    
    # tie each mean explicitly to the *filename*, not the layer name
    fnames_noext <- tools::file_path_sans_ext(basename(tif_list))
    tab <- data.frame(
      rn = fnames_noext,   # we'll treat this as "layer name"
      V1 = as.numeric(means),
      stringsAsFactors = FALSE
    )
    setDT(tab, keep.rownames = FALSE)[]
    
    # remove excluded groups (based on filename)
    tab <- tab[!grepl("ART", tab$rn, ignore.case = TRUE)]
    tab <- tab[!grepl("HAR", tab$rn, ignore.case = TRUE)]
    tab <- tab[!grepl("LOG", tab$rn, ignore.case = TRUE)]
    
    # map rn (filename) -> Species
    tab$Species <- vapply(tab$rn, map_species, character(1))
    # drop rows that we couldn't map
    tab <- tab[!is.na(tab$Species), , drop = FALSE]
    
    if (!nrow(tab)){
      message(sprintf("[warn] No mappable species in %s", fp))
      next
    }
    
    # one row per species: if duplicates exist, keep the first (or average if you prefer)
    tab <- tab[order(tab$Species), ]
    tab <- tab[!duplicated(tab$Species), ]
    
    # site code
    nb <- stringi::stri_extract_last(basename(fp), regex = "(\\d+)")
    tab$Code   <- nb
    tab$lookup <- paste0(tab$Species, tab$Code)
    
    out_df <- tab |>
      dplyr::select(lookup, Species, Code, V1) |>
      dplyr::rename(!!col_TPI := V1)   # mean TPI
    
    site_name <- basename(fp)
    out_path  <- file.path(out_dir_res, paste0(site_name, "_", col_TPI, ".csv"))
    write.csv(out_df, out_path, row.names = FALSE)
    site_csvs <- c(site_csvs, out_path)
    message(sprintf("[ok] %s", out_path))
  }
  
  if (length(site_csvs)){
    master_res <- site_csvs %>% 
      lapply(readr::read_csv, show_col_types = FALSE) %>% 
      dplyr::bind_rows()
    
    write.csv(master_res, file.path(out_dir_res, paste0("Master_", col_TPI, ".csv")), row.names = FALSE)
    
    # Merge into global master
    if (file.exists(GLOBAL_MASTER_PATH)){
      global_master <- read.csv(GLOBAL_MASTER_PATH, check.names = FALSE)
    } else {
      write.csv(master_res, GLOBAL_MASTER_PATH, row.names = FALSE)
      message(sprintf("[init] Created %s with %s", GLOBAL_MASTER_PATH, col_TPI))
      return(invisible(NULL))
    }
    
    drop_cols <- c("X", "...1")
    look <- subset(master_res, select = setdiff(names(master_res), intersect(drop_cols, names(master_res))))
    
    merged <- merge(global_master, look, by = c("lookup","Species","Code"), all = TRUE) |>
      dplyr::arrange(lookup, Species, Code) |>
      dplyr::distinct(lookup, Species, Code, .keep_all = TRUE)
    
    write.csv(merged, GLOBAL_MASTER_PATH, row.names = FALSE)
    message(sprintf("[merge] Updated %s with %s", basename(GLOBAL_MASTER_PATH), col_TPI))
  }
}

# -------------------------
# Run over all resolutions (only parsable TPI*)
# -------------------------
res_dirs_all <- list.dirs(BASE_RES_DIR, full.names = TRUE, recursive = FALSE)
cand <- res_dirs_all[grepl("^TPI", basename(res_dirs_all), ignore.case = TRUE)]

can_parse <- function(name){
  tryCatch({ resolution_token(name); TRUE }, error = function(e) FALSE)
}
res_dirs <- cand[sapply(basename(cand), can_parse)]

if (!length(res_dirs)) stop("No parsable TPI* resolution folders found in BASE_RES_DIR")
for (rd in res_dirs) process_resolution_dir(rd)
message("[done] All resolutions processed.")

#------------------------------------------------------------------------------#
#---------------- REORGANIZE OUTPUT CSV NUMERICALLY ---------------------------#
#------------------------------------------------------------------------------#

# ==== Reorder TPI columns in Resample_TPI_mean.csv by increasing resolution ====
library(dplyr)
library(readr)

# --- Config ---
master_path <- GLOBAL_MASTER_PATH  # "E:/PhD/Chap 4 - Resolution/Resample_TPI_mean.csv"

# --- Load ---
df <- read_csv(master_path, show_col_types = FALSE)

# --- Identify and sort TPI columns numerically ---
tpi_cols <- grep("^TPI", names(df), value = TRUE)

# Extract numeric part (convert cm to meters where needed, if you ever use cm)
extract_res <- function(colname){
  val <- sub("^TPI", "", colname)
  val <- gsub(",", ".", val, fixed = TRUE)
  if (grepl("cm$", val, ignore.case = TRUE)) {
    val <- as.numeric(sub("cm$", "", val, ignore.case = TRUE)) / 100
  } else {
    val <- as.numeric(sub("m$", "", val, ignore.case = TRUE))
  }
  val
}

tpi_order <- tpi_cols[order(sapply(tpi_cols, extract_res))]

# --- Reorder columns: first three ("lookup","Species","Code") stay first ---
ordered_cols <- c("lookup", "Species", "Code", tpi_order)

# Keep only columns that actually exist in df
ordered_cols <- intersect(ordered_cols, names(df))

df_out <- df[, ordered_cols]

# --- Write back to file (overwriting existing) ---
write_csv(df_out, master_path)
message("[ok] Reordered TPI columns written to ", master_path)

# --- Optional: preview new column order ---
print(names(df_out))


################################################################################
#------ Get species averaged patch TPI ----------#
################################################################################

library(dplyr)
library(readr)

# --- Paths ---
base_dir      <- "E:/PhD/Chap 4 - Resolution"
in_master_fp  <- file.path(base_dir, "Resample_TPI_mean.csv")
out_code_fp   <- file.path(base_dir, "Resample_TPI_mean_byCode.csv")

# --- Load master (species-level) table ---
df <- read_csv(in_master_fp, show_col_types = FALSE)

# --- Identify TPI columns ---
tpi_cols <- grep("^TPI", names(df), value = TRUE)

# --- Aggregate to Code level: mean across species for each resolution ---
df_code <- df %>%
  group_by(Code) %>%
  dplyr::summarise(
    across(all_of(tpi_cols), ~ mean(.x, na.rm = TRUE))
  ) %>%
  ungroup()

# (optional) sort by numeric Code if possible
suppressWarnings({
  code_num <- as.numeric(df_code$Code)
})
if (!all(is.na(code_num))) {
  df_code <- df_code[order(code_num), ]
} else {
  df_code <- df_code[order(df_code$Code), ]
}

# --- Write out ---
write_csv(df_code, out_code_fp)
message("[ok] Wrote Code-level mean TPI to: ", out_code_fp)


################################################################################
################################# OI ###########################################
################################################################################
#------------------------------------------------------------------------------#
#-------------------- CREATE OI BY SITE AND RESOLUTION CSV --------------------#
#------------------------------------------------------------------------------#
# --- Libraries (lean) ---
library(dplyr)
library(readr)
library(data.table)
library(raster)
library(stringi)
library(tools)

# -------------------------
# Config: input / output
# -------------------------
BASE_RES_DIR <- "E:/PhD/Rasters/MSc Rasters/TPI/RESOLUTION"   # contains TPI0.1, TPI2, TPI32, ...
BASE_OUT     <- "E:/PhD/Chap 4 - Resolution"                  # <- all CSVs go here now
GLOBAL_MASTER_PATH <- file.path(BASE_OUT, "Resample_TPI_OI.csv")
dir.create(BASE_OUT, showWarnings = FALSE, recursive = TRUE)

# -------------------------
# Helpers
# -------------------------
OI <- function(x){
  n_valid <- length(x) - sum(is.na(x))
  if (n_valid == 0) return(NA_real_)
  100 * sum(x >= 1, na.rm = TRUE) / n_valid
}

# Parse token after "TPI" (supports integers, decimals, optional trailing 'm')
# Examples: TPI6 -> "6", TPI0.05 -> "0.05", TPI2m -> "2"
resolution_token <- function(res_folder_name){
  nm <- gsub("^TPI[ _-]*", "", res_folder_name, ignore.case = TRUE)
  nm <- trimws(nm)
  nm <- gsub(",", ".", nm, fixed = TRUE)
  m  <- regexec("^(\\d+(?:\\.\\d+)?)(?:m)?$", nm, perl = TRUE)
  r  <- regmatches(nm, m)[[1]]
  if (length(r) == 0) stop(sprintf("Could not parse resolution from '%s'", res_folder_name))
  r[2]  # numeric string in meters
}

# Convert token to meters label used in column name (strip any trailing 'm')
#   "6" -> "6"; "0.05" -> "0.05"; "2" -> "2"
meters_label <- function(token){
  token <- trimws(token)
  token <- gsub(",", ".", token, fixed = TRUE)
  sub("m$", "", token, ignore.case = TRUE)
}

oi_colname <- function(m_label) paste0("OI", m_label)

species_vec <- c(
  "antiopa","aphrodite","archippus","atalanta","augustinus","baptisiae",
  "bellona","bimacula","cocyta","comyntas","cybele","cymela","dion",
  "egeremet","eurytheme","henrici","lineola","lucia","lucilius","lygdamus",
  "origenes","pegala","philodice","plexippus","polios","progne","rapae",
  "sassacus","selene","tharos","themistocles","tullia","vestris"
)

# -------------------------
# Process one resolution dir (OI)
# -------------------------
process_resolution_dir <- function(res_dir){
  res_name  <- basename(res_dir)
  token     <- resolution_token(res_name)      # "6", "0.05", "2", ...
  m_label   <- meters_label(token)             # same as token here
  col_OI    <- oi_colname(m_label)            # "OI6", "OI0.05", ...
  
  site_root <- file.path(res_dir, "SITE")
  if (!dir.exists(site_root)){ 
    message(sprintf("[skip] No SITE/ in %s", res_dir)); 
    return(invisible(NULL)) 
  }
  
  out_dir_res <- file.path(BASE_OUT, "OI", col_OI)
  dir.create(out_dir_res, showWarnings = FALSE, recursive = TRUE)
  
  site_folders <- list.dirs(site_root, full.names = TRUE, recursive = FALSE)
  if (!length(site_folders)){ 
    message(sprintf("[warn] No site folders in %s", site_root)); 
    return(invisible(NULL)) 
  }
  
  # helper to map a *filename* (no path, no extension) to a species
  # expects strings like:
  #   AT0.1_1_5_18_ANTIOPA
  #   TPI10_2_7_18_SELENE
  map_species <- function(fname_noext){
    x <- tools::file_path_sans_ext(basename(fname_noext))
    x <- sub("\\.\\d+$", "", x)  # drop any .1, .2 suffixes raster might add
    x <- gsub("^LUC_", "LUCILIUS_", x, ignore.case = TRUE)  # legacy special case
    
    parts <- strsplit(x, "_", fixed = TRUE)[[1]]
    if (length(parts) == 0) return(NA_character_)
    
    cand  <- parts[length(parts)]         # last token, e.g. "ANTIOPA"
    cand  <- gsub("^LUC$", "LUCILIUS", cand, ignore.case = TRUE)
    
    sv     <- species_vec
    sv_l   <- tolower(sv)
    cand_l <- tolower(cand)
    
    hit <- match(cand_l, sv_l)
    if (!is.na(hit)) return(sv[hit])
    
    # fallback: whole string is species
    x_l <- tolower(x)
    hit2 <- match(x_l, sv_l)
    if (!is.na(hit2)) return(sv[hit2])
    
    # fallback: species name is prefix of full string
    hit3 <- which(startsWith(x_l, sv_l))
    if (length(hit3) == 1) return(sv[hit3])
    
    NA_character_
  }
  
  site_csvs <- character(0)
  
  for (fp in site_folders){
    tif_list <- list.files(fp, pattern = "\\.tif$", full.names = TRUE)
    if (!length(tif_list)){ 
      message(sprintf("[warn] No .tif in %s", fp)); 
      next 
    }
    
    # stack the rasters and force layer names to be the filename stems
    stk <- try(stack(tif_list), silent = TRUE)
    if (inherits(stk, "try-error")){ 
      message(sprintf("[error] Could not stack rasters in %s", fp)); 
      next 
    }
    fname_stems <- tools::file_path_sans_ext(basename(tif_list))
    names(stk)  <- fname_stems
    
    dts <- as.data.frame(stk)
    
    # compute OI per layer/column
    vals <- vapply(dts, OI, numeric(1))
    
    tab <- data.frame(
      rn = names(vals),   # these are filename stems
      V1 = as.numeric(vals),
      stringsAsFactors = FALSE
    )
    setDT(tab, keep.rownames = FALSE)[]
    
    # remove excluded groups (based on filename)
    tab <- tab[!grepl("ART", tab$rn, ignore.case = TRUE)]
    tab <- tab[!grepl("HAR", tab$rn, ignore.case = TRUE)]
    tab <- tab[!grepl("LOG", tab$rn, ignore.case = TRUE)]
    
    # map rn (filename stem) -> Species
    tab$Species <- vapply(tab$rn, map_species, character(1))
    # drop rows that we couldn't map
    tab <- tab[!is.na(tab$Species), , drop = FALSE]
    
    if (!nrow(tab)){
      message(sprintf("[warn] No mappable species in %s", fp))
      next
    }
    
    # one row per species: if duplicates exist, keep the first
    tab <- tab[order(tab$Species), ]
    tab <- tab[!duplicated(tab$Species), ]
    
    # site code
    nb <- stringi::stri_extract_last(basename(fp), regex = "(\\d+)")
    tab$Code   <- nb
    tab$lookup <- paste0(tab$Species, tab$Code)
    
    out_df <- tab |>
      dplyr::select(lookup, Species, Code, V1) |>
      dplyr::rename(!!col_OI := V1)
    
    site_name <- basename(fp)
    out_path  <- file.path(out_dir_res, paste0(site_name, "_", col_OI, ".csv"))
    write.csv(out_df, out_path, row.names = FALSE)
    site_csvs <- c(site_csvs, out_path)
    message(sprintf("[ok] %s", out_path))
  }
  
  if (length(site_csvs)){
    master_res <- site_csvs %>% 
      lapply(readr::read_csv, show_col_types = FALSE) %>% 
      dplyr::bind_rows()
    
    write.csv(master_res, file.path(out_dir_res, paste0("Master_", col_OI, ".csv")), row.names = FALSE)
    
    # Merge into global master
    if (file.exists(GLOBAL_MASTER_PATH)){
      global_master <- read.csv(GLOBAL_MASTER_PATH, check.names = FALSE)
    } else {
      write.csv(master_res, GLOBAL_MASTER_PATH, row.names = FALSE)
      message(sprintf("[init] Created %s with %s", GLOBAL_MASTER_PATH, col_OI))
      return(invisible(NULL))
    }
    
    drop_cols <- c("X", "...1")
    look <- subset(master_res, select = setdiff(names(master_res), intersect(drop_cols, names(master_res))))
    
    merged <- merge(global_master, look, by = c("lookup","Species","Code"), all = TRUE) |>
      dplyr::arrange(lookup, Species, Code) |>
      dplyr::distinct(lookup, Species, Code, .keep_all = TRUE)
    
    write.csv(merged, GLOBAL_MASTER_PATH, row.names = FALSE)
    message(sprintf("[merge] Updated %s with %s", basename(GLOBAL_MASTER_PATH), col_OI))
  }
}

# -------------------------
# Run over all resolutions (only parsable TPI*)
# -------------------------
res_dirs_all <- list.dirs(BASE_RES_DIR, full.names = TRUE, recursive = FALSE)
#res_dirs_all <- res_dirs_all[c(6,14,18,24)]
cand <- res_dirs_all[grepl("^TPI", basename(res_dirs_all), ignore.case = TRUE)]

can_parse <- function(name){
  tryCatch({ resolution_token(name); TRUE }, error = function(e) FALSE)
}
res_dirs <- cand[sapply(basename(cand), can_parse)]

if (!length(res_dirs)) stop("No parsable TPI* resolution folders found in BASE_RES_DIR")
for (rd in res_dirs) process_resolution_dir(rd)
message("[done] All resolutions processed.")

#------------------------------------------------------------------------------#
#---------------- REORGANIZE OUTPUT CSV NUMERICALLY ---------------------------#
#------------------------------------------------------------------------------#

# ==== Reorder OI columns in Resample_TPI_OI.csv by increasing resolution ====
library(dplyr)
library(readr)

# --- Config ---
master_path <- GLOBAL_MASTER_PATH  # use the same path as above

# --- Load ---
df <- read_csv(master_path, show_col_types = FALSE)

# --- Identify and sort OI columns numerically ---
oi_cols <- grep("^OI", names(df), value = TRUE)

# Extract numeric part (convert cm to meters where needed)
extract_res <- function(colname){
  val <- sub("^OI", "", colname)
  val <- gsub(",", ".", val, fixed = TRUE)
  if (grepl("cm$", val, ignore.case = TRUE)) {
    val <- as.numeric(sub("cm$", "", val, ignore.case = TRUE)) / 100
  } else {
    val <- as.numeric(sub("m$", "", val, ignore.case = TRUE))
  }
  val
}

oi_order <- oi_cols[order(sapply(oi_cols, extract_res))]

# --- Reorder columns: first three ("lookup","Species","Code") stay first ---
ordered_cols <- c("lookup", "Species", "Code", oi_order)

# Keep only columns that actually exist in df
ordered_cols <- intersect(ordered_cols, names(df))

df_out <- df[, ordered_cols]

# --- Write back to file (overwriting existing) ---
write_csv(df_out, master_path)
message("[ok] Reordered OI columns written to ", master_path)

# --- Optional: preview new column order ---
print(names(df_out))

#------------------------------------------------------------------------------#
#--------------------- Merge with master file ---------------------------------#
#------------------------------------------------------------------------------#

library(dplyr)

# --- File paths ---
base_csv  <- "D:/Gabby/MSc/Data/TPI_OI_Limestone.csv"
merge_csv <- GLOBAL_MASTER_PATH
out_csv   <- file.path(BASE_OUT, "TPI_OI_Resample_final.csv")

# --- Read data ---
df_base  <- read.csv(base_csv, stringsAsFactors = FALSE)
df_merge <- read.csv(merge_csv, stringsAsFactors = FALSE)

# --- Identify columns to add (exclude shared except lookup) ---
shared_cols <- c("Code", "Species", "lookup")
cols_to_add <- setdiff(names(df_merge), shared_cols)

# --- Merge on 'lookup' only ---
df_out2 <- df_base %>%
  left_join(df_merge[, c("lookup", cols_to_add)], by = "lookup")

# --- Write result ---
write.csv(df_out2, out_csv, row.names = FALSE)

cat("Merged CSV saved to:", out_csv, "\n")

################################################################################
################################################################################
################################################################################
#                             Results                                          #
################################################################################
################################################################################
################################################################################

#------------------------------------------------------------------------------#
#------------- Build model comparisons ----------------------------------------#
#------------------------------------------------------------------------------#

# ============================================================
# Compare model performance across OI (thermal) resolutions
# ============================================================

# --- Packages ---
library(glmmTMB)
library(lme4)
library(performance)   # for r2()
library(dplyr)
library(stringr)
library(purrr)
library(readr)

# --- Load merged dataset ---
tpioi <- read.csv("E:/PhD/Chap 4 - Resolution/TPI_OI_Resample_final.csv")

# --- Identify OI columns (start with "OI" but not "OI_Site") ---
oi_cols <- names(tpioi)
oi_cols <- oi_cols[grepl("^OI", oi_cols) & oi_cols != "OI_Site"]
if (length(oi_cols) == 0L) stop("No OI resolution columns found (e.g., OI0.15, OI1, OI10, ...).")

# --- Helper: extract numeric resolution from column name (e.g., 'OI0.15' -> 0.15) ---
parse_resolution <- function(x) as.numeric(sub("^OI", "", x))

# --- Helper: extract p-value for a fixed effect term ---
get_pval <- function(model, term, model_type = c("glmmTMB","glmer")) {
  model_type <- match.arg(model_type)
  if (model_type == "glmmTMB") {
    co <- summary(model)$coefficients$cond
  } else {
    co <- summary(model)$coefficients
  }
  if (!term %in% rownames(co)) return(NA_real_)
  pcol <- grep("Pr\\(>|[tz]\\|\\)", colnames(co), value = TRUE)
  if (length(pcol) == 0L) return(NA_real_)
  as.numeric(co[term, pcol[1]])
}

# --- Fit abundance model for one OI column ---
fit_abundance_one <- function(oi_col) {
  req <- c("Abundance","FHD_Simp","Floral_Richness","Size","Site","Species",oi_col)
  dat <- tpioi %>% dplyr::select(all_of(req)) %>% na.omit()
  if (nrow(dat) < 10) return(NULL)
  
  fml <- as.formula(
    paste0("Abundance ~ ", oi_col, " + FHD_Simp + Floral_Richness + Size + (1|Site) + (1|Species)")
  )
  
  mdl <- glmmTMB(fml, family = nbinom2, data = dat)
  r2v <- performance::r2(mdl)
  
  tibble(
    Model          = "abundance",
    OI_column      = oi_col,
    Resolution     = parse_resolution(oi_col),
    AIC            = AIC(mdl),
    OI_p_value     = get_pval(mdl, oi_col, "glmmTMB"),
    R2_marginal    = as.numeric(r2v$R2_marginal),
    R2_conditional = as.numeric(r2v$R2_conditional),
    N_used         = nrow(dat)
  )
}

# --- Fit occurrence model for one OI column ---
fit_occurrence_one <- function(oi_col) {
  req <- c("Occurrence","FHD_Simp","Floral_Richness","Size","Site","Species",oi_col)
  dat <- tpioi %>% dplyr::select(all_of(req)) %>% na.omit()
  if (nrow(dat) < 10) return(NULL)
  
  fml <- as.formula(
    paste0("Occurrence ~ ", oi_col, " + FHD_Simp + Floral_Richness + Size + (1|Species) + (1|Site)")
  )
  
  mdl <- glmer(fml, family = binomial(link="logit"), data = dat,
               control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
  r2v <- performance::r2(mdl)
  
  tibble(
    Model          = "occurrence",
    OI_column      = oi_col,
    Resolution     = parse_resolution(oi_col),
    AIC            = AIC(mdl),
    OI_p_value     = get_pval(mdl, oi_col, "glmer"),
    R2_marginal    = as.numeric(r2v$R2_marginal),
    R2_conditional = as.numeric(r2v$R2_conditional),
    N_used         = nrow(dat)
  )
}

# --- Run models across all OI resolutions ---
res_abund <- map_dfr(oi_cols, ~ tryCatch(fit_abundance_one(.x), error = function(e) NULL))
res_occ   <- map_dfr(oi_cols, ~ tryCatch(fit_occurrence_one(.x), error = function(e) NULL))
results   <- bind_rows(res_abund, res_occ) %>% arrange(Model, Resolution)

# --- Save results to CSV ---
out_csv <- "E:/PhD/Chap 4 - Resolution/OI_resolution_model_performance.csv"
readr::write_csv(results, out_csv)
message("Saved results to: ", out_csv)

# --- Optional: quick preview ---
print(dplyr::glimpse(results))


#------------------------------------------------------------------------------#
#------------ Same for TPI ----------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#      TPI: Track p-values & effect sizes across resolutions for              #
#      richness (Richness) and diversity (Hill_q1) models                     #
#------------------------------------------------------------------------------#

# --- Packages ---
library(lme4)
library(performance)   # for r2()
library(dplyr)
library(stringr)
library(purrr)
library(readr)
library(tibble)

# --- Load survey-level dataset (formerly `hill`) ---
hill <- read.csv("E:/PhD/Chap 4 - Resolution/MASTER.csv", stringsAsFactors = FALSE)

# --- Identify TPI columns: TPIx where x is numeric (e.g., TPI0.05, TPI1, TPI52) ---
all_tpi <- grep("^TPI", names(hill), value = TRUE)

# keep only TPI columns whose suffix is numeric (to avoid TPI_S_New, TPI_Site, etc.)
suffix_num <- suppressWarnings(as.numeric(sub("^TPI", "", all_tpi)))
tpi_cols   <- all_tpi[!is.na(suffix_num)]

# In case TPI_S_New happens to match the pattern, explicitly drop it
tpi_cols <- setdiff(tpi_cols, "TPI_S_New")

if (length(tpi_cols) == 0L) {
  stop("No TPI resolution columns found (e.g., TPI0.05, TPI1, TPI10, ...).")
}

# --- Helper: extract numeric resolution from column name (e.g., 'TPI0.15' -> 0.15) ---
parse_resolution <- function(x) as.numeric(sub("^TPI", "", x))

# --- Helper: extract approximate p-value for a fixed-effect term from an lmer model ---
# lmer doesn't provide p-values by default, so we approximate using the t-statistic
get_pval_lmer <- function(model, term) {
  co <- summary(model)$coefficients
  if (!term %in% rownames(co)) return(NA_real_)
  tval <- co[term, "t value"]
  # Approximate two-sided p using normal distribution
  2 * pnorm(abs(tval), lower.tail = FALSE)
}

#------------------------------------------------------------------------------#
# 1) Diversity model: Hill_q1 ~ TPI_S_New + FHD_Simp + Floral_Richness + Size
#------------------------------------------------------------------------------#
fit_diversity_one <- function(tpi_col) {
  # Build dataset for this resolution:
  # rename the current TPI column as TPI_S_New so we can reuse the same formula
  req <- c("Hill_q1", "FHD_Simp", "Floral_Richness", "Size", "Site", tpi_col)
  dat <- hill %>%
    dplyr::select(all_of(req)) %>%
    dplyr::rename(TPI_S_New = all_of(tpi_col)) %>%
    na.omit()
  
  if (nrow(dat) < 10) return(NULL)
  
  fml <- Hill_q1 ~ TPI_S_New + FHD_Simp + Floral_Richness + Size + (1 | Site)
  
  mdl <- lmer(fml, data = dat, REML = FALSE)
  r2v <- performance::r2(mdl)
  co  <- summary(mdl)$coefficients
  
  est <- co["TPI_S_New", "Estimate"]
  se  <- co["TPI_S_New", "Std. Error"]
  tval <- co["TPI_S_New", "t value"]
  pval <- get_pval_lmer(mdl, "TPI_S_New")
  
  tibble(
    Model           = "diversity",
    TPI_column      = tpi_col,
    Resolution      = parse_resolution(tpi_col),
    AIC             = AIC(mdl),
    TPI_estimate    = est,
    TPI_se          = se,
    TPI_t_value     = tval,
    TPI_p_value     = pval,
    TPI_abs_t       = abs(tval),                 # effect-size proxy
    R2_marginal     = as.numeric(r2v$R2_marginal),
    R2_conditional  = as.numeric(r2v$R2_conditional),
    N_used          = nrow(dat)
  )
}

#------------------------------------------------------------------------------#
# 2) Richness model: Richness ~ TPI_S_New + FHD_Simp + Floral_Richness + Size
#------------------------------------------------------------------------------#
fit_richness_one <- function(tpi_col) {
  req <- c("Richness", "FHD_Simp", "Floral_Richness", "Size", "Site", tpi_col)
  dat <- hill %>%
    dplyr::select(all_of(req)) %>%
    dplyr::rename(TPI_S_New = all_of(tpi_col)) %>%
    na.omit()
  
  if (nrow(dat) < 10) return(NULL)
  
  fml <- Richness ~ TPI_S_New + FHD_Simp + Floral_Richness + Size + (1 | Site)
  
  mdl <- lmer(fml, data = dat, REML = FALSE)
  r2v <- performance::r2(mdl)
  co  <- summary(mdl)$coefficients
  
  est  <- co["TPI_S_New", "Estimate"]
  se   <- co["TPI_S_New", "Std. Error"]
  tval <- co["TPI_S_New", "t value"]
  pval <- get_pval_lmer(mdl, "TPI_S_New")
  
  tibble(
    Model           = "richness",
    TPI_column      = tpi_col,
    Resolution      = parse_resolution(tpi_col),
    AIC             = AIC(mdl),
    TPI_estimate    = est,
    TPI_se          = se,
    TPI_t_value     = tval,
    TPI_p_value     = pval,
    TPI_abs_t       = abs(tval),                 # effect-size proxy
    R2_marginal     = as.numeric(r2v$R2_marginal),
    R2_conditional  = as.numeric(r2v$R2_conditional),
    N_used          = nrow(dat)
  )
}

#------------------------------------------------------------------------------#
# 3) Run models across all TPI resolutions for both responses
#------------------------------------------------------------------------------#

# Diversity (Hill_q1)
res_div <- purrr::map_dfr(
  tpi_cols,
  ~ tryCatch(fit_diversity_one(.x), error = function(e) NULL)
)

# Richness
res_rich <- purrr::map_dfr(
  tpi_cols,
  ~ tryCatch(fit_richness_one(.x), error = function(e) NULL)
)

# Combine and order by model and resolution
tpi_results <- bind_rows(res_div, res_rich) %>%
  arrange(Model, Resolution)

#------------------------------------------------------------------------------#
# 4) Save results to CSV
#------------------------------------------------------------------------------#
out_csv <- "E:/PhD/Chap 4 - Resolution/TPI_resolution_model_performance.csv"
readr::write_csv(tpi_results, out_csv)
message("Saved TPI results to: ", out_csv)

# Optional quick preview
print(dplyr::glimpse(tpi_results))

#------------------------------------------------------------------------------#
#---- Max/Min Air Temperature per Site-Date-Resolution (with Code) -----------#
#------------------------------------------------------------------------------#

library(raster)
library(dplyr)

# Root directory containing subfolders like AT0.1, AT1, AT2, ..., AT60, etc.
base_at_dir <- "E:/PhD/Rasters/MSc Rasters/Air Temperature/RESOLUTION"

# Output CSV path
out_csv <- "E:/PhD/Chap 4 - Resolution/AirTemp_MinMax_bySiteResolution.csv"

#-------------------------------------------------------------------------------
# 1) Find only resolution folders whose names start with "AT"
#-------------------------------------------------------------------------------

# First-level subdirectories under base_at_dir
res_dirs <- list.dirs(base_at_dir, full.names = TRUE, recursive = FALSE)

# Keep only those whose basename starts with "AT"
res_dirs <- res_dirs[grepl("^AT", basename(res_dirs), ignore.case = FALSE)]

message("Resolution folders detected: ", paste(basename(res_dirs), collapse = ", "))

# List all .tif files inside those AT* folders
at_files <- unlist(lapply(res_dirs, function(d) {
  list.files(d, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)
}))

message("Total .tif files found in AT* folders: ", length(at_files))

if (!length(at_files)) {
  stop("No .tif files found inside AT* folders under: ", base_at_dir)
}

#-------------------------------------------------------------------------------
# 2) Parse metadata from filenames: AT<resolution>_<site>_<month>_<day>.tif
#-------------------------------------------------------------------------------

parse_meta <- function(fp) {
  bn <- basename(fp)
  m  <- regexec("^AT([^_]+)_(\\d+)_(\\d+)_(\\d+)\\.tif$", bn)
  mm <- regmatches(bn, m)[[1]]
  if (length(mm) != 5) {
    warning("Filename did not match expected pattern ATres_site_month_day.tif: ", bn)
    return(NULL)
  }
  data.frame(
    file       = fp,
    Resolution = mm[2],                  # e.g., "1", "0.05", "60"
    Site       = as.integer(mm[3]),      # site ID
    Month      = as.integer(mm[4]),      # month (1–12)
    Day        = as.integer(mm[5]),      # day (1–31)
    stringsAsFactors = FALSE
  )
}

meta_list <- lapply(at_files, parse_meta)
meta <- dplyr::bind_rows(meta_list)

if (nrow(meta) == 0) {
  stop("No AT*.tif files matched expected pattern 'ATres_site_month_day.tif'.")
}

#-------------------------------------------------------------------------------
# 3) Build survey-level Code:
#    - A "survey" is a unique (Site, Month, Day) combo.
#    - For each Site, surveys are ordered chronologically.
#      Site 3: first survey -> Code 31, second -> 32, ...
#      Site 12: first survey -> 121, second -> 122, ...
#-------------------------------------------------------------------------------

survey_codes <- meta %>%
  dplyr::distinct(Site, Month, Day) %>%
  dplyr::arrange(Site, Month, Day) %>%
  dplyr::group_by(Site) %>%
  dplyr::mutate(
    Visit = dplyr::row_number(),         # chronological survey index within site
    Code  = paste0(Site, Visit)          # e.g., Site 3 Visit 1 -> "31"; Site 12 Visit 1 -> "121"
  ) %>%
  dplyr::ungroup()

# Join the survey-level Code back onto every raster (all resolutions) for that survey
meta <- meta %>%
  dplyr::left_join(
    survey_codes %>% dplyr::select(Site, Month, Day, Visit, Code),
    by = c("Site", "Month", "Day")
  )

#-------------------------------------------------------------------------------
# 4) Compute Min and Max for each raster (handle all-NA rasters gracefully)
#-------------------------------------------------------------------------------

meta$Min <- NA_real_
meta$Max <- NA_real_

for (i in seq_len(nrow(meta))) {
  r <- raster(meta$file[i])
  vals <- raster::getValues(r)
  vals <- vals[is.finite(vals)]  # drop NA/NaN/Inf
  
  if (!length(vals)) {
    # all values NA / non-finite
    meta$Min[i] <- NA_real_
    meta$Max[i] <- NA_real_
  } else {
    meta$Min[i] <- min(vals, na.rm = TRUE)
    meta$Max[i] <- max(vals, na.rm = TRUE)
  }
}

#-------------------------------------------------------------------------------
# 5) Final table: one row per raster (site-date-resolution), with survey Code
#-------------------------------------------------------------------------------

airtemp_mm <- meta %>%
  dplyr::select(Code, Site, Month, Day, Resolution, Max, Min) %>%
  dplyr::arrange(Site, Code, Resolution)

# Write CSV
write.csv(airtemp_mm, out_csv, row.names = FALSE)
message("✅ Saved air temperature min/max summary to: ", out_csv)

################################################################################
########################## PLOTTING ############################################
################################################################################

#------------------------------------------------------------------------------#
#  Plotting Air Temperature Min/Max vs Resolution (≤ 100 m only)               #
#  + Fix duplicate 15 m rows (keep lower-temp entry per site)                  #
#  + No grid, no box; only axis lines                                         #
#------------------------------------------------------------------------------#

library(dplyr)
library(tidyr)
library(ggplot2)

# Path to your CSV
csv_path <- "E:/PhD/Chap 4 - Resolution/AirTemp_MinMax_bySiteResolution.csv"

# Read data
airtemp_mm <- read.csv(csv_path, stringsAsFactors = FALSE)

# Ensure Resolution is numeric and keep <= 100 m
airtemp_mm <- airtemp_mm %>%
  mutate(Resolution = as.numeric(Resolution)) %>%
  filter(!is.na(Resolution) & Resolution <= 100)

# Long format for Min/Max points
df_long <- airtemp_mm %>%
  tidyr::pivot_longer(
    cols = c(Min, Max),
    names_to = "Type",
    values_to = "Temp"
  )

# Per-resolution global min & max (across all sites/dates)
# Force use of dplyr versions to avoid plyr conflicts
df_range <- airtemp_mm %>%
  dplyr::group_by(Resolution) %>%
  dplyr::summarise(
    ymin = min(Min, na.rm = TRUE),
    ymax = max(Max, na.rm = TRUE),
    .groups = "drop"
  )

# (Optional) sanity check
# print(names(df_range))
# print(head(df_range))

# Plot: points + "errorbar-style" range line for each resolution
p1 <- ggplot() +
  geom_point(
    data = df_long,
    aes(x = Resolution, y = Temp, color = Type),
    position = position_jitter(width = 0.1, height = 0),
    alpha = 0.6,
    size = 1.8
  ) +
  geom_errorbar(
    data = df_range,
    aes(x = Resolution, ymin = ymin, ymax = ymax),
    width = 0.15,          # cap width
    color = "black",
    linewidth = 0.6
  ) +
  labs(
    x = "Resolution (m)",
    y = "Air Temperature (°C)",
    color = "Value type",
    title = "Site-level Min/Max Air Temperature by Resolution (≤ 100 m)"
  ) +
  theme_classic() +
  theme(
    plot.title      = element_text(hjust = 0.5),
    legend.position = "top"
  )

print(p1)


#-----------------------------------------------
# 2) Plot: floating bars of mean Min/Max per resolution
#-----------------------------------------------

library(dplyr)
library(ggplot2)

# Make sure Resolution is numeric and limited to <= 100 m
airtemp_mm <- airtemp_mm %>%
  mutate(Resolution = as.numeric(Resolution)) %>%
  filter(!is.na(Resolution) & Resolution <= 100)

# Per-resolution mean of Min and Max
df_avg <- airtemp_mm %>%
  dplyr::group_by(Resolution) %>%
  dplyr::summarise(
    Mean_Min = mean(Min, na.rm = TRUE),
    Mean_Max = mean(Max, na.rm = TRUE),
    .groups  = "drop"
  )

# Plot 2: "floating bar" (errorbar) from mean min to mean max
p2 <- ggplot(df_avg, aes(x = Resolution)) +
  geom_errorbar(
    aes(ymin = Mean_Min, ymax = Mean_Max),
    width     = 0.3,   # cap width
    linewidth = 0.8
  ) +
  labs(
    x     = "Resolution (m)",
    y     = "Air Temperature (°C)",
    title = "Mean Min/Max Air Temperature by Resolution"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

print(p2)


#------------------------------------------------------------------------------#
#  Do Max and Min Air Temperatures Decrease Significantly with Resolution?     #
#------------------------------------------------------------------------------#

library(dplyr)

# Path to your CSV from previous steps
csv_path <- "E:/PhD/Chap 4 - Resolution/AirTemp_MinMax_bySiteResolution.csv"

# Read data
airtemp_mm <- read.csv(csv_path, stringsAsFactors = FALSE)

# Ensure Resolution is numeric
airtemp_mm <- airtemp_mm %>%
  mutate(Resolution = as.numeric(Resolution))

#------------------------------------------------------------------------------#
# Optional: restrict to resolutions ≤ 100 m (to match your plots)              #
#------------------------------------------------------------------------------#

airtemp_mm <- airtemp_mm %>%
  filter(!is.na(Resolution) & Resolution <= 100)

#------------------------------------------------------------------------------#
# 1) Simple linear models: Min ~ Resolution and Max ~ Resolution               #
#------------------------------------------------------------------------------#

# Min vs Resolution
lm_min  <- lm(Min ~ Resolution, data = airtemp_mm)
cat("\n================ Min ~ Resolution ================\n")
print(summary(lm_min))

# Max vs Resolution
lm_max  <- lm(Max ~ Resolution, data = airtemp_mm)
cat("\n================ Max ~ Resolution ================\n")
print(summary(lm_max))

# Interpretation:
#  - Check the coefficient for 'Resolution' in each summary:
#      * Estimate < 0 and p-value < 0.05  -> significant DECREASE with resolution.
#      * Estimate > 0 and p-value < 0.05  -> significant INCREASE with resolution.
#      * p-value >= 0.05                  -> no strong evidence for a linear trend.

library(dplyr)

#------------------------------------------------------------------------------#
# 2) Site-level models (average per Site–Resolution before testing)            #
#------------------------------------------------------------------------------#

# Make sure Resolution is numeric and not NA, then collapse to site–resolution means
airtemp_site <- airtemp_mm %>%
  dplyr::mutate(
    Resolution = as.numeric(Resolution)
  ) %>%
  dplyr::filter(!is.na(Resolution)) %>%
  dplyr::group_by(Site, Resolution) %>%
  dplyr::summarise(
    Min = mean(Min, na.rm = TRUE),
    Max = mean(Max, na.rm = TRUE),
    .groups = "drop"
  )

# Quick sanity check (optional)
print(head(airtemp_site))
print(str(airtemp_site))

# Linear models: site-mean Min/Max ~ Resolution
lm_min_site <- lm(Min ~ Resolution, data = airtemp_site)
lm_max_site <- lm(Max ~ Resolution, data = airtemp_site)

cat("\n========== Site-mean Min ~ Resolution ==========\n")
print(summary(lm_min_site))

cat("\n========== Site-mean Max ~ Resolution ==========\n")
print(summary(lm_max_site))

# Same interpretation: look at the 'Resolution' coefficient sign and p-value
# in each model summary to decide whether Min and Max decrease significantly
# with increasing resolution.

#==============================================================================#
#  Linear models: log(p-value) ~ Resolution for all models                     #
#    Models: abundance, occurrence (OI); richness, diversity (TPI)             #
#    Output: CSV with summary-style stats for each model                       #
#==============================================================================#

library(dplyr)
library(readr)

#-----------------------------#
# Paths                                                           
#-----------------------------#
oi_perf_path  <- "E:/PhD/Chap 4 - Resolution/OI_resolution_model_performance.csv"
tpi_perf_path <- "E:/PhD/Chap 4 - Resolution/TPI_resolution_model_performance.csv"
out_path      <- "E:/PhD/Chap 4 - Resolution/logp_vs_resolution_model_summaries.csv"

#-----------------------------#
# 1) Load and harmonize p-values and resolutions                  
#-----------------------------#

# OI-based models: abundance & occurrence
oi_df <- read_csv(oi_perf_path, show_col_types = FALSE) %>%
  mutate(
    Source     = "OI",
    Model      = tolower(Model),          # "abundance", "occurrence"
    Resolution = as.numeric(Resolution),
    p_value    = as.numeric(OI_p_value)
  ) %>%
  filter(!is.na(Resolution), !is.na(p_value)) %>%
  dplyr::select(Source, Model, Resolution, p_value)

# TPI-based models: richness & diversity
tpi_df <- read_csv(tpi_perf_path, show_col_types = FALSE) %>%
  mutate(
    Source     = "TPI",
    Model      = tolower(Model),          # "richness", "diversity"
    Resolution = as.numeric(Resolution),
    p_value    = as.numeric(TPI_p_value)
  ) %>%
  filter(!is.na(Resolution), !is.na(p_value)) %>%
  dplyr::select(Source, Model, Resolution, p_value)

# Combine all
p_all <- bind_rows(oi_df, tpi_df) %>%
  arrange(Source, Model, Resolution)

#-----------------------------#
# 2) Helper to fit lm(logp ~ Resolution) and extract summary-like stats       
#-----------------------------#

fit_logp_model <- function(df_sub, model_name, source_name) {
  if (nrow(df_sub) < 3) return(NULL)  # not enough points to fit
  
  # Use -log10(p) so higher values = stronger signal
  df_sub <- df_sub %>%
    mutate(logp = -log10(p_value + 1e-12))
  
  m <- lm(logp ~ Resolution, data = df_sub)
  s <- summary(m)
  
  # Coefficients
  coefs <- s$coefficients
  # rows: "(Intercept)", "Resolution"
  # cols: Estimate, Std. Error, t value, Pr(>|t|)
  
  # F-statistic info
  fstat <- unname(s$fstatistic)  # c(F, df1, df2)
  F_value <- fstat[1]
  df1     <- fstat[2]
  df2     <- fstat[3]
  F_p     <- pf(F_value, df1, df2, lower.tail = FALSE)
  
  tibble(
    Source          = source_name,
    Model           = model_name,
    N_points        = nrow(df_sub),
    min_resolution  = min(df_sub$Resolution, na.rm = TRUE),
    max_resolution  = max(df_sub$Resolution, na.rm = TRUE),
    
    # Intercept
    Intercept_est   = coefs["(Intercept)", "Estimate"],
    Intercept_se    = coefs["(Intercept)", "Std. Error"],
    Intercept_t     = coefs["(Intercept)", "t value"],
    Intercept_p     = coefs["(Intercept)", "Pr(>|t|)"],
    
    # Slope for Resolution
    Slope_est       = coefs["Resolution", "Estimate"],
    Slope_se        = coefs["Resolution", "Std. Error"],
    Slope_t         = coefs["Resolution", "t value"],
    Slope_p         = coefs["Resolution", "Pr(>|t|)"],
    
    # Goodness of fit
    R2              = s$r.squared,
    R2_adj          = s$adj.r.squared,
    Residual_SE     = s$sigma,
    
    # F-statistic for overall regression
    F_stat          = F_value,
    F_df1           = df1,
    F_df2           = df2,
    F_p             = F_p
  )
}

#-----------------------------#
# 3) Fit models for each of the four cases                                    #
#-----------------------------#

models_to_run <- tibble::tribble(
  ~Source, ~Model,
  "OI",    "abundance",
  "OI",    "occurrence",
  "TPI",   "richness",
  "TPI",   "diversity"
)

summary_rows <- lapply(seq_len(nrow(models_to_run)), function(i) {
  src  <- models_to_run$Source[i]
  mdl  <- models_to_run$Model[i]
  
  df_sub <- p_all %>%
    dplyr::filter(Source == src, Model == mdl) %>%
    dplyr::arrange(Resolution)
  
  fit_logp_model(df_sub, mdl, src)
})

summary_table <- bind_rows(summary_rows)

#-----------------------------#
# 4) Save and inspect                                            #
#-----------------------------#
write_csv(summary_table, out_path)
cat("Saved log(p-value) ~ Resolution model summaries to:\n  ", out_path, "\n")

print(summary_table)


#==============================================================================#
#  AIC vs Resolution — wide-format output CSV                                  #
#    Columns: Resolution, abundance, occurrence, richness, diversity           #
#==============================================================================#

library(dplyr)
library(readr)
library(tidyr)

#-----------------------------#
# Paths                                                           
#-----------------------------#
oi_perf_path  <- "E:/PhD/Chap 4 - Resolution/OI_resolution_model_performance.csv"
tpi_perf_path <- "E:/PhD/Chap 4 - Resolution/TPI_resolution_model_performance.csv"
out_path      <- "E:/PhD/Chap 4 - Resolution/aic_vs_resolution.csv"

#-----------------------------#
# 1) Load AIC values per resolution
#    NOTE: replace "AIC" below with your actual AIC column name if different
#-----------------------------#

oi_df <- read_csv(oi_perf_path, show_col_types = FALSE) %>%
  mutate(
    Model      = tolower(Model),
    Resolution = as.numeric(Resolution),
    AIC        = as.numeric(AIC)          # <-- adjust column name if needed
  ) %>%
  filter(!is.na(Resolution), !is.na(AIC)) %>%
  dplyr::select(Resolution, Model, AIC)

tpi_df <- read_csv(tpi_perf_path, show_col_types = FALSE) %>%
  mutate(
    Model      = tolower(Model),
    Resolution = as.numeric(Resolution),
    AIC        = as.numeric(AIC)          # <-- adjust column name if needed
  ) %>%
  filter(!is.na(Resolution), !is.na(AIC)) %>%
  dplyr::select(Resolution, Model, AIC)

#-----------------------------#
# 2) Pivot each to wide, then join on Resolution
#-----------------------------#

oi_wide <- oi_df %>%
  pivot_wider(names_from = Model, values_from = AIC)   # cols: abundance, occurrence

tpi_wide <- tpi_df %>%
  pivot_wider(names_from = Model, values_from = AIC)   # cols: richness, diversity

aic_wide <- full_join(oi_wide, tpi_wide, by = "Resolution") %>%
  arrange(Resolution) %>%
  dplyr::select(Resolution, abundance, occurrence, richness, diversity)  # enforce column order

#-----------------------------#
# 3) Save and inspect
#-----------------------------#

write_csv(aic_wide, out_path)
cat("Saved AIC ~ Resolution wide table to:\n  ", out_path, "\n")
print(aic_wide)

#==============================================================================#
#  AIC change across resolutions — lm(AIC ~ Resolution) for all 4 models      #
#    Mirrors the logp ~ Resolution analysis above                              #
#    Output: CSV with regression stats for each model                         #
#==============================================================================#

library(dplyr)
library(readr)

#-----------------------------#
# Paths
#-----------------------------#
oi_perf_path  <- "E:/PhD/Chap 4 - Resolution/OI_resolution_model_performance.csv"
tpi_perf_path <- "E:/PhD/Chap 4 - Resolution/TPI_resolution_model_performance.csv"
out_path      <- "E:/PhD/Chap 4 - Resolution/aic_vs_resolution_model_summaries.csv"

#-----------------------------#
# 1) Load and harmonize AIC and resolutions
#-----------------------------#

oi_aic_df <- read_csv(oi_perf_path, show_col_types = FALSE) %>%
  mutate(
    Source     = "OI",
    Model      = tolower(Model),
    Resolution = as.numeric(Resolution),
    AIC        = as.numeric(AIC)
  ) %>%
  filter(!is.na(Resolution), !is.na(AIC)) %>%
  dplyr::select(Source, Model, Resolution, AIC)

tpi_aic_df <- read_csv(tpi_perf_path, show_col_types = FALSE) %>%
  mutate(
    Source     = "TPI",
    Model      = tolower(Model),
    Resolution = as.numeric(Resolution),
    AIC        = as.numeric(AIC)
  ) %>%
  filter(!is.na(Resolution), !is.na(AIC)) %>%
  dplyr::select(Source, Model, Resolution, AIC)

aic_all_lm <- bind_rows(oi_aic_df, tpi_aic_df) %>%
  arrange(Source, Model, Resolution)

#-----------------------------#
# 2) Helper: fit lm(AIC ~ Resolution) and extract summary stats
#-----------------------------#

fit_aic_lm <- function(df_sub, model_name, source_name) {
  if (nrow(df_sub) < 3) return(NULL)
  
  m <- lm(AIC ~ Resolution, data = df_sub)
  s <- summary(m)
  
  coefs   <- s$coefficients
  fstat   <- unname(s$fstatistic)
  F_value <- fstat[1]; df1 <- fstat[2]; df2 <- fstat[3]
  F_p     <- pf(F_value, df1, df2, lower.tail = FALSE)
  
  tibble(
    Source         = source_name,
    Model          = model_name,
    N_points       = nrow(df_sub),
    min_resolution = min(df_sub$Resolution, na.rm = TRUE),
    max_resolution = max(df_sub$Resolution, na.rm = TRUE),
    
    Intercept_est  = coefs["(Intercept)", "Estimate"],
    Intercept_se   = coefs["(Intercept)", "Std. Error"],
    Intercept_t    = coefs["(Intercept)", "t value"],
    Intercept_p    = coefs["(Intercept)", "Pr(>|t|)"],
    
    Slope_est      = coefs["Resolution", "Estimate"],
    Slope_se       = coefs["Resolution", "Std. Error"],
    Slope_t        = coefs["Resolution", "t value"],
    Slope_p        = coefs["Resolution", "Pr(>|t|)"],
    
    R2             = s$r.squared,
    R2_adj         = s$adj.r.squared,
    Residual_SE    = s$sigma,
    
    F_stat         = F_value,
    F_df1          = df1,
    F_df2          = df2,
    F_p            = F_p
  )
}

#-----------------------------#
# 3) Fit for all four models
#-----------------------------#

aic_models_to_run <- tibble::tribble(
  ~Source, ~Model,
  "OI",    "abundance",
  "OI",    "occurrence",
  "TPI",   "richness",
  "TPI",   "diversity"
)

aic_summary_rows <- lapply(seq_len(nrow(aic_models_to_run)), function(i) {
  src <- aic_models_to_run$Source[i]
  mdl <- aic_models_to_run$Model[i]
  
  df_sub <- aic_all_lm %>%
    dplyr::filter(Source == src, Model == mdl) %>%
    dplyr::arrange(Resolution)
  
  fit_aic_lm(df_sub, mdl, src)
})

aic_summary_table <- bind_rows(aic_summary_rows)

#-----------------------------#
# 4) Save and inspect
#-----------------------------#

write_csv(aic_summary_table, out_path)
cat("Saved AIC ~ Resolution model summaries to:\n  ", out_path, "\n")
print(aic_summary_table)

#==============================================================================#
#  Delta-AIC relative to 0.05m resolution                                     #
#==============================================================================#

library(dplyr)
library(readr)
library(ggplot2)

#-----------------------------#
# Paths                                                           
#-----------------------------#
oi_perf_path  <- "E:/PhD/Chap 4 - Resolution/OI_resolution_model_performance.csv"
tpi_perf_path <- "E:/PhD/Chap 4 - Resolution/TPI_resolution_model_performance.csv"

#-----------------------------#
# 1) Load data
#-----------------------------#

oi_df <- read_csv(oi_perf_path, show_col_types = FALSE) %>%
  mutate(
    Source     = "OI",
    Model      = tolower(Model),
    Resolution = as.numeric(Resolution),
    AIC        = as.numeric(AIC)
  ) %>%
  filter(!is.na(Resolution), !is.na(AIC)) %>%
  dplyr::select(Source, Model, Resolution, AIC)

tpi_df <- read_csv(tpi_perf_path, show_col_types = FALSE) %>%
  mutate(
    Source     = "TPI",
    Model      = tolower(Model),
    Resolution = as.numeric(Resolution),
    AIC        = as.numeric(AIC)
  ) %>%
  filter(!is.na(Resolution), !is.na(AIC)) %>%
  dplyr::select(Source, Model, Resolution, AIC)

aic_all <- bind_rows(oi_df, tpi_df) %>%
  arrange(Source, Model, Resolution)

#-----------------------------#
# 2) Compute delta-AIC relative to AIC at 0.05m resolution
#-----------------------------#

aic_delta <- aic_all %>%
  group_by(Source, Model) %>%
  mutate(
    AIC_baseline = AIC[Resolution == 0.05],
    delta_AIC    = AIC - AIC_baseline
  ) %>%
  ungroup()

#-----------------------------#
# 3) Console summaries
#-----------------------------#

# Best AIC resolution per model
best_aic <- aic_all %>%
  group_by(Source, Model) %>%
  slice_min(AIC, n = 1) %>%
  ungroup() %>%
  dplyr::select(Source, Model, Resolution, AIC)

cat("\n--- Resolution with best (lowest) AIC per model ---\n")
print(best_aic)

# First resolution where delta-AIC exceeds 2 (excluding the 0.05m baseline itself)
threshold_crossings <- aic_delta %>%
  filter(Resolution > 0.05, delta_AIC > 2) %>%
  arrange(Source, Model, Resolution) %>%
  group_by(Source, Model) %>%
  slice_min(Resolution, n = 1) %>%
  ungroup() %>%
  dplyr::select(Source, Model, Resolution, delta_AIC)

cat("\n--- First resolution where delta-AIC exceeds 2 (relative to 0.05m) ---\n")
if (nrow(threshold_crossings) == 0) {
  cat("No model crosses the delta-AIC = 2 threshold across tested resolutions.\n")
} else {
  print(threshold_crossings)
}

#-----------------------------#
# 4) Plot
#-----------------------------#

ggplot(aic_delta, aes(x = Resolution, y = delta_AIC, color = Model)) +
  geom_line(size = 1, alpha = 0.8) +
  geom_point(size = 2) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "black") +
  facet_wrap(~ Source, scales = "free_x") +
  labs(
    x     = "Resolution (m)",
    y     = "ΔAIC (relative to 0.05m resolution)",
    color = "Model",
    title = "ΔAIC across resolutions (dashed = ΔAIC threshold of 2)"
  ) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

#---------------------------------------------------------------#
# Plot p-value vs OI resolution with 0.05 threshold line        #
#---------------------------------------------------------------#

library(dplyr)
library(ggplot2)

# Path to your model performance CSV
perf_path <- "E:/PhD/Chap 4 - Resolution/OI_resolution_model_performance.csv"

# Read data
oi_perf <- read.csv(perf_path, stringsAsFactors = FALSE)

# Make sure Resolution and OI_p_value are numeric and remove NAs
oi_perf <- oi_perf %>%
  mutate(
    Resolution = as.numeric(Resolution),
    OI_p_value = as.numeric(OI_p_value)
  ) %>%
  filter(!is.na(Resolution), !is.na(OI_p_value))

# Basic check
# head(oi_perf)

#---------------------------------------------------------------#
# Plot: p-value for the OI term vs resolution                   #
# - X: Resolution                                               #
# - Y: OI_p_value                                               #
# - Points only (no connecting lines)                           #
# - 0.05 threshold line on y-axis                               #
# - Color by model type (abundance vs occurrence)               #
#---------------------------------------------------------------#

ggplot(oi_perf, aes(x = Resolution, y = OI_p_value, color = Model)) +
  geom_line(size = 1, alpha = 0.8) +          # points only, no lines
  geom_hline(yintercept = 0.05,                  # p = 0.05 threshold
             linetype = "dashed",
             color = "black") +
  labs(
    x = "OI Resolution (m)",
    y = "p-value for OI term",
    color = "Model",
    title = "Change in p-value for OI across resolutions"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )



#------------------------------------------------------------------#
#  Table-plots for abundance & occurrence models, 15–26 m          #
#  Drops Model and OI_column, adds titles to each table plot       #
#------------------------------------------------------------------#

# install.packages("gridExtra")  # if not already installed

library(gridExtra)
library(grid)

# Path to your performance CSV
perf_path <- "E:/PhD/Chap 4 - Resolution/OI_resolution_model_performance.csv"

# Read data
perf <- read.csv(perf_path, stringsAsFactors = FALSE)

# Ensure Resolution is numeric
perf$Resolution <- as.numeric(perf$Resolution)

# Columns to keep: drop Model and OI_column
cols_keep <- setdiff(colnames(perf), c("Model", "OI_column"))

#------------------------------------------------------------------#
# 1) Abundance model table (15–26 m)                               #
#------------------------------------------------------------------#

# Logical index for abundance rows in 15–26 m
idx_abund <- perf$Model == "abundance" &
  !is.na(perf$Resolution) &
  perf$Resolution >= 20 &
  perf$Resolution <= 30

abund_tab <- perf[idx_abund, cols_keep, drop = FALSE]
abund_tab <- abund_tab[order(abund_tab$Resolution), ]

# Table grob
abund_table_grob <- tableGrob(abund_tab, rows = NULL)

# Title grob
abund_title <- textGrob(
  "OI performance (Abundance model, 20–30 m)",
  gp = gpar(fontsize = 14, fontface = "bold")
)

# Combine title + table into one grob
abund_plot <- arrangeGrob(
  abund_title,
  abund_table_grob,
  ncol = 1,
  heights = c(0.1, 0.9)
)

# Draw abundance table-plot (you can copy/export as image)
grid.newpage()
grid.draw(abund_plot)

#------------------------------------------------------------------#
# 2) Occurrence model table (15–26 m)                              #
#------------------------------------------------------------------#

idx_occ <- perf$Model == "occurrence" &
  !is.na(perf$Resolution) &
  perf$Resolution >= 10 &
  perf$Resolution <= 26

occ_tab <- perf[idx_occ, cols_keep, drop = FALSE]
occ_tab <- occ_tab[order(occ_tab$Resolution), ]

occ_table_grob <- tableGrob(occ_tab, rows = NULL)

occ_title <- textGrob(
  "OI performance (Occurrence model, 10–26 m)",
  gp = gpar(fontsize = 14, fontface = "bold")
)

occ_plot <- arrangeGrob(
  occ_title,
  occ_table_grob,
  ncol = 1,
  heights = c(0.1, 0.9)
)

grid.newpage()
grid.draw(occ_plot)

#------------------------------------------------------------------------------#
#  Breakpoint & "knee" analysis for ALL MODELS (OI + TPI)                      #
#  Outputs detailed summaries to CSV files                                      #
#------------------------------------------------------------------------------#

library(dplyr)
library(readr)
library(purrr)
library(changepoint)
library(strucchange)
library(segmented)

#--------------------------- Paths & output dir --------------------------------#
base_dir      <- "E:/PhD/Chap 4 - Resolution"
oi_perf_path  <- file.path(base_dir, "OI_resolution_model_performance.csv")
tpi_perf_path <- file.path(base_dir, "TPI_resolution_model_performance.csv")

out_dir <- file.path(base_dir, "Breakpoint_Analysis")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#--------------------------- Helper: knee finder -------------------------------#
find_knee <- function(res, p) {
  ord <- order(res)
  x   <- res[ord]
  y   <- -log10(p[ord])
  
  # If all y equal, no knee
  if (length(unique(y)) <= 1 || max(x) == min(x)) {
    return(list(
      idx = NA_integer_,
      knee_res = NA_real_,
      knee_p = NA_real_,
      knee_dist = NA_real_
    ))
  }
  
  # Normalize to [0,1]
  x_norm <- (x - min(x)) / (max(x) - min(x))
  y_norm <- (y - min(y)) / (max(y) - min(y))
  
  p1 <- c(x_norm[1], y_norm[1])
  p2 <- c(x_norm[length(x_norm)], y_norm[length(x_norm)])
  v  <- p2 - p1
  
  num <- abs(v[2] * x_norm - v[1] * y_norm + p2[1] * p1[2] - p2[2] * p1[1])
  den <- sqrt(v[1]^2 + v[2]^2)
  d   <- num / den
  
  i <- which.max(d)
  list(
    idx       = ord[i],
    knee_res  = res[ord[i]],
    knee_p    = p[ord[i]],
    knee_dist = d[i]
  )
}

#--------------------------- 1. KNEE ANALYSIS ---------------------------------#
# OI models (abundance, occurrence)
oi_perf <- read_csv(oi_perf_path, show_col_types = FALSE) %>%
  dplyr::filter(!is.na(Resolution), !is.na(OI_p_value)) %>%
  dplyr::mutate(
    Resolution = as.numeric(Resolution),
    OI_p_value = as.numeric(OI_p_value)
  ) %>%
  dplyr::arrange(Model, Resolution)

knee_OI <- oi_perf %>%
  dplyr::group_by(Model) %>%
  dplyr::group_modify(~{
    d  <- .x %>% dplyr::arrange(Resolution)
    r0 <- d$Resolution[1]
    p0 <- d$OI_p_value[1]
    
    d <- d %>%
      dplyr::mutate(
        delta_p    = OI_p_value - p0,
        logp       = -log10(OI_p_value),
        logp_base  = -log10(p0),
        delta_logp = logp - logp_base
      )
    
    knee <- find_knee(d$Resolution, d$OI_p_value)
    idx_max_delta <- which.max(d$delta_p)
    res_max_delta <- d$Resolution[idx_max_delta]
    
    tibble::tibble(
      Source              = "OI",
      Model               = unique(.x$Model),
      baseline_res        = r0,
      baseline_p          = p0,
      max_delta_res       = res_max_delta,
      max_delta_p         = d$delta_p[idx_max_delta],
      max_delta_logp      = d$delta_logp[idx_max_delta],
      knee_resolution     = knee$knee_res,
      knee_p              = knee$knee_p,
      knee_distance_score = knee$knee_dist
    )
  }) %>%
  dplyr::ungroup()

# TPI models (richness, diversity)
tpi_perf <- read_csv(tpi_perf_path, show_col_types = FALSE) %>%
  dplyr::filter(!is.na(Resolution), !is.na(TPI_p_value)) %>%
  dplyr::mutate(
    Resolution  = as.numeric(Resolution),
    TPI_p_value = as.numeric(TPI_p_value)
  ) %>%
  dplyr::arrange(Model, Resolution)

knee_TPI <- tpi_perf %>%
  dplyr::group_by(Model) %>%
  dplyr::group_modify(~{
    d  <- .x %>% dplyr::arrange(Resolution)
    r0 <- d$Resolution[1]
    p0 <- d$TPI_p_value[1]
    
    d <- d %>%
      dplyr::mutate(
        delta_p    = TPI_p_value - p0,
        logp       = -log10(TPI_p_value),
        logp_base  = -log10(p0),
        delta_logp = logp - logp_base
      )
    
    knee <- find_knee(d$Resolution, d$TPI_p_value)
    idx_max_delta <- which.max(d$delta_p)
    res_max_delta <- d$Resolution[idx_max_delta]
    
    tibble::tibble(
      Source              = "TPI",
      Model               = unique(.x$Model),
      baseline_res        = r0,
      baseline_p          = p0,
      max_delta_res       = res_max_delta,
      max_delta_p         = d$delta_p[idx_max_delta],
      max_delta_logp      = d$delta_logp[idx_max_delta],
      knee_resolution     = knee$knee_res,
      knee_p              = knee$knee_p,
      knee_distance_score = knee$knee_dist
    )
  }) %>%
  dplyr::ungroup()

knee_all <- dplyr::bind_rows(knee_OI, knee_TPI)
readr::write_csv(knee_all, file.path(out_dir, "knee_summary_all_models.csv"))


#--------------------------- 2. CHANGEPOINT ANALYSIS --------------------------#
run_changepoint <- function(perf_path, p_col, model_name, source_label) {
  df <- read_csv(perf_path, show_col_types = FALSE) %>%
    dplyr::filter(Model == model_name) %>%
    dplyr::mutate(
      Resolution = as.numeric(Resolution),
      p_val      = as.numeric(.data[[p_col]])
    ) %>%
    dplyr::filter(!is.na(Resolution), !is.na(p_val)) %>%
    dplyr::arrange(Resolution)
  
  if (nrow(df) == 0) {
    return(tibble::tibble(
      Source            = source_label,
      Model             = model_name,
      n                 = 0L,
      cp_count          = 0L,
      cp_indices        = NA_character_,
      cp_resolutions    = NA_character_,
      note              = "No data rows"
    ))
  }
  
  logp_vec <- -log10(df$p_val + 1e-12)
  cp_obj   <- cpt.mean(logp_vec, method = "PELT", penalty = "AIC", Q = 1)
  idx_vec  <- cpts(cp_obj)
  
  if (length(idx_vec) == 0) {
    tibble::tibble(
      Source         = source_label,
      Model          = model_name,
      n              = nrow(df),
      cp_count       = 0L,
      cp_indices     = NA_character_,
      cp_resolutions = NA_character_,
      note           = "No changepoints found"
    )
  } else {
    tibble::tibble(
      Source         = source_label,
      Model          = model_name,
      n              = nrow(df),
      cp_count       = length(idx_vec),
      cp_indices     = paste(idx_vec, collapse = ";"),
      cp_resolutions = paste(df$Resolution[idx_vec], collapse = ";"),
      note           = "Changepoints detected"
    )
  }
}

cp_OI <- dplyr::bind_rows(
  run_changepoint(oi_perf_path,  "OI_p_value",  "abundance",  "OI"),
  run_changepoint(oi_perf_path,  "OI_p_value",  "occurrence", "OI")
)

cp_TPI <- dplyr::bind_rows(
  run_changepoint(tpi_perf_path, "TPI_p_value", "richness",   "TPI"),
  run_changepoint(tpi_perf_path, "TPI_p_value", "diversity",  "TPI")
)

cp_all <- dplyr::bind_rows(cp_OI, cp_TPI)
readr::write_csv(cp_all, file.path(out_dir, "changepoint_summary_all_models.csv"))


#--------------------------- 3. STRUCCHANGE ANALYSIS --------------------------#
run_strucchange <- function(perf_path, p_col, model_name, source_label) {
  df <- read_csv(perf_path, show_col_types = FALSE) %>%
    dplyr::filter(Model == model_name) %>%
    dplyr::mutate(
      Resolution = as.numeric(Resolution),
      p_val      = as.numeric(.data[[p_col]])
    ) %>%
    dplyr::filter(!is.na(Resolution), !is.na(p_val)) %>%
    dplyr::arrange(Resolution) %>%
    dplyr::mutate(logp = -log10(p_val + 1e-12))
  
  if (nrow(df) == 0) {
    return(tibble::tibble(
      Source          = source_label,
      Model           = model_name,
      n               = 0L,
      lm_intercept    = NA_real_,
      lm_slope        = NA_real_,
      lm_slope_p      = NA_real_,
      lm_R2           = NA_real_,
      lm_adjR2        = NA_real_,
      lm_sigma        = NA_real_,
      BIC_values      = NA_character_,
      m_opt           = NA_integer_,
      break_indices   = NA_character_,
      break_res       = NA_character_
    ))
  }
  
  fm <- lm(logp ~ Resolution, data = df)
  sm <- summary(fm)
  
  slope_row <- sm$coefficients["Resolution", , drop = FALSE]
  slope_est <- slope_row[1, "Estimate"]
  slope_p   <- slope_row[1, "Pr(>|t|)"]
  
  bp_full  <- breakpoints(logp ~ Resolution, data = df)
  bic_vals <- BIC(bp_full)
  m_opt    <- which.min(bic_vals) - 1
  
  bic_str <- paste(round(bic_vals, 3), collapse = ";")
  
  if (m_opt == 0) {
    tibble::tibble(
      Source        = source_label,
      Model         = model_name,
      n             = nrow(df),
      lm_intercept  = sm$coefficients["(Intercept)", "Estimate"],
      lm_slope      = slope_est,
      lm_slope_p    = slope_p,
      lm_R2         = sm$r.squared,
      lm_adjR2      = sm$adj.r.squared,
      lm_sigma      = sm$sigma,
      BIC_values    = bic_str,
      m_opt         = m_opt,
      break_indices = NA_character_,
      break_res     = NA_character_
    )
  } else {
    bp_opt <- breakpoints(bp_full, breaks = m_opt)
    idx    <- bp_opt$breakpoints
    idx    <- idx[!is.na(idx)]
    res    <- df$Resolution[idx]
    
    tibble::tibble(
      Source        = source_label,
      Model         = model_name,
      n             = nrow(df),
      lm_intercept  = sm$coefficients["(Intercept)", "Estimate"],
      lm_slope      = slope_est,
      lm_slope_p    = slope_p,
      lm_R2         = sm$r.squared,
      lm_adjR2      = sm$adj.r.squared,
      lm_sigma      = sm$sigma,
      BIC_values    = bic_str,
      m_opt         = m_opt,
      break_indices = paste(idx, collapse = ";"),
      break_res     = paste(res, collapse = ";")
    )
  }
}

sc_OI <- dplyr::bind_rows(
  run_strucchange(oi_perf_path,  "OI_p_value",  "abundance",  "OI"),
  run_strucchange(oi_perf_path,  "OI_p_value",  "occurrence", "OI")
)

sc_TPI <- dplyr::bind_rows(
  run_strucchange(tpi_perf_path, "TPI_p_value", "richness",   "TPI"),
  run_strucchange(tpi_perf_path, "TPI_p_value", "diversity",  "TPI")
)

sc_all <- dplyr::bind_rows(sc_OI, sc_TPI)
readr::write_csv(sc_all, file.path(out_dir, "strucchange_summary_all_models.csv"))


#--------------------------- 4. SEGMENTED ANALYSIS ----------------------------#
run_segmented <- function(perf_path, p_col, model_name, source_label) {
  df <- read_csv(perf_path, show_col_types = FALSE) %>%
    dplyr::filter(Model == model_name) %>%
    dplyr::mutate(
      Resolution = as.numeric(Resolution),
      p_val      = as.numeric(.data[[p_col]])
    ) %>%
    dplyr::filter(!is.na(Resolution), !is.na(p_val)) %>%
    dplyr::arrange(Resolution) %>%
    dplyr::mutate(logp = -log10(p_val + 1e-12))
  
  if (nrow(df) == 0) {
    return(tibble::tibble(
      Source          = source_label,
      Model           = model_name,
      n               = 0L,
      init_psi        = NA_real_,
      n_breaks        = NA_integer_,
      break_res       = NA_character_,
      seg_R2          = NA_real_,
      seg_adjR2       = NA_real_,
      status          = "No data rows"
    ))
  }
  
  fm <- lm(logp ~ Resolution, data = df)
  init_psi <- stats::median(df$Resolution, na.rm = TRUE)
  
  seg_fit <- try(segmented(fm, seg.Z = ~Resolution, psi = init_psi), silent = TRUE)
  
  if (inherits(seg_fit, "try-error")) {
    return(tibble::tibble(
      Source          = source_label,
      Model           = model_name,
      n               = nrow(df),
      init_psi        = init_psi,
      n_breaks        = NA_integer_,
      break_res       = NA_character_,
      seg_R2          = NA_real_,
      seg_adjR2       = NA_real_,
      status          = paste("segmented error:", as.character(seg_fit))
    ))
  }
  
  sm <- summary(seg_fit)
  br <- seg_fit$psi[, "Est."]
  tibble::tibble(
    Source    = source_label,
    Model     = model_name,
    n         = nrow(df),
    init_psi  = init_psi,
    n_breaks  = length(br),
    break_res = paste(round(br, 3), collapse = ";"),
    seg_R2    = sm$r.squared,
    seg_adjR2 = sm$adj.r.squared,
    status    = "OK"
  )
}

seg_OI <- dplyr::bind_rows(
  run_segmented(oi_perf_path,  "OI_p_value",  "abundance",  "OI"),
  run_segmented(oi_perf_path,  "OI_p_value",  "occurrence", "OI")
)

seg_TPI <- dplyr::bind_rows(
  run_segmented(tpi_perf_path, "TPI_p_value", "richness",   "TPI"),
  run_segmented(tpi_perf_path, "TPI_p_value", "diversity",  "TPI")
)

seg_all <- dplyr::bind_rows(seg_OI, seg_TPI)
readr::write_csv(seg_all, file.path(out_dir, "segmented_summary_all_models.csv"))

cat("\nAll breakpoint and knee summaries written to:\n", out_dir, "\n")

#---------------------------------------------------------------#
# Plot -log10(p-values) vs resolution for all 4 models          #
#   (abundance, occurrence, richness, diversity)                #
#---------------------------------------------------------------#

library(dplyr)
library(readr)
library(ggplot2)

# Paths to performance CSVs
oi_perf_path  <- "E:/PhD/Chap 4 - Resolution/OI_resolution_model_performance.csv"
tpi_perf_path <- "E:/PhD/Chap 4 - Resolution/TPI_resolution_model_performance.csv"

#-----------------------------
# 1) Load OI-based models
#-----------------------------
oi_perf <- read_csv(oi_perf_path, show_col_types = FALSE) %>%
  dplyr::mutate(
    Model      = tolower(Model),
    Resolution = as.numeric(Resolution),
    p_value    = as.numeric(OI_p_value),
    Source     = "OI"
  ) %>%
  dplyr::filter(
    Model %in% c("abundance", "occurrence"),
    !is.na(Resolution),
    !is.na(p_value)
  ) %>%
  dplyr::select(Source, Model, Resolution, p_value)

#-----------------------------
# 2) Load TPI-based models
#-----------------------------
tpi_perf <- read_csv(tpi_perf_path, show_col_types = FALSE) %>%
  dplyr::mutate(
    Model      = tolower(Model),
    Resolution = as.numeric(Resolution),
    p_value    = as.numeric(TPI_p_value),
    Source     = "TPI"
  ) %>%
  dplyr::filter(
    Model %in% c("richness", "diversity"),
    !is.na(Resolution),
    !is.na(p_value)
  ) %>%
  dplyr::select(Source, Model, Resolution, p_value)

#-----------------------------
# 3) Combine and compute -log10(p)
#-----------------------------
all_perf <- bind_rows(oi_perf, tpi_perf) %>%
  dplyr::mutate(
    logp = -log10(p_value + 1e-12),  # avoid log(0)
    Model_plot = dplyr::case_when(
      Source == "OI"  & Model == "abundance"  ~ "Abundance (OI)",
      Source == "OI"  & Model == "occurrence" ~ "Occurrence (OI)",
      Source == "TPI" & Model == "richness"   ~ "Richness (TPI)",
      Source == "TPI" & Model == "diversity"  ~ "Diversity (TPI)",
      TRUE ~ paste(Source, Model)
    ),
    Model_plot = factor(
      Model_plot,
      levels = c("Abundance (OI)", "Occurrence (OI)",
                 "Richness (TPI)", "Diversity (TPI)")
    )
  )

# Optional quick check
# head(all_perf)

#-----------------------------
# 4) Plot: -log10(p) vs Resolution
#-----------------------------
ggplot(all_perf, aes(x = Resolution, y = logp, color = Model_plot)) +
  geom_line(size = 1, alpha = 0.8) +
  geom_point(size = 1.5, alpha = 0.7) +
  geom_hline(
    yintercept = -log10(0.05),   # significance threshold
    linetype   = "dashed"
  ) +
  labs(
    x     = "Resolution (m)",
    y     = "-log10(p-value) for thermal term",
    color = "Model",
    title = "Change in thermal-term significance across resolutions\n(OI and TPI models)"
  ) +
  theme_classic() +
  theme(
    plot.title      = element_text(hjust = 0.5),
    legend.position = "top"
  )

# ============================================================
# PLOT: mean_div vs resolution + breakpoints (all three)
# ============================================================
summary_res <- read.csv("F:/Chap 4 - Resolution/ThermalDiversity_Resolution_Analysis/ThermalDiversity_vs_Resolution_summary.csv")

plot(summary_res$resolution, summary_res$mean_div,
     pch = 19, cex = 0.6,
     xlab = "Resolution (m)",
     ylab = "Mean Simpson Diversity",
     main = "Thermal Diversity vs Resolution with Breakpoints")

# ----- Add strucchange model with 3 breakpoints -----
bp_mod3 <- breakpoints(mean_div ~ resolution, data = summary_res, breaks = 3)
fit_sc3 <- fitted(bp_mod3, breaks = 3)

lines(summary_res$resolution, fit_sc3, col = "red", lwd = 2)

bp_idx3 <- bp_mod3$breakpoints
bp_res3 <- summary_res$resolution[bp_idx3]

abline(v = bp_res3, col = "blue", lty = 2, lwd = 2)

# ----- Segmented regression breakpoint -----
seg_fit <- segmented(lm(mean_div ~ resolution, data = summary_res),
                     seg.Z = ~ resolution)
seg_bp <- seg_fit$psi[2]

abline(v = seg_bp, col = "darkgreen", lty = 3, lwd = 2)

# ----- Knee point -----
curv <- diff(summary_res$mean_div, differences = 2)
idx_all <- 2:(nrow(summary_res) - 1)
knee_pos <- which.max(abs(curv))
knee_res <- summary_res$resolution[idx_all[knee_pos]]

abline(v = knee_res, col = "purple", lty = 4, lwd = 2)

# ------------------------------------------------------------
# Legend FIX: bottom aligns at y = 0.7
# ------------------------------------------------------------

legend_height <- 0.12   # approximate height of legend box
legend_bottom <- 0.7
legend_top <- legend_bottom + legend_height

legend(x = min(summary_res$resolution),
       y = legend_top,   # TOP of legend goes here
       legend = c(
         "Mean diversity",
         #"Strucchange fit (3 segments)",
         #"Strucchange breakpoints (3)",
         "Segmented breakpoint",
         "Knee point"
       ),
       pch = c(19, NA, NA, NA, NA),
       lty = c(NA, 1, 2, 3, 4),
       col = c("black", "darkgreen", "purple"),
       #col = c("black", "red", "blue", "darkgreen", "purple"),
       bty = "o",
       bg = "white",
       cex = 0.9
)

# -----------------------------------------------------------
# Inset scatter: Pixel Count vs Resolution, highlighting
# the high-resolution range where pixel counts collapse
# -----------------------------------------------------------

library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)

# 1) Load data
pix_path <- "F:/Chap 4 - Resolution/AirTemp_pixel_counts_1_100_by_Code.csv"

pix_df <- read_csv(pix_path, show_col_types = FALSE)

# 2) Decide what "collapse" means: threshold on median pixel count
#    (e.g., when median count per resolution <= 10)
collapse_threshold <- 10

collapse_info <- pix_df %>%
  group_by(resolution) %>%
  summarise(
    median_pixels = median(`Pixel Count`, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(median_pixels))

# Resolutions where median_pixels is at or below the threshold
collapse_res <- collapse_info %>%
  filter(median_pixels <= collapse_threshold)

if (nrow(collapse_res) > 0) {
  inset_x_min <- min(collapse_res$resolution)
} else {
  # Fallback: just zoom in on the last 20 m of the range
  inset_x_min <- max(pix_df$resolution, na.rm = TRUE) - 20
}

inset_x_max <- max(pix_df$resolution, na.rm = TRUE)

# Data for inset: only the "collapse" x-range
pix_inset <- pix_df %>%
  filter(resolution >= inset_x_min)

# 3) Main plot: full resolution range, with rectangle showing inset region
p_main <- ggplot(pix_df,
                 aes(x = resolution, y = `Pixel Count`)) +
  geom_point(alpha = 0.5, size = 1.5) +
  # Rectangle to show the inset x-range
  annotate(
    "rect",
    xmin  = inset_x_min,
    xmax  = inset_x_max,
    ymin  = -Inf,
    ymax  = Inf,
    fill  = NA,
    colour = "black",
    linetype = "dashed",
    linewidth = 0.4
  ) +
  labs(
    x = "Resolution (m)",
    y = "Number of pixels",
    title = "Pixel count across resolutions"
  ) +
  theme_classic()

# 4) Inset plot: zoomed in on the collapse range
p_inset <- ggplot(pix_inset,
                  aes(x = resolution, y = `Pixel Count`)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_x_continuous(limits = c(inset_x_min, inset_x_max)) +
  labs(
    x = "Resolution (m)",
    y = "Number of pixels"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.background = element_rect(colour = "black", fill = "white", linewidth = 0.4),
    axis.title = element_text(size = 10),
    axis.text  = element_text(size = 10),
    plot.margin = margin(2, 2, 2, 2)
  )

# 5) Combine: main + inset
combined_plot <- ggdraw(p_main) +
  # Adjust x, y, width, height (0–1) to position inset
  draw_plot(p_inset, x = 0.39, y = 0.4, width = 0.59, height = 0.58)

# This will show in the RStudio Plots pane
print(combined_plot)

# (Optional) Save to file
# ggsave("F:/Chap 4 - Resolution/PixelCount_with_inset.png",
#        combined_plot, width = 8, height = 5, dpi = 300)
