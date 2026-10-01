# =============================================================================
# Landsat 5 (TM), 7 (ETM+) and 8 (OLI) cross-sensor harmonization of EVI, NDMI, NBR
# End-to-end workflow accompanying the methodological report
# Study area: Alta Albegna (Roccalbegna, Semproniano, Castell'Azzara, Sorano), Tuscany
#
# Dependencies: terra, rstac, ggplot2, hexbin
# Input (one file): <ROOT>/cartography/study_area.shp
# Optional input  : <ROOT>/cartography/land_use.shp (field ucs07, regional land-cover codes)
# Output          : <ROOT>/harmonization/{01_raw,02_grid,03_indices,04_stats,05_figures,06_invariance,07_L7vsL8}
#
# Stages are independent and restartable (existing files are skipped).
# Set RUN_* flags below to choose what to execute.
# =============================================================================

library(terra)
library(rstac)

# ------------------------------- CONFIGURATION -------------------------------
ROOT        <- "I:/R_Workspace/ALBEGNA_VALLEY"   # change for your machine
OUT         <- file.path(ROOT, "harmonization")
CRS_OUT     <- "EPSG:32632"
CC_MAX      <- 40            # scene-level cloud cover threshold (%), screening only
MAX_GAP     <- 16            # max days between the two scenes of a pair
BUFFER_FRAC <- 0.20          # bounding-box enlargement on each side
B_BOOT      <- 1000          # bootstrap replicates (resampling of pairs)
N_HEX       <- 200000        # pixels per panel for hexbin colour scale
SEED        <- 1

RUN_DOMAIN     <- TRUE
RUN_DISCOVERY  <- FALSE      # candidate pairs query (pairs below are the selection used)
RUN_DOWNLOAD   <- TRUE
RUN_GRID       <- TRUE
RUN_INDICES    <- TRUE
RUN_STATS      <- TRUE
RUN_FIGURES    <- TRUE
RUN_INVARIANCE <- TRUE
RUN_L7_vs_L8   <- TRUE

IX  <- c("EVI", "NDMI", "NBR")
GO  <- c("COMPRESS=DEFLATE", "TILED=YES")
GOF <- c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES")
dirs <- c("01_raw", "02_grid", "03_indices", "04_stats", "05_figures", "06_invariance", "07_L7vsL8")
for (d in dirs) dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)
D <- function(...) file.path(OUT, ...)
Sys.setenv(GDAL_HTTP_MAX_RETRY = "5", GDAL_HTTP_RETRY_DELAY = "2",
           CPL_VSIL_CURL_ALLOWED_EXTENSIONS = ".tif", GDAL_DISABLE_READDIR_ON_OPEN = "EMPTY_DIR")

# ------------------------- SELECTED PAIRS (20) -------------------------------
# ref = older-generation sensor (x), other = newer-generation sensor (y)
P <- function(chain, season, n, ref, other)
  data.frame(chain, season, pair = sprintf("%s_%s_%02d", chain, season, n), ref, other)
pairs <- rbind(
  P("TM_ETM","spring",1,"LT05_L2SP_191030_20000504_02_T1","LE07_L2SP_191030_20000426_02_T1"),
  P("TM_ETM","spring",2,"LT05_L2SP_192030_20020501_02_T2","LE07_L2SP_192030_20020423_02_T1"),
  P("TM_ETM","spring",3,"LT05_L2SP_191031_20000504_02_T1","LE07_L2SP_191031_20000426_02_T1"),
  P("TM_ETM","spring",4,"LT05_L2SP_191031_20010523_02_T1","LE07_L2SP_191031_20010531_02_T1"),
  P("TM_ETM","spring",5,"LT05_L2SP_191030_20020526_02_T1","LE07_L2SP_191030_20020518_02_T1"),
  P("TM_ETM","summer",1,"LT05_L2SP_191031_20000723_02_T1","LE07_L2SP_191031_20000731_02_T1"),
  P("TM_ETM","summer",2,"LT05_L2SP_191031_20000808_02_T1","LE07_L2SP_191031_20000816_02_T1"),
  P("TM_ETM","summer",3,"LT05_L2SP_191031_20010726_02_T1","LE07_L2SP_191031_20010803_02_T1"),
  P("TM_ETM","summer",4,"LT05_L2SP_191030_20020627_02_T1","LE07_L2SP_191030_20020619_02_T1"),
  P("TM_ETM","summer",5,"LT05_L2SP_192030_19990914_02_T1","LE07_L2SP_192030_19990906_02_T1"),
  P("ETM_OLI","spring",1,"LE07_L2SP_191031_20220520_02_T1","LC08_L2SP_191031_20220517_02_T1"),
  P("ETM_OLI","spring",2,"LE07_L2SP_192030_20180419_02_T1","LC08_L2SP_192030_20180427_02_T1"),
  P("ETM_OLI","spring",3,"LE07_L2SP_191031_20170409_02_T1","LC08_L2SP_191031_20170417_02_T1"),
  P("ETM_OLI","spring",4,"LE07_L2SP_191030_20200417_02_T1","LC08_L2SP_191030_20200409_02_T1"),
  P("ETM_OLI","spring",5,"LE07_L2SP_191030_20150420_02_T1","LC08_L2SP_191030_20150412_02_T1"),
  P("ETM_OLI","summer",1,"LE07_L2SP_191031_20220705_02_T1","LC08_L2SP_191031_20220704_02_T1"),
  P("ETM_OLI","summer",2,"LE07_L2SP_191031_20170730_02_T1","LC08_L2SP_191031_20170807_02_T1"),
  P("ETM_OLI","summer",3,"LE07_L2SP_191030_20200722_02_T1","LC08_L2SP_191030_20200730_02_T1"),
  P("ETM_OLI","summer",4,"LE07_L2SP_192030_20170822_02_T1","LC08_L2SP_192030_20170830_02_T1"),
  P("ETM_OLI","summer",5,"LE07_L2SP_191030_20170815_02_T1","LC08_L2SP_191030_20170823_02_T1"))
pairs$pathrow <- paste0(substr(pairs$ref, 11, 13), "/", substr(pairs$ref, 14, 16))
write.csv(pairs, D("01_raw/pairs_selected.csv"), row.names = FALSE)

# Pairs excluded from regression (reasons documented in the report)
EXCLUDE <- c("ETM_OLI_spring_01",   # SWIR1 band of the ETM+ scene unreadable at source
             "TM_ETM_spring_02")    # TM Tier 2 scene, empty opacity band, R2 < 0.07

# ------------------------------ HELPERS --------------------------------------
sensor_of <- function(id) substr(id, 1, 4)                     # LT05, LE07, LC08
aux_of    <- function(id) switch(sensor_of(id), LT05 = NULL, LE07 = "atmos_opacity", "qa_aerosol")
bands_of  <- function(id) c("blue", "green", "red", "nir08", "swir16", "swir22", "qa_pixel", aux_of(id))
sc        <- function(x) x * 0.0000275 - 0.2                  # Collection 2 L2 scaling

# Pixel-level QC. TM: QA_PIXEL only (opacity threshold does not transfer to TM, see report 2.3)
qc_mask <- function(r, id) {
  qa <- r$qa_pixel
  m <- switch(sensor_of(id),
    LT05 = qa == 5440,
    LE07 = (qa == 5440) & (r$atmos_opacity < 250),
    (qa == 21824) & (r$qa_aerosol %in% c(64, 66, 96)))
  rho <- sc(r[[c("blue", "red", "nir08", "swir16", "swir22")]])
  m & (sum((rho >= 0) & (rho <= 1)) == 5)
}
calc_idx <- function(r) {
  b <- sc(r$blue); rd <- sc(r$red); n <- sc(r$nir08); s1 <- sc(r$swir16); s2 <- sc(r$swir22)
  x <- c(2.5 * (n - rd) / (n + 6 * rd - 7.5 * b + 1), (n - s1) / (n + s1), (n - s2) / (n + s2))
  names(x) <- IX
  ifel(x >= -1 & x <= 1, x, NA)
}
download_scene <- function(stac, id, crop_v, raw_dir) {
  d <- file.path(raw_dir, id); dir.create(d, showWarnings = FALSE)
  f <- NULL
  for (k in 1:3) {
    f <- tryCatch((stac_search(stac, collections = "landsat-c2-l2", ids = id) |> get_request() |>
                   items_sign(sign_planetary_computer()))$features[[1]], error = function(e) NULL)
    if (!is.null(f)) break
  }
  if (is.null(f)) return(data.frame(scene = id, asset = NA, status = "NOT_FOUND"))
  out <- list()
  for (a in bands_of(id)) {
    dst <- file.path(d, paste0(id, "_", a, ".tif"))
    if (file.exists(dst)) { out[[a]] <- data.frame(scene = id, asset = a, status = "EXISTS"); next }
    if (is.null(f$assets[[a]])) { out[[a]] <- data.frame(scene = id, asset = a, status = "NO_ASSET"); next }
    st <- "ERROR"
    for (k in 1:3) {
      st <- tryCatch({
        r  <- rast(paste0("/vsicurl/", f$assets[[a]]$href))
        cv <- project(crop_v, crs(r))
        if (is.null(terra::intersect(ext(r), ext(cv)))) "NO_OVERLAP" else {
          writeRaster(crop(r, ext(cv), snap = "out"), dst, datatype = datatype(r), overwrite = TRUE, gdal = GO); "OK" }
      }, error = function(e) paste("ERROR:", substr(conditionMessage(e), 1, 60)))
      if (!grepl("^ERROR", st)) break
    }
    out[[a]] <- data.frame(scene = id, asset = a, status = st)
    if (st == "NO_OVERLAP") break
  }
  do.call(rbind, out)
}
# Scenes -> common grid (nearest neighbour, no interpolation of reflectance)
to_grid <- function(id, g) {
  fs <- file.path(D("01_raw"), id, paste0(id, "_", bands_of(id), ".tif"))
  if (!all(file.exists(fs))) return(NULL)
  r <- tryCatch({ x <- rast(fs); names(x) <- bands_of(id); project(x, g, method = "near") }, error = function(e) NULL)
  r
}
# Sufficient-statistics regression
fit_ss <- function(d) {
  n <- sum(d$n); Sx <- sum(d$Sx); Sy <- sum(d$Sy)
  sxx <- sum(d$Sxx) - Sx^2 / n; sxy <- sum(d$Sxy) - Sx * Sy / n; syy <- sum(d$Syy) - Sy^2 / n
  b <- sxy / sxx; a <- (Sy - b * Sx) / n
  c(n = n, slope = b, int = a, r2 = sxy^2 / (sxx * syy), rmse = sqrt(max(syy - b * sxy, 0) / (n - 2)))
}
x_file <- function(p) D("03_indices", paste0(p, "_x.tif"))
y_file <- function(p) D("03_indices", paste0(p, "_y.tif"))

# ============================ STAGE 1: DOMAIN ================================
# Bounding box of the study area enlarged by 20% on each side
crop_shp <- file.path(ROOT, "cartography/crop_area.shp")
if (RUN_DOMAIN || !file.exists(crop_shp)) {
  sa <- project(vect(file.path(ROOT, "cartography/study_area.shp")), CRS_OUT)
  e  <- ext(sa); w <- xmax(e) - xmin(e); h <- ymax(e) - ymin(e)
  ce <- ext(xmin(e) - BUFFER_FRAC * w, xmax(e) + BUFFER_FRAC * w,
            ymin(e) - BUFFER_FRAC * h, ymax(e) + BUFFER_FRAC * h)
  writeVector(as.polygons(ce, crs = CRS_OUT), crop_shp, overwrite = TRUE)
}
crop_v <- vect(crop_shp)

# ========================= STAGE 2: CANDIDATE PAIRS ==========================
# Candidate discovery (not needed to reproduce the analysis: selection is hard-coded above)
if (RUN_DISCOVERY) {
  bb   <- as.vector(ext(project(crop_v, "EPSG:4326")))[c(1, 3, 2, 4)]
  stac <- stac("https://planetarycomputer.microsoft.com/api/stac/v1")
  getc <- function(dt) {
    it <- stac_search(stac, collections = "landsat-c2-l2", bbox = bb, datetime = dt, limit = 500) |>
      get_request() |> items_fetch()
    x <- do.call(rbind, lapply(it$features, function(f) { p <- f$properties
      data.frame(id = f$id, plat = p$platform, date = as.Date(substr(p$datetime, 1, 10)),
                 path = p$`landsat:wrs_path`, row = p$`landsat:wrs_row`, cc = p$`eo:cloud_cover`) }))
    x[!is.na(x$cc) & x$cc <= CC_MAX, ]
  }
  season <- function(d) { m <- as.integer(format(d, "%m")); ifelse(m %in% 4:5, "spring", ifelse(m %in% 6:9, "summer", NA)) }
  mk <- function(x, p1, p2, chain) {
    g <- merge(x[x$plat == p1, ], x[x$plat == p2, ], by = c("path", "row"), suffixes = c("_1", "_2"))
    g$gap <- abs(as.integer(g$date_1 - g$date_2)); g$season <- season(g$date_1)
    g <- g[g$gap <= MAX_GAP & !is.na(g$season) & g$season == season(g$date_2), ]
    g$chain <- chain; g$cc_tot <- g$cc_1 + g$cc_2; g
  }
  cand <- rbind(mk(getc("1999-04-01T00:00:00Z/2003-05-30T23:59:59Z"), "landsat-5", "landsat-7", "TM_ETM"),
                mk(getc("2013-04-01T00:00:00Z/2022-12-31T23:59:59Z"), "landsat-7", "landsat-8", "ETM_OLI"))
  cand <- cand[order(cand$chain, cand$season, cand$cc_tot), ]
  write.csv(cand, D("01_raw/candidate_pairs.csv"), row.names = FALSE)
}

# =========================== STAGE 3: DOWNLOAD ===============================
if (RUN_DOWNLOAD) {
  stac <- stac("https://planetarycomputer.microsoft.com/api/stac/v1")
  ids  <- unique(c(pairs$ref, pairs$other)); lg <- list()
  for (id in ids) lg[[id]] <- download_scene(stac, id, crop_v, D("01_raw"))
  lg <- do.call(rbind, lg); write.csv(lg, D("01_raw/download_log.csv"), row.names = FALSE)
  print(table(lg$status))
}

# ======================= STAGE 4: COMMON 30 m GRID ===========================
# Cell edges aligned to the Landsat lattice (15 m offset)
e <- ext(crop_v)
GRID <- rast(ext(floor((xmin(e) - 15) / 30) * 30 + 15, ceiling((xmax(e) - 15) / 30) * 30 + 15,
                 floor((ymin(e) - 15) / 30) * 30 + 15, ceiling((ymax(e) - 15) / 30) * 30 + 15),
             resolution = 30, crs = CRS_OUT)

# ===================== STAGE 5: QC MASK AND INDICES PER PAIR =================
# A pixel is kept only if valid in both scenes and for all three indices.
if (RUN_INDICES) {
  log <- list()
  for (i in seq_len(nrow(pairs))) {
    p <- pairs$pair[i]
    if (file.exists(x_file(p))) next
    a <- to_grid(pairs$ref[i], GRID); b <- to_grid(pairs$other[i], GRID)
    if (is.null(a) || is.null(b)) { log[[p]] <- data.frame(pair = p, status = "SCENE_MISSING", n = NA); next }
    m  <- qc_mask(a, pairs$ref[i]) & qc_mask(b, pairs$other[i])
    xa <- calc_idx(a); xb <- calc_idx(b)
    v  <- ifel(m & !is.na(xa$EVI + xa$NDMI + xa$NBR + xb$EVI + xb$NDMI + xb$NBR), 1, NA)
    xa <- xa * v; xb <- xb * v
    writeRaster(xa, x_file(p), datatype = "FLT4S", NAflag = -9999, overwrite = TRUE, gdal = GOF)
    writeRaster(xb, y_file(p), datatype = "FLT4S", NAflag = -9999, overwrite = TRUE, gdal = GOF)
    log[[p]] <- data.frame(pair = p, status = "OK", n = global(!is.na(v), "sum")[1, 1])
  }
  if (length(log)) { lg <- do.call(rbind, log); write.csv(lg, D("03_indices/qc_summary.csv"), row.names = FALSE); print(lg) }
}

# ====================== STAGE 6: REGRESSION AND POOLING ======================
# OLS y ~ x per pair; pooled by summing sufficient statistics; bootstrap over pairs
use <- pairs[!pairs$pair %in% EXCLUDE & file.exists(x_file(pairs$pair)), ]
if (RUN_STATS) {
  set.seed(SEED); S <- list(); sub_xy <- list()
  for (i in seq_len(nrow(use))) {
    X <- rast(x_file(use$pair[i])); names(X) <- IX
    Y <- rast(y_file(use$pair[i])); names(Y) <- IX
    vx <- values(X); vy <- values(Y)
    for (k in IX) {
      x <- vx[, k]; y <- vy[, k]; ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; n <- length(x)
      if (n < 1000) next
      S[[paste(use$pair[i], k)]] <- data.frame(pair = use$pair[i], chain = use$chain[i], season = use$season[i],
        pr = use$pathrow[i], index = k, n = n, Sx = sum(x), Sy = sum(y), Sxx = sum(x^2), Sxy = sum(x * y), Syy = sum(y^2))
      j <- sample.int(n, min(n, 20000))
      sub_xy[[paste(use$pair[i], k)]] <- data.frame(chain = use$chain[i], season = use$season[i], index = k, x = x[j], y = y[j])
    }
  }
  S <- do.call(rbind, S); rownames(S) <- NULL
  pf <- do.call(rbind, lapply(seq_len(nrow(S)), function(i)
    cbind(S[i, c("pair", "chain", "season", "pr", "index")], t(fit_ss(S[i, ])))))
  write.csv(pf, D("04_stats/pair_fits.csv"), row.names = FALSE)

  cells <- unique(S[, c("chain", "season", "index")]); pooled <- list()
  for (c in seq_len(nrow(cells))) {
    sel <- S$chain == cells$chain[c] & S$season == cells$season[c] & S$index == cells$index[c]
    d <- S[sel, ]; f0 <- fit_ss(d); np <- nrow(d)
    bs <- t(replicate(B_BOOT, fit_ss(d[sample.int(np, np, replace = TRUE), ])[c("slope", "int")]))
    p1 <- pf[pf$chain == cells$chain[c] & pf$season == cells$season[c] & pf$index == cells$index[c], ]
    pooled[[c]] <- data.frame(cells[c, ], n_pairs = np, n_pix = f0[["n"]], slope = f0[["slope"]], int = f0[["int"]],
      r2 = f0[["r2"]], rmse = f0[["rmse"]],
      slope_lo = unname(quantile(bs[, 1], .025)), slope_hi = unname(quantile(bs[, 1], .975)),
      int_lo = unname(quantile(bs[, 2], .025)), int_hi = unname(quantile(bs[, 2], .975)),
      slope_sd_pairs = sd(p1$slope), slope_min = min(p1$slope), slope_max = max(p1$slope))
  }
  POOL <- do.call(rbind, pooled); rownames(POOL) <- NULL
  write.csv(POOL, D("04_stats/pooled_by_cell.csv"), row.names = FALSE)

  # EVI by path/row (Table 3)
  E <- S[S$index == "EVI", ]
  prt <- do.call(rbind, lapply(split(E, list(E$chain, E$season, E$pr), drop = TRUE), function(d)
    data.frame(chain = d$chain[1], season = d$season[1], pathrow = d$pr[1], n_pairs = nrow(d), t(round(fit_ss(d), 3)))))
  write.csv(prt, D("04_stats/evi_by_pathrow.csv"), row.names = FALSE)

  # Quadratic diagnostic (20,000-pixel subsample per pair)
  sx <- do.call(rbind, sub_xy)
  quad <- do.call(rbind, lapply(split(sx, list(sx$chain, sx$season, sx$index), drop = TRUE), function(d) {
    m1 <- lm(y ~ x, d); m2 <- lm(y ~ x + I(x^2), d)
    data.frame(chain = d$chain[1], season = d$season[1], index = d$index[1],
               r2_lin = summary(m1)$r.squared, r2_quad = summary(m2)$r.squared,
               d_r2 = summary(m2)$r.squared - summary(m1)$r.squared) }))
  write.csv(quad, D("04_stats/quad_check.csv"), row.names = FALSE)

  # Final summer functions. TM -> OLI by composition: OLI = (a2 + b2 a1) + (b2 b1) TM
  sm <- POOL[POOL$season == "summer", ]
  fin <- do.call(rbind, lapply(IX, function(k) {
    t <- sm[sm$chain == "TM_ETM" & sm$index == k, ]; e <- sm[sm$chain == "ETM_OLI" & sm$index == k, ]
    data.frame(index = k, TM_ETM_a = t$int, TM_ETM_b = t$slope, ETM_OLI_a = e$int, ETM_OLI_b = e$slope,
               TM_OLI_a = e$int + e$slope * t$int, TM_OLI_b = e$slope * t$slope) }))
  write.csv(fin, D("04_stats/final_summer_functions.csv"), row.names = FALSE)
  print(round(fin[, -1], 4)); print(POOL)
}

# =========================== STAGE 7: HEXBIN FIGURES =========================
# One figure per index, summer pairs only (top row TM-ETM+, bottom row ETM+-OLI)
if (RUN_FIGURES) {
  library(ggplot2)
  set.seed(SEED); sp <- use[use$season == "summer", ]; Dh <- list(); Sh <- list()
  for (i in seq_len(nrow(sp))) {
    X <- rast(x_file(sp$pair[i])); names(X) <- IX; Y <- rast(y_file(sp$pair[i])); names(Y) <- IX
    vx <- values(X); vy <- values(Y)
    for (k in IX) {
      x <- vx[, k]; y <- vy[, k]; ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
      m <- lm(y ~ x)
      Sh[[paste(sp$pair[i], k)]] <- data.frame(pair = sp$pair[i], index = k,
        txt = sprintf("n=%s\nslope=%.3f\nint=%.3f\nR2=%.3f", format(length(x), big.mark = ","),
                      coef(m)[2], coef(m)[1], summary(m)$r.squared))
      j <- sample.int(length(x), min(length(x), N_HEX))
      Dh[[paste(sp$pair[i], k)]] <- data.frame(pair = sp$pair[i], index = k, x = x[j], y = y[j])
    }
  }
  Dh <- do.call(rbind, Dh); Sh <- do.call(rbind, Sh)
  Dh$pair <- factor(Dh$pair, levels = sp$pair); Sh$pair <- factor(Sh$pair, levels = sp$pair)
  for (k in IX) {
    d <- Dh[Dh$index == k, ]; s <- Sh[Sh$index == k, ]
    lim <- quantile(c(d$x, d$y), c(0.001, 0.999))
    p <- ggplot(d, aes(x, y)) + geom_hex(bins = 60) +
      scale_fill_viridis_c(trans = "log10", name = "pixel count\n(log10)") +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "white", linewidth = 0.4) +
      geom_smooth(method = "lm", se = FALSE, colour = "red", linewidth = 0.5, formula = y ~ x) +
      geom_text(data = s, aes(x = -Inf, y = Inf, label = txt), inherit.aes = FALSE,
                hjust = -0.05, vjust = 1.1, size = 2, colour = "white", lineheight = 0.9) +
      facet_wrap(~pair, ncol = 5) + coord_equal(xlim = lim, ylim = lim) +
      labs(x = paste(k, "older sensor (TM for TM_ETM, ETM+ for ETM_OLI)"),
           y = paste(k, "newer sensor (ETM+ for TM_ETM, OLI for ETM_OLI)")) +
      theme_bw(base_size = 8) + theme(panel.grid = element_blank(), strip.text = element_text(size = 6))
    ggsave(D("05_figures", paste0("hexbin_summer_", k, ".png")), p, width = 14, height = 6.5, dpi = 300)
  }
}

# ========================= STAGE 8: INVARIANCE TEST ==========================
# Mean difference (y - x) before/after the pooled summer transformation on stable classes
lu_shp <- file.path(ROOT, "cartography/land_use.shp")
if (RUN_INVARIANCE && file.exists(lu_shp)) {
  sp <- use[use$season == "summer", ]; g0 <- rast(x_file(sp$pair[1]))[[1]]
  lu  <- crop(project(vect(lu_shp), crs(g0)), ext(g0))
  cod <- suppressWarnings(as.numeric(as.character(as.data.frame(lu)[["ucs07"]])))
  grp <- c("111" = "urban", "121" = "urban", "311" = "forest_broad", "312" = "forest_conif", "332" = "bare_rock")
  gid <- c(urban = 1, forest_broad = 2, forest_conif = 3, bare_rock = 4)
  lab <- grp[as.character(cod)]; keep <- !is.na(lab)
  v <- lu[keep]; v$g <- gid[lab[keep]]
  v <- buffer(v, width = -30); v <- v[!is.na(expanse(v)) & expanse(v) > 0]     # erosion: no mixed pixels
  cv <- values(rasterize(v, g0, field = "g"), mat = FALSE)
  sm <- read.csv(D("04_stats/pooled_by_cell.csv")); sm <- sm[sm$season == "summer", ]
  R <- list()
  for (i in seq_len(nrow(sp))) {
    X <- rast(x_file(sp$pair[i])); names(X) <- IX; Y <- rast(y_file(sp$pair[i])); names(Y) <- IX
    for (k in IX) {
      cf <- sm[sm$chain == sp$chain[i] & sm$index == k, ]
      x <- values(X[[k]], mat = FALSE); y <- values(Y[[k]], mat = FALSE)
      ok <- is.finite(x) & is.finite(y) & !is.na(cv)
      for (cl in names(gid)) {
        s <- ok & cv == gid[[cl]]; if (sum(s) < 50) next
        R[[paste(i, k, cl)]] <- data.frame(pair = sp$pair[i], chain = sp$chain[i], index = k, cls = cl, n = sum(s),
          bias_raw = mean(y[s] - x[s]), bias_adj = mean(y[s] - (cf$int + cf$slope * x[s])),
          sd_adj = sd(y[s] - (cf$int + cf$slope * x[s])))
      }
    }
  }
  R <- do.call(rbind, R); write.csv(R, D("06_invariance/invariance_by_pair.csv"), row.names = FALSE)
  A <- do.call(rbind, lapply(split(R, list(R$chain, R$index, R$cls), drop = TRUE), function(d) data.frame(
    chain = d$chain[1], index = d$index[1], cls = d$cls[1], pairs = nrow(d), n = sum(d$n),
    bias_raw = weighted.mean(d$bias_raw, d$n), bias_adj = weighted.mean(d$bias_adj, d$n),
    bias_adj_sd_pairs = if (nrow(d) > 1) sd(d$bias_adj) else NA, sd_adj = weighted.mean(d$sd_adj, d$n))))
  write.csv(A, D("06_invariance/invariance_summary.csv"), row.names = FALSE); print(A)
}

# ====================== STAGE 9: INDEPENDENT L7 vs L8 CHECK ==================
# Summer 2013-14 median composites of harmonized L7 (ETM+ -> OLI applied) and L8; difference L8 - L7
if (RUN_L7_vs_L8) {
  bb   <- as.vector(ext(project(crop_v, "EPSG:4326")))[c(1, 3, 2, 4)]
  stac <- stac("https://planetarycomputer.microsoft.com/api/stac/v1")
  inv <- list()
  for (y in 2013:2014) {
    it <- stac_search(stac, collections = "landsat-c2-l2", bbox = bb, limit = 500,
                      datetime = sprintf("%d-06-01T00:00:00Z/%d-09-30T23:59:59Z", y, y)) |> get_request() |> items_fetch()
    inv[[as.character(y)]] <- do.call(rbind, lapply(it$features, function(f) { p <- f$properties
      data.frame(id = f$id, plat = p$platform, cc = p$`eo:cloud_cover`) }))
  }
  inv <- do.call(rbind, inv)
  inv <- inv[!is.na(inv$cc) & inv$cc <= CC_MAX & grepl("_T1$", inv$id) & inv$plat %in% c("landsat-7", "landsat-8"), ]
  raw2 <- D("07_L7vsL8/raw"); dir.create(raw2, showWarnings = FALSE)
  for (id in inv$id) download_scene(stac, id, crop_v, raw2)

  fin <- read.csv(D("04_stats/final_summer_functions.csv"))
  idx_scene <- function(id) {
    fs <- file.path(raw2, id, paste0(id, "_", bands_of(id), ".tif"))
    if (!all(file.exists(fs))) return(NULL)
    r <- tryCatch({ x <- rast(fs); names(x) <- bands_of(id); project(x, GRID, method = "near") }, error = function(e) NULL)
    if (is.null(r)) return(NULL)
    x <- ifel(qc_mask(r, id), calc_idx(r), NA)
    if (sensor_of(id) == "LE07")
      for (k in seq_along(IX)) { f <- fin[fin$index == IX[k], ]; x[[k]] <- f$ETM_OLI_a + f$ETM_OLI_b * x[[k]] }
    x
  }
  L <- setNames(lapply(inv$id, idx_scene), inv$id); L <- L[!sapply(L, is.null)]
  l7 <- names(L)[sensor_of(names(L)) == "LE07"]; l8 <- names(L)[sensor_of(names(L)) == "LC08"]
  med <- function(ids, k) app(rast(lapply(ids, function(i) L[[i]][[k]])), "median", na.rm = TRUE)
  nob <- function(ids, k) app(rast(lapply(ids, function(i) L[[i]][[k]])), function(v) sum(!is.na(v)))
  res <- list()
  for (k in seq_along(IX)) {
    a <- med(l7, k); b <- med(l8, k); ok <- !is.na(b - a) & nob(l7, k) >= 3 & nob(l8, k) >= 3
    v <- values(b - a, mat = FALSE)[values(ok, mat = FALSE)]
    res[[k]] <- data.frame(index = IX[k], scenes_L7 = length(l7), scenes_L8 = length(l8), n_pix = length(v),
      mean_diff_L8_minus_L7 = mean(v), median = median(v), sd = sd(v),
      p05 = unname(quantile(v, .05)), p95 = unname(quantile(v, .95)))
  }
  res <- do.call(rbind, res); write.csv(res, D("07_L7vsL8/L7_vs_L8_2013_14.csv"), row.names = FALSE); print(res)
}

message("Done. Outputs in: ", OUT)