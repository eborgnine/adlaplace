# India 0.25° grid-cell polygons (country-specific; not in canexposure)
#
# Builds analysis cells from a CAMS-style 0.25° template over India, attaches
# agro-eco-region IDs, and caches a GeoPackage with geo_id = terra cell index.

india_cells_gpkg <- function(path) {
  file.path(path, "india_cells.gpkg")
}

india_worldpop_tif <- function(path, year = 2020) {
  file.path(path, sprintf("ind_pd_%d_1km_UNadj.tif", as.integer(year)))
}

india_eco_regions_rds <- function(path) {
  file.path(path, "india_eco_regions.rds")
}

india_eco_regions_shp <- function() {
  candidates <- c(
    path.expand("~/Downloads/Agro_Ecological_regions/Agro_Ecological_regions.shp"),
    "/store/patrick/copernicus/agg/indiaEcoRegions_source.shp"
  )
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) {
    stop(
      "Agro-ecological regions shapefile not found. ",
      "Place Agro_Ecological_regions.shp under ~/Downloads/Agro_Ecological_regions/"
    )
  }
  hit
}

india_geo_id_ok <- function(cells) {
  gid <- as.character(cells$geo_id)
  length(unique(gid[!is.na(gid)])) > 100L &&
    !any(grepl("NA", gid, fixed = TRUE), na.rm = TRUE)
}

#' Map aggregation polygon index `id` → `geo_id` for India cells
attach_geo_id <- function(df, cell_geo_id) {
  data.table::setDT(df)
  ids <- as.character(df$id)
  if (length(ids) && all(ids %in% cell_geo_id)) {
    df[, geo_id := ids]
  } else {
    idx <- as.integer(ids)
    if (anyNA(idx) || any(idx < 1L) || any(idx > length(cell_geo_id))) {
      stop("Cannot map aggregation id to india cell geo_id")
    }
    df[, geo_id := cell_geo_id[idx]]
  }
  df
}

remap_geo_from_id <- function(x, cells) {
  x <- data.table::as.data.table(x)
  if (!"id" %in% names(x)) {
    if (!"geo_id" %in% names(x)) {
      stop("exposure table needs an id or geo_id column")
    }
    x[, geo_id := as.character(geo_id)]
    return(x)
  }
  attach_geo_id(x, as.character(cells$geo_id))
}

get_india_eco_regions <- function(path) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  rds <- india_eco_regions_rds(path)
  if (file.exists(rds)) {
    return(terra::unwrap(readRDS(rds)))
  }
  eco <- terra::vect(india_eco_regions_shp())
  saveRDS(terra::wrap(eco), rds)
  message("Wrote ", rds)
  eco
}

#' CAMS-style 0.25° template raster over India (cell values = 1:ncell)
india_template_rast <- function() {
  ext <- terra::ext(67.875, 97.625, 7.875, 35.625)
  r <- terra::rast(ext = ext, res = 0.25, crs = "EPSG:4326")
  terra::values(r) <- seq_len(terra::ncell(r))
  names(r) <- "cell"
  r
}

get_india_cells <- function(path) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  gpkg <- india_cells_gpkg(path)
  if (file.exists(gpkg) && file.info(gpkg)$size > 1e5) {
    cells <- terra::vect(gpkg)
    if (india_geo_id_ok(cells)) {
      message("Loading India cells from ", gpkg)
      return(cells)
    }
    message("Cached cells have broken geo_id; rebuilding.")
  }

  message("Building India cells GeoPackage under ", path)
  india <- geodata::gadm("IND", level = 0, path = path)
  template <- india_template_rast()
  india_ll <- terra::project(india, terra::crs(template))
  masked <- terra::mask(template, india_ll)
  cells <- terra::as.polygons(masked, dissolve = FALSE, values = TRUE)
  cells$geo_id <- as.character(cells$cell)

  eco <- get_india_eco_regions(path)
  eco_ll <- terra::project(eco, terra::crs(cells))
  cent <- terra::centroids(cells)
  rel <- terra::relate(cent, eco_ll, "within")
  eco_idx <- apply(rel, 1, function(xx) c(which(xx), NA_integer_)[1L])
  cells$eco_region <- as.integer(eco_idx)
  na_eco <- which(is.na(cells$eco_region))
  if (length(na_eco)) {
    nearest <- terra::nearest(cent[na_eco, ], eco_ll)
    cells$eco_region[na_eco] <- terra::values(nearest)[, "to_id"]
  }

  terra::writeVector(cells, gpkg, overwrite = TRUE)
  message("Wrote ", gpkg, " (", nrow(cells), " cells)")
  cells
}

#' Build year-week3-DOW stratum id used in India CC models
add_year_week3_dow <- function(dt, date_col = "date") {
  data.table::setDT(dt)
  if (!inherits(dt[[date_col]], "Date")) {
    dt[, (date_col) := as.Date(get(date_col))]
  }
  week <- as.integer(format(dt[[date_col]], "%U"))
  week3 <- 3L * ceiling(week / 3)
  dt[, year_week3_dow := paste(
    format(get(date_col), "%Y"),
    formatC(week3, width = 2, format = "d", flag = "0"),
    format(get(date_col), "%a"),
    sep = "_"
  )]
  invisible(dt)
}

#' Expand daily deaths to full case-crossover calendar within `date_range`
expand_india_cc_calendar <- function(
  deaths_daily,
  cells,
  date_range = NULL
) {
  data.table::setDT(deaths_daily)
  if (is.null(date_range)) {
    date_range <- range(deaths_daily$date, na.rm = TRUE)
  }
  date_range <- as.Date(date_range)
  all_dates <- data.frame(
    date = seq(date_range[1], date_range[2], by = "day")
  )
  all_dates$year_week3_dow <- add_year_week3_dow(all_dates)$year_week3_dow

  cell_df <- data.table::data.table(
    geo_id = as.character(cells$geo_id),
    eco_region = as.integer(cells$eco_region)
  )
  active_cells <- unique(deaths_daily[, .(geo_id, eco_region)])
  cell_df <- merge(cell_df, active_cells, by = c("geo_id", "eco_region"))
  all_dates <- data.table::as.data.table(all_dates)
  add_year_week3_dow(all_dates)
  grid <- data.table::CJ(
    geo_id = cell_df$geo_id,
    date = all_dates$date,
    unique = TRUE
  )
  grid <- merge(grid, cell_df, by = "geo_id", all.x = TRUE)
  grid <- merge(grid, all_dates[, .(date, year_week3_dow)], by = "date", all.x = TRUE)

  mort_cols <- intersect(
    c("mort_all", "mort_cir", "mort_resp"),
    names(deaths_daily)
  )
  if (!length(mort_cols)) {
    stop("deaths_daily needs mort_* columns")
  }

  deaths_daily[, geo_id := as.character(geo_id)]
  agg <- deaths_daily[, lapply(.SD, sum, na.rm = TRUE),
    by = .(date, geo_id, eco_region),
    .SDcols = mort_cols
  ]
  out <- merge(
    grid, agg,
    by = c("date", "geo_id", "eco_region"),
    all.x = TRUE
  )
  for (col in mort_cols) {
    out[is.na(get(col)), (col) := 0L]
  }
  out[, cell := geo_id]
  add_year_week3_dow(out)
  data.table::setorderv(out, c("cell", "eco_region", "date"))
  out
}

#' Join exposures and build model columns for India hiwp fits
prepare_india_cc <- function(
  deaths_cc,
  exposures,
  date_range = NULL,
  pm_cap = 1000,
  temp_floor = -5
) {
  data.table::setDT(deaths_cc)
  data.table::setDT(exposures)
  exposures[, geo_id := as.character(geo_id)]
  if (!inherits(exposures$date, "Date")) {
    exposures[, date := as.Date(date)]
  }
  if (!is.null(date_range)) {
    date_range <- as.Date(date_range)
    deaths_cc <- deaths_cc[date >= date_range[1] & date <= date_range[2]]
    exposures <- exposures[date >= date_range[1] & date <= date_range[2]]
  }

  x <- merge(
    deaths_cc,
    exposures[, .(date, geo_id, pm, temperature)],
    by = c("date", "geo_id"),
    all.x = FALSE
  )
  x[, pm := fifelse(is.na(pm), NA_real_, pm)]
  x[, temperature := fifelse(is.na(temperature), NA_real_, temperature)]
  x <- x[!is.na(pm) & !is.na(temperature)]

  x[, sqrt_pm := sqrt(pmin(pm, pm_cap))]
  x[, temp_div10 := pmax(temperature, temp_floor) / 10]
  x[, dom := as.integer(format(date, "%d"))]
  if (requireNamespace("Hmisc", quietly = TRUE)) {
    x[dom == Hmisc::monthDays(date), dom := 31L]
  }
  x[format(date, "%m%d") == "0101", dom := 0L]
  x[, dom := factor(dom, levels = as.character(0:31))]
  x[, dom := stats::relevel(dom, ref = "1")]
  x[, eco_region := factor(as.character(eco_region))]

  strata <- c("cell", "eco_region", "year_week3_dow")
  mort_cols <- grep("^mort_", names(x), value = TRUE)
  keep <- x[, .(keep = sum(unlist(.SD)) > 0L), by = strata, .SDcols = mort_cols]
  x <- x[keep[keep == TRUE], on = strata, nomatch = NULL]
  x
}

pick_example_cells <- function(cells, lon_lat = list(
  delhi = c(77.2, 28.6),
  chennai = c(80.2, 13.1)
)) {
  cent <- terra::centroids(cells)
  xy <- terra::crds(cent)
  out <- list()
  for (nm in names(lon_lat)) {
    target <- lon_lat[[nm]]
    d <- (xy[, 1] - target[1])^2 + (xy[, 2] - target[2])^2
    out[[nm]] <- as.character(cells$geo_id[which.min(d)])
  }
  out
}
