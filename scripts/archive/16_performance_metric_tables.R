#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(fs)
  library(openxlsx)
  library(stringr)
})

# =============================================================================
# Create performance metric tables from aggregated confusion matrices
#
# INPUT
#   data/out/09_performance/confusion_matrices.csv
#
# OUTPUT DIRECTORY
#   data/out/09_performance/performance_metric_tables/
#
# OUTPUT WORKBOOKS
#   precision.xlsx
#   recall.xlsx
#   specificity.xlsx
#   accuracy.xlsx
#   f1.xlsx
#   f2.xlsx
#   confusion_matrix_counts.xlsx
#
# Each metric workbook contains:
#   - heatindex
#   - precipitation
#
# Rows:
#   Native
#   STM1 bilinear
#   STM1 max
#   STM2 bilinear
#   STM2 max
#
# Columns:
#   2-year RP
#   5-year RP
#   10-year RP
#   25-year RP
# =============================================================================

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------

ROOT <- Sys.getenv(
  "CISC_ROOT",
  "/proj/sebastianlab/users/hquintal/cisc"
)

IN_FILE <- file.path(
  ROOT,
  "data/out/09_performance/confusion_matrices.csv"
)

OUT_DIR <- file.path(
  ROOT,
  "data/out/09_performance/performance_metric_tables"
)

dir_create(OUT_DIR, recurse = TRUE)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

RETURN_PERIODS <- c(2L, 5L, 10L, 25L)

METRICS <- c(
  "precision",
  "recall",
  "specificity",
  "accuracy",
  "f1",
  "f2"
)

COUNT_METRICS <- c(
  "TP",
  "FP",
  "FN",
  "TN",
  "total"
)

RESOLUTION_ORDER <- c(
  "Native",
  "STM1 bilinear",
  "STM1 max",
  "STM2 bilinear",
  "STM2 max"
)

OVERWRITE <- TRUE
METRIC_DIGITS <- 4L

# -----------------------------------------------------------------------------
# Logging
# -----------------------------------------------------------------------------

.log <- function(fmt, ...) {
  message(
    sprintf(
      "[%s] %s",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      sprintf(fmt, ...)
    )
  )
}

# -----------------------------------------------------------------------------
# Dataset labels
# -----------------------------------------------------------------------------

.dataset_label <- function(dataset_tag, variable) {
  
  dataset_tag <- as.character(dataset_tag)
  variable <- tolower(as.character(variable))
  
  out <- rep(NA_character_, length(dataset_tag))
  
  # Heat index
  out[
    variable == "heatindex" &
      dataset_tag == "heatindex_0.25"
  ] <- "Native"
  
  out[
    variable == "heatindex" &
      dataset_tag == "heatindex_0.30_bilinear"
  ] <- "STM1 bilinear"
  
  out[
    variable == "heatindex" &
      dataset_tag == "heatindex_0.30_max"
  ] <- "STM1 max"
  
  out[
    variable == "heatindex" &
      dataset_tag == "heatindex_0.36_bilinear"
  ] <- "STM2 bilinear"
  
  out[
    variable == "heatindex" &
      dataset_tag == "heatindex_0.36_max"
  ] <- "STM2 max"
  
  # Precipitation
  out[
    variable == "precipitation" &
      dataset_tag == "precipitation_0.25"
  ] <- "Native"
  
  out[
    variable == "precipitation" &
      dataset_tag == "precipitation_0.11_bilinear"
  ] <- "STM1 bilinear"
  
  out[
    variable == "precipitation" &
      dataset_tag == "precipitation_0.11_max"
  ] <- "STM1 max"
  
  out[
    variable == "precipitation" &
      dataset_tag == "precipitation_0.14_bilinear"
  ] <- "STM2 bilinear"
  
  out[
    variable == "precipitation" &
      dataset_tag == "precipitation_0.14_max"
  ] <- "STM2 max"
  
  out
}

# -----------------------------------------------------------------------------
# Read and validate
# -----------------------------------------------------------------------------

if (!file_exists(IN_FILE)) {
  stop(
    "Missing input file:\n",
    IN_FILE,
    "\n\nRun summarise_confusion_matrices_by_variable_rp.R first."
  )
}

dt <- fread(IN_FILE)

.log(
  "Input columns: %s",
  paste(names(dt), collapse = ", ")
)

print(dt[1:5])

required_columns <- c(
  "variable",
  "dataset_tag",
  "return_period",
  "TP",
  "FP",
  "FN",
  "TN"
)

missing_columns <- setdiff(required_columns, names(dt))

if (length(missing_columns)) {
  stop(
    "Input file is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

dt[, `:=`(
  variable = tolower(trimws(as.character(variable))),
  dataset_tag = as.character(dataset_tag),
  return_period = as.integer(return_period),
  TP = as.numeric(TP),
  FP = as.numeric(FP),
  FN = as.numeric(FN),
  TN = as.numeric(TN)
)]

dt[
  variable %in% c("heat index", "heat_index"),
  variable := "heatindex"
]

dt[
  variable %in% c("precip", "precipitation"),
  variable := "precipitation"
]

dt <- dt[
  variable %in% c("heatindex", "precipitation") &
    return_period %in% RETURN_PERIODS
]

# -----------------------------------------------------------------------------
# Recalculate all metrics from the summed confusion matrices
#
# This avoids depending on whether older copies of confusion_matrices.csv
# contain every desired derived metric.
# -----------------------------------------------------------------------------

dt[, total := TP + FP + FN + TN]

dt[, precision := fifelse(
  TP + FP > 0,
  TP / (TP + FP),
  NA_real_
)]

dt[, recall := fifelse(
  TP + FN > 0,
  TP / (TP + FN),
  NA_real_
)]

dt[, specificity := fifelse(
  TN + FP > 0,
  TN / (TN + FP),
  NA_real_
)]

dt[, accuracy := fifelse(
  total > 0,
  (TP + TN) / total,
  NA_real_
)]

dt[, f1 := fifelse(
  2 * TP + FP + FN > 0,
  2 * TP / (2 * TP + FP + FN),
  NA_real_
)]

dt[, f2 := fifelse(
  5 * TP + 4 * FN + FP > 0,
  5 * TP / (5 * TP + 4 * FN + FP),
  NA_real_
)]

dt[, resolution := .dataset_label(dataset_tag, variable)]

unknown_datasets <- unique(
  dt[is.na(resolution), .(variable, dataset_tag)]
)

if (nrow(unknown_datasets)) {
  warning(
    paste0(
      "The following dataset tags were not recognized and will be excluded:\n",
      paste(
        sprintf(
          "  %s: %s",
          unknown_datasets$variable,
          unknown_datasets$dataset_tag
        ),
        collapse = "\n"
      )
    )
  )
}

dt <- dt[!is.na(resolution)]

# Confirm uniqueness before reshaping.
duplicate_rows <- dt[
  ,
  .N,
  by = .(
    variable,
    resolution,
    return_period
  )
][N > 1L]

if (nrow(duplicate_rows)) {
  stop(
    paste0(
      "Multiple rows exist for at least one variable-resolution-RP combination.\n",
      "The input confusion_matrices.csv should contain one row per ",
      "variable, dataset_tag, and return_period.\n\n",
      paste(
        capture.output(print(duplicate_rows)),
        collapse = "\n"
      )
    )
  )
}

# -----------------------------------------------------------------------------
# Styling
# -----------------------------------------------------------------------------

title_style <- createStyle(
  fontName = "Calibri",
  fontSize = 12,
  textDecoration = "bold",
  halign = "left",
  valign = "center"
)

header_style <- createStyle(
  fontName = "Calibri",
  fontSize = 10,
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  border = "Bottom",
  borderStyle = "thin"
)

row_header_style <- createStyle(
  fontName = "Calibri",
  fontSize = 10,
  textDecoration = "bold",
  halign = "left",
  valign = "center"
)

body_style <- createStyle(
  fontName = "Calibri",
  fontSize = 10,
  halign = "center",
  valign = "center"
)

metric_number_style <- createStyle(
  fontName = "Calibri",
  fontSize = 10,
  halign = "center",
  valign = "center",
  numFmt = paste0("0.", strrep("0", METRIC_DIGITS))
)

integer_number_style <- createStyle(
  fontName = "Calibri",
  fontSize = 10,
  halign = "center",
  valign = "center",
  numFmt = "#,##0"
)

# -----------------------------------------------------------------------------
# Build one resolution-by-return-period table
# -----------------------------------------------------------------------------

.make_metric_table <- function(
    source_dt,
    variable_name,
    metric_name
) {
  
  sub <- source_dt[
    variable == variable_name,
    .(
      resolution,
      return_period,
      value = get(metric_name)
    )
  ]
  
  template <- CJ(
    resolution = RESOLUTION_ORDER,
    return_period = RETURN_PERIODS,
    unique = TRUE
  )
  
  sub <- merge(
    template,
    sub,
    by = c("resolution", "return_period"),
    all.x = TRUE
  )
  
  # Explicit numeric ordering before reshaping
  sub[
    ,
    resolution_order := match(
      resolution,
      RESOLUTION_ORDER
    )
  ]
  
  setorder(
    sub,
    resolution_order,
    return_period
  )
  
  sub[, resolution_order := NULL]
  
  wide <- dcast(
    sub,
    resolution ~ return_period,
    value.var = "value"
  )
  
  # Restore the requested row order after dcast()
  wide[
    ,
    resolution_order := match(
      resolution,
      RESOLUTION_ORDER
    )
  ]
  
  setorder(
    wide,
    resolution_order
  )
  
  wide[, resolution_order := NULL]
  
  setnames(
    wide,
    old = as.character(RETURN_PERIODS),
    new = paste0(RETURN_PERIODS, "-year RP"),
    skip_absent = TRUE
  )
  
  wide[]
}

# -----------------------------------------------------------------------------
# Write one metric workbook
# -----------------------------------------------------------------------------

.write_metric_workbook <- function(metric_name) {
  
  out_file <- file.path(
    OUT_DIR,
    paste0(metric_name, ".xlsx")
  )
  
  if (file_exists(out_file) && !OVERWRITE) {
    .log("Skipping existing workbook: %s", out_file)
    return(invisible(out_file))
  }
  
  wb <- createWorkbook(
    creator = "CISC"
  )
  
  for (variable_name in c("heatindex", "precipitation")) {
    
    sheet_name <- variable_name
    
    addWorksheet(
      wb,
      sheetName = sheet_name,
      gridLines = FALSE,
      tabColour = if (variable_name == "heatindex") {
        "#D55E00"
      } else {
        "#0072B2"
      }
    )
    
    tab <- .make_metric_table(
      source_dt = dt,
      variable_name = variable_name,
      metric_name = metric_name
    )
    
    variable_title <- if (variable_name == "heatindex") {
      "Heat index"
    } else {
      "Precipitation"
    }
    
    writeData(
      wb,
      sheet = sheet_name,
      x = paste0(
        variable_title,
        ": ",
        metric_name,
        " by grid resolution and return period"
      ),
      startRow = 1L,
      startCol = 1L
    )
    
    addStyle(
      wb,
      sheet = sheet_name,
      style = title_style,
      rows = 1L,
      cols = 1L,
      gridExpand = TRUE
    )
    
    writeData(
      wb,
      sheet = sheet_name,
      x = tab,
      startRow = 3L,
      startCol = 1L,
      headerStyle = header_style,
      borders = "rows"
    )
    
    addStyle(
      wb,
      sheet = sheet_name,
      style = row_header_style,
      rows = 4L:(3L + nrow(tab)),
      cols = 1L,
      gridExpand = TRUE,
      stack = TRUE
    )
    
    addStyle(
      wb,
      sheet = sheet_name,
      style = metric_number_style,
      rows = 4L:(3L + nrow(tab)),
      cols = 2L:ncol(tab),
      gridExpand = TRUE,
      stack = TRUE
    )
    
    setColWidths(
      wb,
      sheet = sheet_name,
      cols = 1L,
      widths = 19
    )
    
    setColWidths(
      wb,
      sheet = sheet_name,
      cols = 2L:ncol(tab),
      widths = 14
    )
    
    setRowHeights(
      wb,
      sheet = sheet_name,
      rows = 1L,
      heights = 21
    )
    
    freezePane(
      wb,
      sheet = sheet_name,
      firstActiveRow = 4L,
      firstActiveCol = 2L
    )
    
    addFilter(
      wb,
      sheet = sheet_name,
      rows = 3L,
      cols = 1L:ncol(tab)
    )
  }
  
  saveWorkbook(
    wb,
    file = out_file,
    overwrite = TRUE
  )
  
  .log("Wrote: %s", out_file)
  
  invisible(out_file)
}

# -----------------------------------------------------------------------------
# Write confusion-matrix count workbook
# -----------------------------------------------------------------------------

.write_count_workbook <- function() {
  
  out_file <- file.path(
    OUT_DIR,
    "confusion_matrix_counts.xlsx"
  )
  
  if (file_exists(out_file) && !OVERWRITE) {
    .log("Skipping existing workbook: %s", out_file)
    return(invisible(out_file))
  }
  
  wb <- createWorkbook(
    creator = "CISC"
  )
  
  for (variable_name in c("heatindex", "precipitation")) {
    
    for (count_name in COUNT_METRICS) {
      
      sheet_name <- paste0(
        if (variable_name == "heatindex") "heat" else "precip",
        "_",
        count_name
      )
      
      addWorksheet(
        wb,
        sheetName = sheet_name,
        gridLines = FALSE
      )
      
      tab <- .make_metric_table(
        source_dt = dt,
        variable_name = variable_name,
        metric_name = count_name
      )
      
      variable_title <- if (variable_name == "heatindex") {
        "Heat index"
      } else {
        "Precipitation"
      }
      
      writeData(
        wb,
        sheet = sheet_name,
        x = paste0(
          variable_title,
          ": ",
          count_name,
          " by grid resolution and return period"
        ),
        startRow = 1L,
        startCol = 1L
      )
      
      addStyle(
        wb,
        sheet = sheet_name,
        style = title_style,
        rows = 1L,
        cols = 1L,
        gridExpand = TRUE
      )
      
      writeData(
        wb,
        sheet = sheet_name,
        x = tab,
        startRow = 3L,
        startCol = 1L,
        headerStyle = header_style,
        borders = "rows"
      )
      
      addStyle(
        wb,
        sheet = sheet_name,
        style = row_header_style,
        rows = 4L:(3L + nrow(tab)),
        cols = 1L,
        gridExpand = TRUE,
        stack = TRUE
      )
      
      addStyle(
        wb,
        sheet = sheet_name,
        style = integer_number_style,
        rows = 4L:(3L + nrow(tab)),
        cols = 2L:ncol(tab),
        gridExpand = TRUE,
        stack = TRUE
      )
      
      setColWidths(
        wb,
        sheet = sheet_name,
        cols = 1L,
        widths = 19
      )
      
      setColWidths(
        wb,
        sheet = sheet_name,
        cols = 2L:ncol(tab),
        widths = 14
      )
      
      freezePane(
        wb,
        sheet = sheet_name,
        firstActiveRow = 4L,
        firstActiveCol = 2L
      )
      
      addFilter(
        wb,
        sheet = sheet_name,
        rows = 3L,
        cols = 1L:ncol(tab)
      )
    }
  }
  
  saveWorkbook(
    wb,
    file = out_file,
    overwrite = TRUE
  )
  
  .log("Wrote: %s", out_file)
  
  invisible(out_file)
}

# -----------------------------------------------------------------------------
# Also write a combined long-format workbook for auditing
# -----------------------------------------------------------------------------

.write_all_metrics_workbook <- function() {
  
  out_file <- file.path(
    OUT_DIR,
    "all_performance_metrics.xlsx"
  )
  
  if (file_exists(out_file) && !OVERWRITE) {
    .log("Skipping existing workbook: %s", out_file)
    return(invisible(out_file))
  }
  
  wb <- createWorkbook(
    creator = "CISC"
  )
  
  audit_columns <- c(
    "variable",
    "dataset_tag",
    "resolution",
    "return_period",
    "TP",
    "FP",
    "FN",
    "TN",
    "total",
    "precision",
    "recall",
    "specificity",
    "accuracy",
    "f1",
    "f2"
  )
  
  for (variable_name in c("heatindex", "precipitation")) {
    
    sheet_name <- variable_name
    
    addWorksheet(
      wb,
      sheetName = sheet_name,
      gridLines = FALSE
    )
    
    tab <- copy(
      dt[
        variable == variable_name,
        ..audit_columns
      ]
    )
    
    tab[
      ,
      resolution_order := match(
        resolution,
        RESOLUTION_ORDER
      )
    ]
    
    setorder(
      tab,
      return_period,
      resolution_order
    )
    
    tab[, resolution_order := NULL]
    
    writeDataTable(
      wb,
      sheet = sheet_name,
      x = tab,
      startRow = 1L,
      startCol = 1L,
      tableStyle = "TableStyleMedium2",
      withFilter = TRUE
    )
    
    setColWidths(
      wb,
      sheet = sheet_name,
      cols = 1L:ncol(tab),
      widths = "auto"
    )
    
    freezePane(
      wb,
      sheet = sheet_name,
      firstActiveRow = 2L,
      firstActiveCol = 1L
    )
    
    metric_cols <- match(
      METRICS,
      names(tab)
    )
    
    metric_cols <- metric_cols[
      !is.na(metric_cols)
    ]
    
    addStyle(
      wb,
      sheet = sheet_name,
      style = metric_number_style,
      rows = 2L:(nrow(tab) + 1L),
      cols = metric_cols,
      gridExpand = TRUE,
      stack = TRUE
    )
  }
  
  saveWorkbook(
    wb,
    file = out_file,
    overwrite = TRUE
  )
  
  .log("Wrote: %s", out_file)
  
  invisible(out_file)
}

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

.log("Reading: %s", IN_FILE)
.log("Rows retained: %d", nrow(dt))

for (metric_name in METRICS) {
  .write_metric_workbook(metric_name)
}

.write_count_workbook()
.write_all_metrics_workbook()

.log(
  "Completed performance tables in: %s",
  OUT_DIR
)