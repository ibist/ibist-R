# Export the prepared R datasets and their metadata for the Python package.
# Run from the R package root after the relevant data-raw preparation scripts:
# Rscript data-raw/export-python.R [../ibist-py/src/ibist/_data]
# This step never changes the R data files or duplicates their cleaning rules.

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Install jsonlite to run the development data exporter.")
}

args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[[1L]] else "../ibist-py/src/ibist/_data"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
dataset_names <- c("asah", "depress", "fluorosis", "midbp", "mihdl",
                   "nrs", "rds")
package_version <- read.dcf("DESCRIPTION", fields = "Version")[[1L]]
revision <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
status_args <- c("status", "--porcelain", "--untracked-files=normal",
                 "--", "data", "data-raw", "R/data-*.R", "man")
status <- system2("git", status_args, stdout = TRUE)
dirty <- length(status) > 0L

rd_text <- function(x) {
  if (identical(attr(x, "Rd_tag"), "USERMACRO")) {
    prefix <- if (identical(attr(x, "macro"), "\\doi")) {
      "https://doi.org/"
    } else ""
    return(paste0(prefix, paste(x[-1L], collapse = "")))
  }
  if (identical(attr(x, "Rd_tag"), "\\Sexpr")) return("")
  if (is.list(x)) return(paste(vapply(x, rd_text, ""), collapse = ""))
  paste(x, collapse = "")
}

plain_text <- function(x) {
  trimws(gsub("[[:space:]]+", " ", rd_text(x)))
}

rd_nodes <- function(x, tag) {
  found <- list()
  if (is.list(x)) {
    if (identical(attr(x, "Rd_tag"), tag)) found <- list(x)
    for (child in x) found <- c(found, rd_nodes(child, tag))
  }
  found
}

hash <- function(path) unname(tools::md5sum(path))
raw_inputs <- list(
  asah = character(),
  depress = c("data-raw/kalmbach2018.zip", "data-raw/S1 dataset.sav"),
  fluorosis = character(),
  midbp = "data-raw/midbp.csv",
  mihdl = "data-raw/mihdl.csv",
  nrs = "data-raw/nrs.csv",
  rds = "data-raw/rds.tab"
)
units <- list(depress = list(age = "years"), nrs = list(age = "years"),
              midbp = list(pyears = "person-years"),
              mihdl = list(pyears = "person-years"),
              rds = list(bwt = "grams"))

for (name in dataset_names) {
  prepared <- file.path("data", paste0(name, ".rda"))
  script <- file.path("data-raw", paste0(name, ".R"))
  documentation <- file.path("man", paste0(name, ".Rd"))
  env <- new.env(parent = baseenv())
  loaded <- load(prepared, envir = env)
  stopifnot(identical(loaded, name), is.data.frame(env[[name]]))
  dat <- env[[name]]
  rd <- tools::parse_Rd(documentation)
  items <- Filter(function(x) length(x) == 2L,
                  rd_nodes(rd_nodes(rd, "\\format"), "\\item"))
  descriptions <- setNames(lapply(items, function(x) plain_text(x[[2L]])),
                            vapply(items, function(x) plain_text(x[[1L]]), ""))
  columns <- lapply(names(dat), function(column_name) {
    x <- dat[[column_name]]
    value_labels <- attr(x, "labels")
    list(name = column_name,
         type = if (is.factor(x)) "category" else {
           switch(typeof(x), integer = "integer", double = "number",
                  character = "string", logical = "boolean",
                  stop("Unsupported column type in ", name, ": ", column_name))
         },
         r_class = I(class(x)), description = descriptions[[column_name]],
         label = attr(x, "label"), units = units[[name]][[column_name]],
         levels = if (is.factor(x)) I(levels(x)) else NULL,
         ordered = is.ordered(x), missing_count = sum(is.na(x)),
         value_labels = unname(lapply(seq_along(value_labels), function(i) {
           list(value = unname(value_labels[[i]]),
                label = names(value_labels)[i])
         })))
  })
  values <- lapply(dat, function(x) {
    if (is.factor(x)) as.character(x) else {
      attributes(x) <- NULL
      x
    }
  })
  data_path <- file.path(output, paste0(name, ".json.gz"))
  con <- gzfile(data_path, "wt", encoding = "UTF-8")
  writeLines(jsonlite::toJSON(values, auto_unbox = FALSE, digits = NA,
                              na = "null", null = "null"), con, useBytes = TRUE)
  close(con)
  con <- gzfile(data_path, "rb")
  roundtrip <- jsonlite::fromJSON(con, simplifyVector = TRUE)
  close(con)
  if (!isTRUE(all.equal(roundtrip, values, check.attributes = FALSE,
                        tolerance = 0))) {
    stop("JSON round-trip changed data values in ", name)
  }
  metadata <- list(
    schema_version = 1L, name = name,
    title = plain_text(rd_nodes(rd, "\\title")[[1L]]),
    description = plain_text(rd_nodes(rd, "\\description")[[1L]]),
    source = lapply(rd_nodes(rd, "\\source"), plain_text),
    notes = lapply(rd_nodes(rd, "\\details"), plain_text),
    nrows = nrow(dat), ncols = ncol(dat), columns = columns,
    provenance = list(
      repository = "https://github.com/ibist/ibist-R", r_revision = revision,
      r_version = as.character(getRversion()),
      package_version = package_version,
      uncommitted_preparation_changes = dirty,
      prepared_file = prepared, prepared_md5 = hash(prepared),
      preparation_script = script, preparation_md5 = hash(script),
      documentation_file = documentation,
      documentation_md5 = hash(documentation),
      exporter = "data-raw/export-python.R",
      exporter_md5 = hash("data-raw/export-python.R"),
      raw_inputs = lapply(raw_inputs[[name]], function(path) {
        list(path = path, available = file.exists(path),
             md5 = if (file.exists(path)) hash(path) else NULL)
      }),
      data_md5 = hash(data_path)
    )
  )
  jsonlite::write_json(metadata, file.path(output, paste0(name, ".json")),
                       auto_unbox = TRUE, pretty = TRUE,
                       null = "null", digits = NA)
  message("Exported ", name, ": ", nrow(dat), " rows, ", ncol(dat), " columns")
}

jsonlite::write_json(list(schema_version = 1L, datasets = I(dataset_names)),
                     file.path(output, "index.json"), auto_unbox = TRUE,
                     pretty = TRUE)
