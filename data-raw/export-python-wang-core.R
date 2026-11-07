# Export the authoritative, language-neutral Wang core. Run from ibist-R.
# --check verifies synchronization without writing either file.
source <- "src/wang_core.h"
output <- "../ibist-py/src/ibist/_native/wang_core.h"
metadata <- "../ibist-py/src/ibist/_native/wang_core.json"
bytes <- readBin(source, "raw", file.info(source)$size)
if ("--check" %in% commandArgs(trailingOnly = TRUE)) {
  if (!file.exists(output) || !identical(bytes, readBin(
    output, "raw", file.info(output)$size))) stop("Wang core copies differ.")
  manifest <- jsonlite::read_json(metadata)
  if (manifest$md5 != unname(tools::md5sum(source))) {
    stop("Wang core provenance hash differs.")
  }
  message("Wang core copies and source hash match.")
} else {
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  writeBin(bytes, output)
  jsonlite::write_json(list(authoritative_source = "ibist-R/src/wang_core.h",
                           md5 = unname(tools::md5sum(source)),
                           license = "GPL-3.0-or-later"), metadata,
                       pretty = TRUE, auto_unbox = TRUE)
  message("Exported Wang core and provenance.")
}
