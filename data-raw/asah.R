## data-raw/asah.R

# Reconstruct the paired ratings from the published aggregate tables.
# Each scale has 103 paired assessments. The generated assessment numbers
# identify expanded rows only; they are not participant identifiers, and the
# rows do not link assessments across scales.
asah_tables <- list(
  WFNS = matrix(
    c(
      37, 9, 0, 0, 0,
       9, 27, 1, 2, 0,
       0, 2, 3, 0, 0,
       0, 3, 1, 4, 0,
       0, 0, 0, 1, 4
    ),
    nrow = 5,
    byrow = TRUE
  ),
  HH = matrix(
    c(
      18, 4, 2, 0, 0,
      12, 32, 3, 0, 0,
       4, 9, 8, 1, 0,
       0, 0, 2, 4, 1,
       0, 0, 0, 1, 2
    ),
    nrow = 5,
    byrow = TRUE
  ),
  PAASH = matrix(
    c(
      37, 8, 0, 0, 0,
      10, 40, 1, 0, 0,
       0, 1, 0, 0, 0,
       0, 0, 2, 2, 0,
       0, 0, 0, 0, 2
    ),
    nrow = 5,
    byrow = TRUE
  )
)

asah <- do.call(
  rbind,
  lapply(names(asah_tables), function(scale_name) {
    tab <- asah_tables[[scale_name]]
    cells <- which(tab > 0, arr.ind = TRUE)
    counts <- tab[cells]

    data.frame(
      scale = scale_name,
      assessment = seq_len(sum(counts)),
      observer1 = factor(
        rep(cells[, 1], counts),
        levels = seq_len(nrow(tab)),
        ordered = TRUE
      ),
      observer2 = factor(
        rep(cells[, 2], counts),
        levels = seq_len(ncol(tab)),
        ordered = TRUE
      )
    )
  })
)

rownames(asah) <- NULL
asah$scale <- factor(asah$scale, levels = names(asah_tables))

save(asah, file = "data/asah.rda", compress = "xz")
