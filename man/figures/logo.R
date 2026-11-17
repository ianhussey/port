# This script draws man/figures/logo.png; run from the package root. The text is
# set in Space Grotesk, an open font. Before drawing, install Space Grotesk, or
# the text falls back to a generic sans.
#
# Containers stacked on the quay of a port. They hold a correlation table: the
# containers are stacked as the lower triangle of a 4 x 4 matrix, the diagonal
# of 1s forming the steps, as in the package's earlier logo. The table is one
# that check_corr_psd() certifies as consistent, so the cargo is fit to ship.

sky <- "#cfe3f1"
water <- "#1d5a8e"
border <- "#0f3a5f"
hull <- "#23272e"
quay <- "#c8b08c"
coping <- "#8c6e4f"
cream <- "#f4ede0"
ink <- "#23272e"
cargo <- c("#e8955a", "#c8553d", "#d9a43c")

r <- matrix(
  c(
    1,
    .42,
    .31,
    -.17,
    .42,
    1,
    .58,
    -.08,
    .31,
    .58,
    1,
    .24,
    -.17,
    -.08,
    .24,
    1
  ),
  4
)

pts <- function(x, y) paste(sprintf("%.1f,%.1f", x, y), collapse = " ")
hexs <- pts(
  c(260, 511.1, 511.1, 260, 8.9, 8.9),
  c(10, 155, 445, 590, 445, 155)
)

# The top of the quay, and the waterline a little below it
deck <- 330
sea <- 362

# One container, its top left corner at (x, y), with ribs and a label. The
# diagonal is cream; the others take turns in the warm colors.
cw <- 66
ch <- 44
container <- function(i, j) {
  x <- 128 + (j - 1) * cw
  y <- deck - (5 - i) * ch
  fill <- if (i == j) cream else cargo[(i + j) %% 3 + 1]
  label <- if (i == j) "1" else sub("0\\.", ".", sprintf("%.2f", r[i, j]))
  ribs <- x + seq(8, cw - 8, length.out = 6)
  c(
    sprintf(
      paste(
        '<rect x="%.1f" y="%.1f" width="%g" height="%g" fill="%s"',
        'stroke="%s" stroke-width="2"/>'
      ),
      x + 1,
      y + 1,
      cw - 2,
      ch - 2,
      fill,
      ink
    ),
    sprintf(
      paste(
        '<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s"',
        'stroke-opacity="0.18" stroke-width="2"/>'
      ),
      ribs,
      y + 5,
      ribs,
      y + ch - 5,
      ink
    ),
    sprintf(
      paste(
        '<text x="%.1f" y="%.1f" text-anchor="middle"',
        'font-family="Space Grotesk, sans-serif" font-weight="bold"',
        'font-size="21" fill="%s">%s</text>'
      ),
      x + cw / 2,
      y + ch / 2 + 7.5,
      ink,
      label
    )
  )
}
cells <- which(lower.tri(r, diag = TRUE), arr.ind = TRUE)
containers <- unlist(Map(container, cells[, "row"], cells[, "col"]))

# A wavy band of water, from y down, that the quay wall sinks into
waves <- function(y, amp, len, fill) {
  x <- seq(0, 520, by = len / 2)
  ys <- y + amp * rep_len(c(-1, 1), length(x))
  path <- paste0(
    sprintf("M 0 %.1f ", y),
    paste(
      sprintf("Q %.1f %.1f %.1f %.1f", x[-1] - len / 4, ys[-1], x[-1], y),
      collapse = " "
    ),
    " L 520 600 L 0 600 Z"
  )
  sprintf('<path d="%s" fill="%s"/>', path, fill)
}

# Two gulls
gull <- function(x, y, s) {
  sprintf(
    paste(
      '<path d="M %.1f %.1f q %g %g %g 0 q %g %g %g 0" fill="none"',
      'stroke="%s" stroke-width="3" stroke-linecap="round"/>'
    ),
    x - s,
    y,
    s / 2,
    -s / 2,
    s,
    s / 2,
    -s / 2,
    s,
    ink
  )
}

font_size <- 84

svg <- c(
  paste(
    '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 520 600"',
    'width="520" height="600">'
  ),
  sprintf('<clipPath id="hex"><polygon points="%s"/></clipPath>', hexs),
  '<g clip-path="url(#hex)">',
  sprintf('<rect width="520" height="600" fill="%s"/>', sky),
  gull(350, 160, 14),
  gull(390, 138, 10),
  containers,
  # the quay: a stone wall into the water, with a darker coping along the
  # top and rubber fenders on its face
  sprintf(
    '<rect x="0" y="%g" width="456" height="%g" fill="%s"/>',
    deck,
    sea + 40 - deck,
    quay
  ),
  sprintf(
    '<rect x="0" y="%g" width="456" height="8" fill="%s"/>',
    deck,
    coping
  ),
  sprintf(
    '<rect x="%g" y="%g" width="12" height="%g" rx="4" fill="%s"/>',
    c(60, 250, 420),
    deck + 12,
    sea - deck,
    hull
  ),
  waves(sea + 8, 5, 52, water),
  # the baseline sits below the hex's lower corners (y = 445) by half the
  # height of the letters, so the word is centered on their level
  sprintf(
    paste(
      '<text x="260" y="%.1f" text-anchor="middle"',
      'font-family="Space Grotesk, sans-serif" font-weight="bold"',
      'font-size="%g" fill="%s">port</text>'
    ),
    445 + 0.349 * font_size,
    font_size,
    cream
  ),
  "</g>",
  sprintf(
    paste(
      '<polygon points="%s" fill="none" stroke="%s" stroke-width="16"',
      'stroke-linejoin="round"/>'
    ),
    hexs,
    border
  ),
  "</svg>"
)

# Write logo to disk
svg_path <- tempfile(fileext = ".svg")
writeLines(svg, svg_path)
rsvg::rsvg_png(svg_path, "man/figures/logo.png", width = 1040)
