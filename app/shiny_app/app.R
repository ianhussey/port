# port: a browser front end ---------------------------------------------------
#
# Given a reported correlation matrix whose off-diagonals have been rounded to a
# fixed number of decimals, decide whether ANY positive semidefinite (PSD)
# matrix lies within the induced rounding box; when none does, localize and
# grade the inconsistency. All computation is done by the port package; this app
# only handles input, layout, and download.
#
# Hussey, I. (2026). port: Positive-semi-definiteness Of Rounded correlation
# Tables. https://doi.org/10.5281/zenodo.21434236

library(shiny)
library(bslib)
library(port)

port_blue <- "#1d5a8e"
port_cream <- "#f4ede0"

port_citation <- paste0(
  "Hussey, I. (2026). port: Positive-semi-definiteness Of Rounded correlation ",
  "Tables. [Computer software] https://github.com/ianhussey/port ",
  "doi:10.5281/zenodo.21434236"
)

# ---- input handling ---------------------------------------------------------

# Read an uploaded csv/xlsx into a data frame, without mangling column names.
read_uploaded <- function(path, name) {
  ext <- tolower(tools::file_ext(name))
  if (ext %in% c("csv", "txt")) {
    utils::read.csv(
      path,
      header = TRUE,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      na.strings = c("NA", "na", "", ".", "-")
    )
  } else if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("Reading Excel files requires the 'readxl' package.", call. = FALSE)
    }
    as.data.frame(readxl::read_excel(path), check.names = FALSE)
  } else {
    stop(
      "Unsupported file type '.", ext, "'. Upload a .csv or .xlsx file.",
      call. = FALSE
    )
  }
}

# Coerce a data frame into a symmetric correlation matrix with unit diagonal.
# Accepts a full matrix or either triangle (the empty half left blank / NA), and
# an optional leading column of variable names. Off-diagonal cells that are
# missing in BOTH triangles stay NA: port frees those to [-1, 1], so verdicts
# then hold whatever the unreported values were.
parse_corr_matrix <- function(df) {
  df <- as.data.frame(df, check.names = FALSE)
  if (nrow(df) == 0L || ncol(df) == 0L) {
    stop("The file contains no data.", call. = FALSE)
  }

  # A leading column of variable names, if there is one. An all-empty first
  # column is the omitted half of a triangular table, not a label column.
  labels <- NULL
  first <- df[[1L]]
  if ((is.character(first) || is.factor(first)) && !all(is.na(first))) {
    labels <- trimws(as.character(first))
    df <- df[, -1L, drop = FALSE]
  }

  p <- ncol(df)
  if (p < 2L) {
    stop("At least two variables are needed.", call. = FALSE)
  }
  if (nrow(df) != p) {
    stop(
      sprintf(
        paste0(
          "Expected a square table: found %d rows and %d correlation columns. ",
          "Supply one row per variable, in the same order as the columns."
        ),
        nrow(df),
        p
      ),
      call. = FALSE
    )
  }

  m <- matrix(NA_real_, p, p)
  for (k in seq_len(p)) {
    col <- df[[k]]
    num <- suppressWarnings(as.numeric(as.character(col)))
    bad <- !is.na(col) & is.na(num) & trimws(as.character(col)) != ""
    if (any(bad)) {
      stop(
        sprintf(
          "Column %d contains a value that is not a number: '%s'.",
          k,
          as.character(col)[which(bad)[1L]]
        ),
        call. = FALSE
      )
    }
    m[, k] <- num
  }

  # Fold the two triangles together.
  R <- matrix(NA_real_, p, p)
  diag(R) <- 1
  for (i in seq_len(p - 1L)) {
    for (j in seq(i + 1L, p)) {
      a <- m[i, j]
      b <- m[j, i]
      v <- if (is.na(a) && is.na(b)) {
        NA_real_
      } else if (is.na(a)) {
        b
      } else if (is.na(b)) {
        a
      } else {
        if (abs(a - b) > 1e-8) {
          stop(
            sprintf(
              paste0(
                "Cell (%d,%d) disagrees across the two triangles: %s vs %s. ",
                "Supply one triangle, or a symmetric matrix."
              ),
              i,
              j,
              format(a),
              format(b)
            ),
            call. = FALSE
          )
        }
        (a + b) / 2
      }
      R[i, j] <- R[j, i] <- v
    }
  }

  off <- R[upper.tri(R)]
  if (all(is.na(off))) {
    stop("No correlations were found off the diagonal.", call. = FALSE)
  }
  if (any(abs(off) > 1, na.rm = TRUE)) {
    stop(
      "Some entries lie outside [-1, 1]; these are not correlations.",
      call. = FALSE
    )
  }

  # Variable names: the column headers, unless they are placeholders.
  nm <- colnames(df)
  placeholder <- is.null(nm) ||
    any(is.na(nm) | trimws(nm) == "") ||
    all(grepl("^(V|X)[0-9]+$", nm))
  if (placeholder) {
    nm <- if (!is.null(labels) && length(labels) == p) {
      labels
    } else {
      paste0("V", seq_len(p))
    }
  }
  dimnames(R) <- list(make.unique(nm), make.unique(nm))
  R
}

# ---- presentation -----------------------------------------------------------

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "—", formatC(x, digits = digits, format = "g"))
}

# The matrix as a display-ready data frame (blank upper triangle: it is a mirror
# of the lower one, and the redundancy makes larger tables harder to scan).
matrix_df <- function(R) {
  out <- matrix("", nrow(R), ncol(R))
  for (i in seq_len(nrow(R))) {
    for (j in seq_len(ncol(R))) {
      out[i, j] <- if (j > i) {
        ""
      } else if (i == j) {
        "1"
      } else if (is.na(R[i, j])) {
        "NA"
      } else {
        format(R[i, j])
      }
    }
  }
  df <- as.data.frame(out, stringsAsFactors = FALSE)
  colnames(df) <- colnames(R)
  cbind(variable = rownames(R), df)
}

verdict_headline <- function(chk) {
  switch(
    chk$verdict,
    inconsistent = "INCONSISTENT",
    consistent = if (isTRUE(chk$certified)) {
      "CONSISTENT (certified)"
    } else {
      "CONSISTENT (presumed)"
    },
    undecided = "UNDECIDED",
    toupper(chk$verdict)
  )
}

verdict_gloss <- function(chk) {
  switch(
    chk$verdict,
    inconsistent = paste(
      "No positive semidefinite matrix lies within the rounding box. The",
      "reported values cannot be the complete-data, single-sample Pearson",
      "correlation matrix of any dataset, up to the stated rounding."
    ),
    consistent = if (isTRUE(chk$certified)) {
      paste(
        "An in-box positive semidefinite matrix was exhibited and verified",
        "(Rump 2006): a rigorous existence certificate."
      )
    } else {
      paste(
        "Presumed consistent: the witness margin is comfortably positive, but",
        "no in-box PSD matrix was exhibited. This means only 'not shown",
        "inconsistent'."
      )
    },
    undecided = paste(
      "The witness margin sits in the precision-limited ambiguous zone and no",
      "in-box PSD matrix could be certified, including by the",
      "alternating-projection escalation."
    ),
    ""
  )
}

verdict_banner <- function(chk) {
  cls <- switch(
    chk$verdict,
    inconsistent = "port-bad",
    consistent = "port-good",
    "port-warn"
  )
  div(
    class = paste("port-banner", cls),
    div(class = "port-verdict", verdict_headline(chk)),
    div(class = "port-gloss", verdict_gloss(chk))
  )
}

stat_row <- function(chk, excusable = NA_real_) {
  cells <- list(
    c("variables", as.character(chk$p)),
    c("rounding half-width δ", fmt(chk$delta)),
    c("tier", chk$tier),
    c("excusable δ", fmt(excusable))
  )
  div(
    class = "port-stats",
    lapply(cells, function(x) {
      div(class = "port-stat", div(class = "k", x[1]), div(class = "v", x[2]))
    })
  )
}

# Feasible interval per cell, against the reported point and its rounding box.
# A row with no feasible interval at all (the cell cannot be rescued on its own)
# is shaded rather than labelled, so the reported point stays readable.
plot_intervals <- function(iv, delta, labels = NULL) {
  n <- nrow(iv)
  lab <- if (!is.null(labels)) {
    paste(labels[iv$i], labels[iv$j], sep = " – ")
  } else {
    sprintf("(%d,%d)", iv$i, iv$j)
  }
  y <- rev(seq_len(n))
  red <- "#b3261e"

  op <- par(mar = c(4, max(8, min(22, max(nchar(lab)) * 0.55)), 3.5, 1))
  on.exit(par(op), add = TRUE)

  plot(
    NA,
    xlim = c(-1, 1),
    ylim = c(0.5, n + 0.5),
    yaxt = "n",
    xlab = "correlation",
    ylab = "",
    bty = "n"
  )
  abline(v = seq(-1, 1, 0.5), col = "grey90")
  abline(v = 0, col = "grey70")
  axis(2, at = y, labels = lab, las = 1, tick = FALSE, cex.axis = 0.85)

  for (k in seq_len(n)) {
    feasible <- identical(iv$status[k], "feasible")
    if (feasible) {
      segments(iv$lo[k], y[k], iv$hi[k], y[k], col = port_blue, lwd = 6)
    } else {
      rect(-1, y[k] - 0.4, 1, y[k] + 0.4, col = "#fdeeed", border = NA)
    }
    r <- iv$reported[k]
    if (!is.na(r)) {
      inside <- feasible && r >= iv$lo[k] && r <= iv$hi[k]
      segments(r - delta, y[k], r + delta, y[k], col = "black", lwd = 2)
      points(
        r,
        y[k],
        pch = 21,
        cex = 1.3,
        bg = if (inside) "white" else red,
        col = "black"
      )
    }
  }

  # Place the legend in the top margin, a fixed number of text lines above the
  # plot region, so it clears the rows whatever the plot's height.
  usr <- par("usr")
  line_h <- diff(usr[3:4]) / par("pin")[2] * par("csi")
  legend(
    x = mean(usr[1:2]),
    y = usr[4] + 2 * line_h,
    xjust = 0.5,
    yjust = 0.5,
    horiz = TRUE,
    bty = "n",
    xpd = NA,
    cex = 0.8,
    legend = c("feasible interval", "reported (in range)", "reported (out of range)"),
    pch = c(NA, 21, 21),
    pt.bg = c(NA, "white", red),
    lwd = c(6, NA, NA),
    col = c(port_blue, "black", "black")
  )
  invisible(NULL)
}

# ---- example data -----------------------------------------------------------

example_path <- "example_correlation_matrix.csv"

# ---- UI ---------------------------------------------------------------------

port_theme <- bs_theme(
  version = 5,
  primary = port_blue,
  bg = "#ffffff",
  fg = "#1c1c1c",
  base_font = font_link(
    "Space Grotesk",
    href = paste0(
      "https://fonts.googleapis.com/css2?",
      "family=Space+Grotesk:wght@400;500;600;700&display=swap"
    )
  )
)

port_css <- sprintf(
  "
  .port-banner { border-radius: 10px; padding: 1rem 1.25rem; margin: 0 0 1rem 0;
                 border-left: 10px solid; }
  .port-banner .port-verdict { font-size: 1.5rem; font-weight: 700;
                               letter-spacing: 0.02em; }
  .port-banner .port-gloss { margin-top: 0.35rem; font-size: 0.95rem; }
  .port-good { background: #eef6ee; border-color: #2e7d32; }
  .port-good .port-verdict { color: #2e7d32; }
  .port-bad  { background: #fdeeed; border-color: #b3261e; }
  .port-bad .port-verdict { color: #b3261e; }
  .port-warn { background: #fdf6e6; border-color: #b26a00; }
  .port-warn .port-verdict { color: #b26a00; }
  .port-stats { display: flex; flex-wrap: wrap; gap: 0.75rem; margin-bottom: 1rem; }
  .port-stat { background: %s; border-radius: 8px; padding: 0.5rem 0.9rem;
               min-width: 8rem; }
  .port-stat .k { font-size: 0.75rem; text-transform: uppercase;
                  letter-spacing: 0.06em; color: #5b5b5b; }
  .port-stat .v { font-size: 1.15rem; font-weight: 600; color: %s; }
  .port-note { font-size: 0.9rem; color: #5b5b5b; }
  .navbar-brand img { margin-right: 0.5rem; }
  .shiny-table-output { overflow-x: auto; }
  .shiny-table-output table td, .shiny-table-output table th {
    white-space: nowrap;
  }
  pre { background: #f7f7f7; border-radius: 8px; padding: 0.9rem;
        overflow-x: auto; }
  ",
  port_cream,
  port_blue
)

matrix_panel <- nav_panel(
  title = "Correlation matrix",
  layout_sidebar(
    fillable = FALSE,
    sidebar = sidebar(
      width = 330,
      fileInput(
        "file",
        "Correlation matrix (.csv or .xlsx)",
        accept = c(".csv", ".xlsx", ".xls")
      ),
      downloadLink("download_example", "Download an example file"),
      tags$hr(),
      numericInput(
        "decimals",
        "Decimal places the correlations were rounded to",
        value = 2,
        min = 0,
        max = 8,
        step = 1
      ),
      selectInput(
        "rounding",
        "Rounding rule used by the authors",
        choices = c(
          "nearest (conventional)" = "nearest",
          "truncate (toward zero)" = "truncate",
          "floor" = "floor",
          "ceiling" = "ceiling"
        )
      ),
      tags$p(
        class = "port-note",
        "Using the rule that actually generated the reported values is a",
        "soundness requirement: a mismatched box can exclude the true value."
      ),
      tags$hr(),
      actionButton("localize", "Localize the fault", class = "btn-primary"),
      tags$p(
        class = "port-note",
        "Fault localization is the slow step (seconds to tens of seconds)."
      ),
      tags$hr(),
      uiOutput("download_ui")
    ),
    uiOutput("matrix_status"),
    navset_underline(
      nav_panel(
        "Verdict",
        uiOutput("verdict_ui"),
        verbatimTextOutput("check_print")
      ),
      nav_panel(
        "Matrix",
        tags$p(
          class = "port-note",
          "As parsed. The upper triangle is left blank (it mirrors the lower);",
          "NA marks a cell that was not reported, which port frees to [-1, 1]."
        ),
        tableOutput("matrix_table")
      ),
      nav_panel(
        "Fault localization",
        uiOutput("fault_ui"),
        verbatimTextOutput("fault_print")
      ),
      nav_panel(
        "Cell intervals",
        tags$p(
          class = "port-note",
          "The interval each cell could occupy for the matrix to be a valid",
          "correlation matrix, and the edit beyond rounding that the reported",
          "value would need. Bars are the feasible interval; the point and",
          "whisker are the reported value and its rounding box."
        ),
        checkboxInput(
          "hold_box",
          "Let the other cells roam their rounding boxes (rigorous, slower)",
          value = FALSE
        ),
        uiOutput("interval_plot_note"),
        plotOutput("interval_plot", height = "auto"),
        tableOutput("interval_table")
      ),
      nav_panel(
        "Inconsistent triples",
        tags$p(
          class = "port-note",
          "Triples of variables whose 3x3 submatrix cannot be positive",
          "semidefinite for any values within the rounding boxes. A sound",
          "interval bound on the determinant: no solver involved."
        ),
        tableOutput("triples_table")
      )
    )
  )
)

trio_panel <- nav_panel(
  title = "Three correlations",
  layout_sidebar(
    fillable = FALSE,
    sidebar = sidebar(
      width = 330,
      tags$p(
        class = "port-note",
        "Three variables A, B and C measured in one sample. Any two of the",
        "correlations constrain the third."
      ),
      numericInput("rAB", "r(A, B)", value = 0.9, min = -1, max = 1, step = 0.05),
      numericInput("rAC", "r(A, C)", value = 0.9, min = -1, max = 1, step = 0.05),
      numericInput("rBC", "r(B, C)", value = -0.9, min = -1, max = 1, step = 0.05),
      tags$hr(),
      numericInput(
        "trio_decimals",
        "Decimal places they were rounded to",
        value = 2,
        min = 0,
        max = 8,
        step = 1
      )
    ),
    uiOutput("trio_verdict_ui"),
    tags$p(
      class = "port-note",
      "Each correlation's feasible interval given the other two, letting",
      "those roam their rounding boxes."
    ),
    plotOutput("trio_plot", height = "260px"),
    tableOutput("trio_table"),
    verbatimTextOutput("trio_print")
  )
)

about_panel <- nav_panel(
  title = "About",
  div(
    class = "container",
    style = "max-width: 52rem;",
    tags$h3("Positive-semi-definiteness Of Rounded correlation Tables"),
    tags$p(
      "port is a forensic meta-science tool for auditing reported Pearson",
      "correlation matrices. Published correlations are rounded, so each true",
      "value lies somewhere within a small range around the printed number.",
      "port asks whether any valid correlation matrix — one that is",
      tags$em("positive semidefinite"),
      "— can lie within all of those ranges at once."
    ),
    tags$h4("The method"),
    tags$p(HTML(paste0(
      "For a reported matrix R rounded to d decimals, let ",
      "&delta; = 0.5&middot;10<sup>&minus;d</sup>. For any test direction v, ",
      "the most PSD-favourable in-box matrix along v attains, exactly, ",
      "<code>v'Rv + &delta;(&#8214;v&#8214;&#8321;&sup2; &minus; ",
      "&#8214;v&#8214;&#8322;&sup2;)</code>. If that maximum is negative, no ",
      "in-box matrix is PSD and v certifies inconsistency. The evaluation is ",
      "bounded by a rigorous a-priori floating-point error model ",
      "(Higham 2002), so an inconsistent verdict does not depend on any ",
      "solver tolerance. A consistent verdict is certified by exhibiting an ",
      "in-box matrix and verifying it PSD with a guaranteed Cholesky test ",
      "(Rump 2006)."
    ))),
    tags$h4("What an inconsistent verdict does and does not prove"),
    tags$p(
      "It proves that the reported values cannot be the complete-data,",
      "single-sample Pearson correlation matrix of any dataset, up to the",
      "stated rounding. Getting from there to a reporting error means ruling",
      "out legitimate generators of non-PSD reported matrices:"
    ),
    tags$ul(
      tags$li(
        tags$b("pairwise deletion:"),
        "each cell computed on different available cases;"
      ),
      tags$li(
        tags$b("polychoric / tetrachoric estimation:"),
        "estimated cell by cell;"
      ),
      tags$li(
        tags$b("meta-analytic assembly:"),
        "each cell pooled from a different set of studies;"
      ),
      tags$li(
        tags$b("upstream disattenuation:"),
        "reliability-corrected matrices can legitimately leave the PSD cone."
      )
    ),
    tags$p(
      "The severity and excusable-δ statistics are the discriminator: those",
      "mechanisms typically produce small violations, so an inconsistency",
      "requiring corrections many times the rounding step survives all of",
      "these excuses, while a near-boundary violation does not."
    ),
    tags$h4("Citation"),
    tags$p(port_citation),
    tags$h4("Links"),
    tags$ul(
      tags$li(a(
        "Source code and documentation on GitHub",
        href = "https://github.com/ianhussey/port",
        target = "_blank",
        rel = "noopener"
      )),
      tags$li(a(
        "Archived release (doi:10.5281/zenodo.21434236)",
        href = "https://doi.org/10.5281/zenodo.21434236",
        target = "_blank",
        rel = "noopener"
      )),
      tags$li(a(
        "Report a bug",
        href = "https://github.com/ianhussey/port/issues",
        target = "_blank",
        rel = "noopener"
      ))
    ),
    tags$p(
      class = "port-note",
      sprintf("This app runs port version %s.", utils::packageVersion("port"))
    ),
    tags$p(
      class = "port-note",
      "Code is MIT licensed © Ian Hussey. Text and images are CC BY 4.0."
    )
  )
)

ui <- page_navbar(
  title = tags$span(
    tags$img(src = "logo.png", height = 34, alt = "port logo"),
    "port"
  ),
  window_title = "port",
  theme = port_theme,
  header = tags$head(tags$style(HTML(port_css))),
  trio_panel,
  matrix_panel,
  about_panel
)

# ---- server -----------------------------------------------------------------

server <- function(input, output, session) {
  decimals <- reactive({
    d <- input$decimals
    validate(need(
      is.numeric(d) && !is.na(d) && d >= 0 && d <= 8,
      "Decimal places must be between 0 and 8."
    ))
    as.integer(d)
  })

  # -- panel 1: a full matrix -------------------------------------------------

  parsed <- reactive({
    req(input$file)
    tryCatch(
      list(R = parse_corr_matrix(read_uploaded(
        input$file$datapath,
        input$file$name
      ))),
      error = function(e) list(error = conditionMessage(e))
    )
  })

  R_mat <- reactive({
    p <- parsed()
    validate(need(is.null(p$error), p$error))
    p$R
  })

  output$matrix_status <- renderUI({
    if (is.null(input$file)) {
      return(div(
        class = "port-note",
        "Upload a correlation matrix to begin, or download the example file.",
        "One row and one column per variable; either triangle, or the full",
        "matrix; an optional leading column of variable names."
      ))
    }
    p <- parsed()
    if (!is.null(p$error)) {
      return(div(class = "port-banner port-bad", div(
        class = "port-verdict",
        "Could not read the file"
      ), div(class = "port-gloss", p$error)))
    }
    NULL
  })

  check <- reactive({
    R <- R_mat()
    check_corr_psd(R, decimals = decimals(), rounding = input$rounding)
  })

  excusable <- reactive({
    R <- R_mat()
    if (anyNA(R)) {
      return(NA_real_)
    }
    withProgress(
      message = "Locating the excusable rounding half-width",
      value = 0.5,
      excusable_delta(R)
    )
  })

  output$verdict_ui <- renderUI({
    chk <- check()
    tagList(verdict_banner(chk), stat_row(chk, excusable()))
  })

  output$check_print <- renderPrint({
    print(check())
  })

  output$matrix_table <- renderTable(
    {
      matrix_df(R_mat())
    },
    striped = TRUE,
    na = "",
    rownames = FALSE
  )

  fault <- eventReactive(input$localize, {
    R <- R_mat()
    withProgress(
      message = "Localizing the fault",
      detail = "per-cell intervals, triples, leave-one-out, sparse correction",
      value = 0.5,
      tryCatch(
        localize_psd_fault(R, decimals = decimals()),
        error = function(e) list(error = conditionMessage(e))
      )
    )
  })

  output$fault_ui <- renderUI({
    if (input$localize == 0L) {
      return(div(
        class = "port-note",
        "Press 'Localize the fault' in the sidebar. This asks where the",
        "non-PSDness comes from and how severe it is: the per-cell intervals,",
        "inconsistent triples, variables whose removal restores feasibility,",
        "the smallest joint correction, and whether the cause looks benign",
        "(a structural dependency near the boundary) or substantive."
      ))
    }
    f <- fault()
    if (!is.null(f$error)) {
      return(div(class = "port-banner port-bad", div(
        class = "port-verdict",
        "Localization failed"
      ), div(class = "port-gloss", f$error)))
    }
    NULL
  })

  output$fault_print <- renderPrint({
    req(input$localize > 0L)
    f <- fault()
    req(is.null(f$error))
    print(f)
  })

  interval_cells <- reactive({
    R <- R_mat()
    p <- nrow(R)
    na_cells <- which(is.na(R) & upper.tri(R), arr.ind = TRUE)
    if (nrow(na_cells) > 1L) {
      return(list(error = sprintf(
        paste0(
          "%d cells were not reported. An interval can only be computed when ",
          "every other cell is present, so this is available for matrices with ",
          "at most one unreported cell."
        ),
        nrow(na_cells)
      )))
    }
    cells <- if (nrow(na_cells) == 1L) {
      matrix(as.integer(na_cells[1L, ]), ncol = 2)
    } else {
      t(utils::combn(p, 2L))
    }
    list(cells = cells)
  })

  intervals <- reactive({
    R <- R_mat()
    ic <- interval_cells()
    validate(need(is.null(ic$error), ic$error))
    hold <- if (isTRUE(input$hold_box)) "box" else "points"
    run <- function() {
      implied_interval(
        R,
        cells = ic$cells,
        decimals = decimals(),
        hold = hold
      )
    }
    if (identical(hold, "box")) {
      withProgress(message = "Computing box-aware intervals", value = 0.5, run())
    } else {
      run()
    }
  })

  # The plot is capped so that a large matrix does not produce a several
  # thousand pixel tall image; the table below it always carries every cell.
  interval_plot_max <- 40L

  interval_plot_data <- reactive({
    iv <- intervals()
    if (nrow(iv) <= interval_plot_max) {
      return(list(iv = iv, dropped = 0L))
    }
    edit <- abs(ifelse(is.na(iv$required_edit), 0, iv$required_edit))
    ord <- order(iv$status == "feasible", -edit)
    list(
      iv = iv[ord[seq_len(interval_plot_max)], , drop = FALSE],
      dropped = nrow(iv) - interval_plot_max
    )
  })

  output$interval_plot_note <- renderUI({
    d <- interval_plot_data()$dropped
    if (d == 0L) {
      return(NULL)
    }
    div(
      class = "port-note",
      sprintf(
        paste0(
          "Plotting the %d cells needing the largest edit; the remaining %d ",
          "are in the table below."
        ),
        interval_plot_max,
        d
      )
    )
  })

  output$interval_plot <- renderPlot(
    {
      plot_intervals(
        interval_plot_data()$iv,
        delta = check()$delta,
        labels = colnames(R_mat())
      )
    },
    height = function() {
      d <- tryCatch(interval_plot_data(), error = function(e) NULL)
      n <- if (is.null(d)) 1L else nrow(d$iv)
      max(200, 80 + 22 * n)
    }
  )

  output$interval_table <- renderTable(
    {
      iv <- intervals()
      nm <- colnames(R_mat())
      data.frame(
        cell = sprintf("%s – %s", nm[iv$i], nm[iv$j]),
        reported = fmt(iv$reported),
        lower = fmt(iv$lo),
        upper = fmt(iv$hi),
        `required edit` = fmt(iv$required_edit),
        status = iv$status,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    },
    striped = TRUE,
    rownames = FALSE
  )

  output$triples_table <- renderTable(
    {
      R <- R_mat()
      tr <- inconsistent_triples(R, decimals = decimals())
      nm <- colnames(R)
      if (nrow(tr) == 0L) {
        return(data.frame(
          result = "No inconsistent triple: no 3x3 submatrix is refuted on its own.",
          check.names = FALSE
        ))
      }
      data.frame(
        triple = sprintf("%s, %s, %s", nm[tr$i], nm[tr$j], nm[tr$k]),
        `box-max determinant` = fmt(tr$box_max_det, digits = 4),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    },
    striped = TRUE,
    rownames = FALSE
  )

  # Offered only once there is something to download.
  output$download_ui <- renderUI({
    ok <- !is.null(input$file) && is.null(parsed()$error)
    if (!ok) {
      return(div(class = "port-note", "Results download appears once a matrix is loaded."))
    }
    downloadButton("download_results", "Download results (.xlsx)")
  })

  output$download_example <- downloadHandler(
    filename = function() "example_correlation_matrix.csv",
    content = function(file) file.copy(example_path, file, overwrite = TRUE)
  )

  output$download_results <- downloadHandler(
    filename = function() {
      paste0("port_results_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      R <- R_mat()
      chk <- check()
      nm <- colnames(R)

      summary_df <- data.frame(
        field = c(
          "verdict",
          "certified",
          "tier",
          "variables",
          "decimals",
          "rounding rule",
          "delta",
          "witness margin B",
          "rigorous bound B_upper",
          "excusable delta",
          "port version",
          "citation"
        ),
        value = c(
          chk$verdict,
          as.character(isTRUE(chk$certified)),
          chk$tier,
          as.character(chk$p),
          as.character(decimals()),
          chk$rounding,
          fmt(chk$delta, 6),
          fmt(chk$margin, 6),
          fmt(chk$b_upper, 6),
          fmt(excusable(), 6),
          as.character(utils::packageVersion("port")),
          port_citation
        ),
        stringsAsFactors = FALSE
      )

      reported <- as.data.frame(R)
      reported <- cbind(variable = nm, reported)

      iv <- tryCatch(intervals(), error = function(e) NULL)
      intervals_df <- if (is.null(iv)) {
        data.frame(note = "not available for this matrix")
      } else {
        data.frame(
          i = nm[iv$i],
          j = nm[iv$j],
          reported = iv$reported,
          lower = iv$lo,
          upper = iv$hi,
          required_edit = iv$required_edit,
          hold = iv$hold,
          status = iv$status,
          stringsAsFactors = FALSE
        )
      }

      tr <- inconsistent_triples(R, decimals = decimals())
      triples_df <- if (nrow(tr) == 0L) {
        data.frame(note = "no inconsistent triples")
      } else {
        data.frame(
          i = nm[tr$i],
          j = nm[tr$j],
          k = nm[tr$k],
          box_max_det = tr$box_max_det,
          stringsAsFactors = FALSE
        )
      }

      sheets <- list(
        summary = summary_df,
        reported_matrix = reported,
        cell_intervals = intervals_df,
        inconsistent_triples = triples_df
      )

      if (isTRUE(input$localize > 0L)) {
        f <- tryCatch(fault(), error = function(e) NULL)
        if (!is.null(f) && is.null(f$error)) {
          sheets$localization <- data.frame(
            field = c(
              "localization verdict",
              "severity (largest single edit)",
              "severity (Frobenius)",
              "best achievable lambda_min",
              "excusable delta",
              "notes"
            ),
            value = c(
              f$localization_verdict,
              fmt(f$severity$severity_max, 6),
              fmt(f$severity$severity_frob, 6),
              fmt(f$severity$best_lambda_min, 6),
              fmt(f$severity$excusable_delta, 6),
              paste(f$notes, collapse = " | ")
            ),
            stringsAsFactors = FALSE
          )
        }
      }

      writexl::write_xlsx(sheets, path = file)
    }
  )

  # -- panel 2: three correlations -------------------------------------------

  trio_R <- reactive({
    v <- c(input$rAB, input$rAC, input$rBC)
    validate(need(
      all(!is.na(v)) && all(abs(v) <= 1),
      "Each correlation must be a number in [-1, 1]."
    ))
    R <- matrix(
      c(
        1, v[1], v[2],
        v[1], 1, v[3],
        v[2], v[3], 1
      ),
      3,
      3
    )
    dimnames(R) <- list(c("A", "B", "C"), c("A", "B", "C"))
    R
  })

  trio_decimals <- reactive({
    d <- input$trio_decimals
    validate(need(
      is.numeric(d) && !is.na(d) && d >= 0 && d <= 8,
      "Decimal places must be between 0 and 8."
    ))
    as.integer(d)
  })

  trio_check <- reactive({
    check_corr_psd(trio_R(), decimals = trio_decimals())
  })

  trio_intervals <- reactive({
    implied_interval(
      trio_R(),
      cells = rbind(c(1L, 2L), c(1L, 3L), c(2L, 3L)),
      decimals = trio_decimals(),
      hold = "box"
    )
  })

  output$trio_verdict_ui <- renderUI({
    chk <- trio_check()
    tagList(
      verdict_banner(chk),
      stat_row(chk, excusable_delta(trio_R()))
    )
  })

  output$trio_print <- renderPrint({
    print(trio_check())
  })

  output$trio_plot <- renderPlot({
    plot_intervals(
      trio_intervals(),
      delta = trio_check()$delta,
      labels = colnames(trio_R())
    )
  })

  output$trio_table <- renderTable(
    {
      iv <- trio_intervals()
      nm <- colnames(trio_R())
      data.frame(
        cell = sprintf("r(%s, %s)", nm[iv$i], nm[iv$j]),
        reported = fmt(iv$reported),
        lower = fmt(iv$lo),
        upper = fmt(iv$hi),
        `required edit` = fmt(iv$required_edit),
        status = iv$status,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    },
    striped = TRUE,
    rownames = FALSE
  )
}

shinyApp(ui = ui, server = server)
