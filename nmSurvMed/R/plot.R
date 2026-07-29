# Plot complete mediation results -----------------------------------------

.plot_effect_panel <- function(table, effect, title, color, ci_color,
                               show_legend, conf_level, ylim) {
  estimate <- table[[effect]]
  lower <- table[[paste0(effect, "_lower")]]
  upper <- table[[paste0(effect, "_upper")]]
  graphics::plot(
    table$time,
    estimate,
    type = "n",
    ylim = ylim,
    xlab = "Follow-up time",
    ylab = "Effect Estimate",
    main = title
  )
  graphics::polygon(
    c(table$time, rev(table$time)),
    c(lower, rev(upper)),
    col = ci_color,
    border = NA
  )
  graphics::lines(table$time, estimate, col = color, lwd = 2.5)
  graphics::abline(h = 0, lty = 2, lwd = 1.5, col = "black")

  if (show_legend) {
    graphics::legend(
      "bottomleft",
      legend = c(
        "Non-Markov based",
        sprintf("%.0f%% pointwise CI", 100 * conf_level)
      ),
      col = c(color, ci_color),
      lty = c(1, NA),
      lwd = c(3, NA),
      pch = c(NA, 15),
      pt.cex = 2,
      bty = "n",
      cex = 1.05
    )
  }
}

#' Plot survival and cumulative-hazard mediation effects
#'
#' Produces a paper-style 2 by 2 figure containing survival direct and
#' indirect effects in the first row and cumulative-hazard direct and indirect
#' effects in the second row. Pointwise bootstrap confidence intervals are
#' shown as shaded regions.
#'
#' @param x An `nmSurvMed` fit created with `scale = "both"`.
#' @param file Optional output filename ending in `.png` or `.pdf`. When
#'   `NULL`, the figure is drawn on the active graphics device. The parent
#'   directory must already exist; this function never creates directories.
#' @param width,height Output dimensions. For PNG these are pixels; for PDF
#'   these are inches.
#' @param res PNG resolution in dots per inch.
#' @param color Color of the estimated effect curves.
#' @param ci_alpha Opacity of the confidence interval shading.
#' @param main_prefix Optional text appended to each panel title, such as a
#'   data-set or subgroup name.
#' @param show_legend Logical; display a legend in each panel. The default is
#'   `FALSE`.
#'
#' @return The input fit, invisibly.
#' @export
#'
#' @examples
#' data <- nmSurvMed:::nm_example_data(100)
#' time <- seq(0, ceiling(max(data$T2) * 10) / 10, by = 0.1)
#' fit <- nmSurvMed(
#'   data,
#'   T1 = "T1_hat", D1 = "D1",
#'   T2 = "T2_hat", D2 = "D2", Z = "Z",
#'   time = time, bandwidth = c(0.2, 0.3),
#'   folds = 2, bootstrap = 3, verbose = FALSE
#' )
#' plot_nmSurvMed(fit)
plot_nmSurvMed <- function(x, file = NULL, width = 4000, height = 3000,
                           res = 300, color = "blue", ci_alpha = 0.2,
                           main_prefix = NULL, show_legend = FALSE) {
  if (!inherits(x, "nmSurvMed")) {
    stop("`x` must be an object returned by `nmSurvMed()`.", call. = FALSE)
  }
  if (!identical(x$scale, "both") || !is.list(x$effects) ||
      is.null(x$effects$survival) ||
      is.null(x$effects$cumulative_hazard)) {
    stop(
      "The fit must be created with `scale = \"both\"`.",
      call. = FALSE
    )
  }
  if (length(ci_alpha) != 1L || !is.finite(ci_alpha) ||
      ci_alpha < 0 || ci_alpha > 1) {
    stop("`ci_alpha` must be between 0 and 1.", call. = FALSE)
  }

  opened_device <- FALSE
  if (!is.null(file)) {
    if (length(file) != 1L || !is.character(file)) {
      stop("`file` must be NULL or one filename.", call. = FALSE)
    }
    parent <- dirname(file)
    if (!dir.exists(parent)) {
      stop(
        "The output directory does not exist: ", parent,
        ". This function does not create directories.",
        call. = FALSE
      )
    }
    extension <- tolower(tools::file_ext(file))
    if (extension == "png") {
      grDevices::png(file, width = width, height = height, res = res)
    } else if (extension == "pdf") {
      grDevices::pdf(file, width = width, height = height)
    } else {
      stop("`file` must end in `.png` or `.pdf`.", call. = FALSE)
    }
    opened_device <- TRUE
  }

  old_parameters <- graphics::par(
    mfrow = c(2, 2),
    mar = c(5, 5, 4, 2),
    cex.lab = 1.35,
    cex.axis = 1.2,
    cex.main = 1.3
  )
  on.exit({
    graphics::par(old_parameters)
    if (opened_device) grDevices::dev.off()
  }, add = TRUE)

  ci_color <- grDevices::adjustcolor(color, alpha.f = ci_alpha)
  suffix <- if (is.null(main_prefix) || !nzchar(main_prefix)) {
    ""
  } else {
    paste0(" - ", main_prefix)
  }
  conf_level <- x$bootstrap$conf_level
  shared_limits <- function(table) {
    limits <- range(
      c(
        table$direct_lower, table$direct_upper,
        table$indirect_lower, table$indirect_upper,
        0
      ),
      finite = TRUE
    )
    if (diff(limits) == 0) {
      limits <- limits +
        c(-1, 1) * max(0.01, abs(limits[1L]) * 0.05)
    }
    padding <- diff(limits) * 0.03
    limits + c(-padding, padding)
  }
  survival_ylim <- shared_limits(x$effects$survival)
  cumulative_hazard_ylim <- shared_limits(
    x$effects$cumulative_hazard
  )

  .plot_effect_panel(
    x$effects$survival, "direct",
    paste0("Survival Direct Effect", suffix),
    color, ci_color, show_legend, conf_level, survival_ylim
  )
  .plot_effect_panel(
    x$effects$survival, "indirect",
    paste0("Survival Indirect Effect", suffix),
    color, ci_color, show_legend, conf_level, survival_ylim
  )
  .plot_effect_panel(
    x$effects$cumulative_hazard, "direct",
    paste0("Cumulative-Hazard Direct Effect", suffix),
    color, ci_color, show_legend, conf_level, cumulative_hazard_ylim
  )
  .plot_effect_panel(
    x$effects$cumulative_hazard, "indirect",
    paste0("Cumulative-Hazard Indirect Effect", suffix),
    color, ci_color, show_legend, conf_level, cumulative_hazard_ylim
  )

  invisible(x)
}

.plot_sensitivity_panel <- function(table, effect, title, colors, ylim,
                                    ci_alpha) {
  settings <- unique(table$setting)
  graphics::plot(
    range(table$time), ylim,
    type = "n", xlab = "Follow-up time", ylab = "Effect Estimate",
    main = title
  )
  for (index in seq_along(settings)) {
    current <- table[table$setting == settings[index], , drop = FALSE]
    lower <- current[[paste0(effect, "_lower")]]
    upper <- current[[paste0(effect, "_upper")]]
    graphics::polygon(
      c(current$time, rev(current$time)),
      c(lower, rev(upper)),
      col = grDevices::adjustcolor(colors[index], alpha.f = ci_alpha),
      border = NA
    )
  }
  for (index in seq_along(settings)) {
    current <- table[table$setting == settings[index], , drop = FALSE]
    graphics::lines(
      current$time, current[[effect]], col = colors[index], lwd = 2.5
    )
    graphics::lines(
      current$time, current[[paste0(effect, "_lower")]],
      col = colors[index], lty = 3, lwd = 1
    )
    graphics::lines(
      current$time, current[[paste0(effect, "_upper")]],
      col = colors[index], lty = 3, lwd = 1
    )
  }
  graphics::abline(h = 0, lty = 2, lwd = 1.5, col = "gray35")
}

#' Plot bandwidth sensitivity results
#'
#' Draws survival effects in the first row and cumulative-hazard effects in
#' the second row. Direct and indirect panels in the same row share a y-axis
#' scale. Each bandwidth setting is shown with its pointwise bootstrap
#' confidence interval.
#'
#' @param x An object returned by [bandwidth_sensitivity()] with
#'   `scale = "both"`.
#' @inheritParams plot_nmSurvMed
#' @param colors Colors used for the bandwidth settings.
#'
#' @return The input sensitivity object, invisibly.
#' @export
plot_bandwidth_sensitivity <- function(
    x, file = NULL, width = 4000, height = 3000, res = 300,
    colors = NULL, ci_alpha = 0.12,
    main_prefix = NULL, show_legend = FALSE) {
  if (!inherits(x, "nm_bandwidth_sensitivity")) {
    stop(
      "`x` must be an object returned by `bandwidth_sensitivity()`.",
      call. = FALSE
    )
  }
  if (!identical(x$scale, "both")) {
    stop("The sensitivity fit must use `scale = \"both\"`.", call. = FALSE)
  }
  settings <- unique(x$effects$survival$setting)
  if (is.null(colors)) {
    colors <- grDevices::hcl.colors(length(settings), palette = "Dark 3")
  }
  if (length(colors) < length(settings)) {
    stop("Provide at least one color per bandwidth setting.", call. = FALSE)
  }
  colors <- colors[seq_along(settings)]
  if (length(ci_alpha) != 1L || !is.finite(ci_alpha) ||
      ci_alpha < 0 || ci_alpha > 1) {
    stop("`ci_alpha` must be between 0 and 1.", call. = FALSE)
  }

  opened_device <- FALSE
  if (!is.null(file)) {
    if (length(file) != 1L || !is.character(file)) {
      stop("`file` must be NULL or one filename.", call. = FALSE)
    }
    if (!dir.exists(dirname(file))) {
      stop(
        "The output directory does not exist: ", dirname(file),
        ". This function does not create directories.",
        call. = FALSE
      )
    }
    extension <- tolower(tools::file_ext(file))
    if (extension == "png") {
      grDevices::png(file, width = width, height = height, res = res)
    } else if (extension == "pdf") {
      grDevices::pdf(file, width = width, height = height)
    } else {
      stop("`file` must end in `.png` or `.pdf`.", call. = FALSE)
    }
    opened_device <- TRUE
  }
  old_parameters <- graphics::par(
    mfrow = c(2, 2), mar = c(5, 5, 4, 2),
    cex.lab = 1.35, cex.axis = 1.2, cex.main = 1.3
  )
  on.exit({
    graphics::par(old_parameters)
    if (opened_device) grDevices::dev.off()
  }, add = TRUE)

  shared_limits <- function(table) {
    limits <- range(c(
      table$direct_lower, table$direct_upper,
      table$indirect_lower, table$indirect_upper, 0
    ), finite = TRUE)
    if (diff(limits) == 0) {
      limits <- limits + c(-1, 1) * max(0.01, abs(limits[1L]) * 0.05)
    }
    limits + c(-1, 1) * diff(limits) * 0.03
  }
  suffix <- if (is.null(main_prefix) || !nzchar(main_prefix)) {
    ""
  } else {
    paste0(" - ", main_prefix)
  }
  survival_ylim <- shared_limits(x$effects$survival)
  hazard_ylim <- shared_limits(x$effects$cumulative_hazard)

  .plot_sensitivity_panel(
    x$effects$survival, "direct",
    paste0("Survival Direct Effect", suffix),
    colors, survival_ylim, ci_alpha
  )
  .plot_sensitivity_panel(
    x$effects$survival, "indirect",
    paste0("Survival Indirect Effect", suffix),
    colors, survival_ylim, ci_alpha
  )
  .plot_sensitivity_panel(
    x$effects$cumulative_hazard, "direct",
    paste0("Cumulative-Hazard Direct Effect", suffix),
    colors, hazard_ylim, ci_alpha
  )
  .plot_sensitivity_panel(
    x$effects$cumulative_hazard, "indirect",
    paste0("Cumulative-Hazard Indirect Effect", suffix),
    colors, hazard_ylim, ci_alpha
  )
  if (show_legend) {
    graphics::legend(
      "bottomleft", legend = settings, col = colors,
      lty = 1, lwd = 2.5, bty = "n",
      cex = if (length(settings) > 6L) 0.7 else 0.9,
      ncol = if (length(settings) > 6L) 2L else 1L
    )
  }
  invisible(x)
}
