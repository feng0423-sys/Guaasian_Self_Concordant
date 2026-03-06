#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

root_dir <- getwd()
out_dir <- file.path(root_dir, "html_pdf")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

summary_files <- list(
  list(id = "A", tol_outer = "1e-04", tol_inner = "1e-05", path = file.path(root_dir, "lambda_summary_case_norm_mode_tol_1e-04_inner_1e-05.csv")),
  list(id = "B", tol_outer = "1e-05", tol_inner = "1e-06", path = file.path(root_dir, "lambda_summary_case_norm_mode_tol_1e-05_inner_1e-06.csv"))
)

for (i in seq_along(summary_files)) {
  if (!file.exists(summary_files[[i]]$path)) {
    stop(sprintf("Missing required summary file: %s", summary_files[[i]]$path))
  }
}

plot_records <- list()
plot_index <- 0L
comp_time_rows <- list()
comp_time_index <- 0L

make_setting_plots <- function(df, setting_label, tol_outer, tol_inner) {
  df$log_lambda <- log(df$lambda_gamma)
  df$method <- factor(df$method, levels = c("Prox-Newton", "Prox-Gradient"))

  subtitle <- sprintf("case_norm_mode | tol_outer=%s | tol_inner_base=%s", tol_outer, tol_inner)

  p_loss <- ggplot(df, aes(x = log_lambda, y = final_loss, color = method, shape = method, group = method)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.0) +
    scale_shape_manual(values = c(`Prox-Newton` = 4, `Prox-Gradient` = 16)) +
    labs(
      title = sprintf("Training Loss vs log(lambda) [%s]", setting_label),
      subtitle = subtitle,
      x = "log(lambda)",
      y = "Training loss",
      color = NULL,
      shape = NULL
    ) +
    theme_minimal(base_size = 12)

  p_iter <- ggplot(df, aes(x = log_lambda, y = outer_iters, color = method, shape = method, group = method)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.0) +
    scale_shape_manual(values = c(`Prox-Newton` = 4, `Prox-Gradient` = 16)) +
    labs(
      title = sprintf("Iterations vs log(lambda) [%s]", setting_label),
      subtitle = subtitle,
      x = "log(lambda)",
      y = "Outer iterations",
      color = NULL,
      shape = NULL
    ) +
    theme_minimal(base_size = 12)

  list(loss = p_loss, iter = p_iter)
}

for (spec in summary_files) {
  df <- read.csv(spec$path, check.names = FALSE)
  plots <- make_setting_plots(
    df = df,
    setting_label = sprintf("Setting %s", spec$id),
    tol_outer = spec$tol_outer,
    tol_inner = spec$tol_inner
  )

  loss_png <- file.path(out_dir, sprintf("case_norm_mode_training_loss_setting_%s.png", spec$id))
  iter_png <- file.path(out_dir, sprintf("case_norm_mode_iterations_setting_%s.png", spec$id))

  ggsave(loss_png, plots$loss, width = 10, height = 5, dpi = 300)
  ggsave(iter_png, plots$iter, width = 10, height = 5, dpi = 300)

  plot_index <- plot_index + 1L
  plot_records[[plot_index]] <- list(type = "Training loss", setting = spec$id, png = loss_png, plot = plots$loss)
  plot_index <- plot_index + 1L
  plot_records[[plot_index]] <- list(type = "Iterations", setting = spec$id, png = iter_png, plot = plots$iter)

  total_by_method <- aggregate(elapsed_sec ~ method, data = df, sum)
  pn_time <- total_by_method$elapsed_sec[total_by_method$method == "Prox-Newton"]
  pg_time <- total_by_method$elapsed_sec[total_by_method$method == "Prox-Gradient"]
  if (length(pn_time) == 0) pn_time <- NA_real_
  if (length(pg_time) == 0) pg_time <- NA_real_

  comp_time_index <- comp_time_index + 1L
  comp_time_rows[[comp_time_index]] <- data.frame(
    setting_id = spec$id,
    tol_outer = spec$tol_outer,
    tol_inner_base = spec$tol_inner,
    tol_inner_floor = "1e-08",
    pn_total_elapsed_sec = pn_time[1],
    pg_total_elapsed_sec = pg_time[1],
    stringsAsFactors = FALSE
  )
}

comp_time_table <- do.call(rbind, comp_time_rows)

# Build PDF report.
pdf_path <- file.path(out_dir, "simulation_study_stop_modes_latest_tol_1e-04_case_norm_mode_report.pdf")
pdf(pdf_path, width = 11, height = 8.5, onefile = TRUE)
for (rec in plot_records) {
  print(rec$plot)
}

# Timing pages as text tables using base graphics.
plot.new()
title("Computation Time + Tolerance Settings", cex.main = 1.2)
text_lines <- capture.output(print(comp_time_table, row.names = FALSE))
text(0.01, 0.95, paste(text_lines, collapse = "\n"), adj = c(0, 1), family = "mono", cex = 0.85)

dev.off()

# Build HTML report with embedded PNG links and timing tables.
html_path <- file.path(out_dir, "simulation_study_stop_modes_latest_tol_1e-04_case_norm_mode_report.html")

escape_html <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

build_html_table <- function(df) {
  cols <- names(df)
  header <- paste0("<tr>", paste0("<th>", escape_html(cols), "</th>", collapse = ""), "</tr>")
  rows <- apply(df, 1, function(r) {
    paste0("<tr>", paste0("<td>", escape_html(as.character(r)), "</td>", collapse = ""), "</tr>")
  })
  paste0("<table>", header, paste0(rows, collapse = ""), "</table>")
}

html_lines <- c(
  "<!DOCTYPE html>",
  "<html><head><meta charset='utf-8'><title>Case Norm Mode Report</title>",
  "<style>",
  "body{font-family:Arial,Helvetica,sans-serif; margin:24px; line-height:1.4;}",
  "h1,h2,h3{margin:0 0 10px 0;}",
  ".meta{margin-bottom:20px;color:#333;}",
  ".plot{margin:20px 0 32px 0;}",
  "img{max-width:1100px;width:100%;border:1px solid #ddd;}",
  "table{border-collapse:collapse; margin:10px 0 24px 0; font-size:13px;}",
  "th,td{border:1px solid #bbb; padding:6px 8px; text-align:left;}",
  "th{background:#f5f5f5;}",
  "</style></head><body>",
  "<h1>Simulation Study Report (case_norm_mode)</h1>",
  sprintf("<div class='meta'>Generated at: %s</div>", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "<h2>Settings</h2>",
  "<ul>",
  "<li>Setting A: tol_outer=1e-04, tol_inner_base=1e-05, tol_inner_floor=1e-08</li>",
  "<li>Setting B: tol_outer=1e-05, tol_inner_base=1e-06, tol_inner_floor=1e-08</li>",
  "<li>Stopping mode: case_norm_mode (PN: local_norm, PG: l2_step_norm)</li>",
  "</ul>",
  "<h2>All Plots</h2>"
)

for (rec in plot_records) {
  rel <- basename(rec$png)
  html_lines <- c(
    html_lines,
    sprintf("<div class='plot'><h3>%s (Setting %s)</h3>", rec$type, rec$setting),
    sprintf("<img src='%s' alt='%s Setting %s'></div>", rel, rec$type, rec$setting)
  )
}

html_lines <- c(
  html_lines,
  "<h2>Computation Time + Tolerance Record</h2>",
  build_html_table(comp_time_table),
  "</body></html>"
)

writeLines(html_lines, con = html_path, useBytes = TRUE)

cat("Generated report files:\n")
cat(sprintf("- %s\n", html_path))
cat(sprintf("- %s\n", pdf_path))
