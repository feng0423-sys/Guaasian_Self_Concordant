#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

root_dir <- getwd()
out_dir <- file.path(root_dir, "html_pdf")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

per_path_file <- file.path(root_dir, "pn_pg_time_per_path_case_norm_mode_tol_1e-04_max_iter_10000.csv")
total_file <- file.path(root_dir, "pn_pg_time_total_case_norm_mode_tol_1e-04_max_iter_10000.csv")
comparison_file <- file.path(root_dir, "pn_pg_target_comparison_case_norm_mode_tol_1e-04_max_iter_10000.csv")

for (f in c(per_path_file, total_file, comparison_file)) {
  if (!file.exists(f)) stop(sprintf("Missing required file: %s", f))
}

per_path <- read.csv(per_path_file, check.names = FALSE)
total_df <- read.csv(total_file, check.names = FALSE)
cmp_df <- read.csv(comparison_file, check.names = FALSE)
METHOD_COLORS <- c(PN = "#1f77b4", PG = "#2ca02c")

# Build loss and iteration long tables.
loss_df <- rbind(
  data.frame(path_index = cmp_df$path_index, lambda_gamma = cmp_df$lambda_gamma, method = "PN", final_loss = cmp_df$pn_final_loss, stringsAsFactors = FALSE),
  data.frame(path_index = cmp_df$path_index, lambda_gamma = cmp_df$lambda_gamma, method = "PG", final_loss = cmp_df$pg_final_loss_at_stop, stringsAsFactors = FALSE)
)
loss_df$log_lambda <- log(loss_df$lambda_gamma)
loss_df$method <- factor(loss_df$method, levels = c("PN", "PG"))

iter_df <- rbind(
  data.frame(path_index = cmp_df$path_index, lambda_gamma = cmp_df$lambda_gamma, method = "PN", outer_iters = cmp_df$pn_outer_iters, stringsAsFactors = FALSE),
  data.frame(path_index = cmp_df$path_index, lambda_gamma = cmp_df$lambda_gamma, method = "PG", outer_iters = cmp_df$pg_outer_iters, stringsAsFactors = FALSE)
)
iter_df$log_lambda <- log(iter_df$lambda_gamma)
iter_df$method <- factor(iter_df$method, levels = c("PN", "PG"))

# Keep total table ordered PN then PG.
total_df$method <- factor(total_df$method, levels = c("PN", "PG"))
total_df <- total_df[order(total_df$method), , drop = FALSE]

# Plot 1: training loss vs log(lambda).
p_loss <- ggplot(loss_df, aes(x = log_lambda, y = final_loss, color = method, shape = method, group = method)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.0) +
  scale_color_manual(values = METHOD_COLORS) +
  labs(
    title = "Training Loss vs log(lambda)",
    subtitle = "case_norm_mode | PN-first + PG-target workflow",
    x = "log(lambda)",
    y = "Training loss",
    color = NULL,
    shape = NULL
  ) +
  theme_minimal(base_size = 12)

# Plot 2: outer iterations vs log(lambda).
p_iter <- ggplot(iter_df, aes(x = log_lambda, y = outer_iters, color = method, shape = method, group = method)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.0) +
  scale_color_manual(values = METHOD_COLORS) +
  labs(
    title = "Outer Iterations vs log(lambda)",
    subtitle = "case_norm_mode | PN-first + PG-target workflow",
    x = "log(lambda)",
    y = "Outer iterations",
    color = NULL,
    shape = NULL
  ) +
  theme_minimal(base_size = 12)

# Plot 3: total simulation time only.
p_total <- ggplot(total_df, aes(x = method, y = total_elapsed_sec, fill = method)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.2f", total_elapsed_sec)), vjust = -0.4, size = 4) +
  scale_fill_manual(values = METHOD_COLORS) +
  labs(
    title = "Total Simulation Time",
    subtitle = "Total elapsed seconds by method",
    x = NULL,
    y = "Total elapsed time (sec)",
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  guides(fill = "none")

# Save PNGs for HTML embedding.
png_loss <- file.path(out_dir, "simulation_training_loss_vs_loglambda.png")
png_iter <- file.path(out_dir, "simulation_outer_iterations_vs_loglambda.png")
png_total <- file.path(out_dir, "simulation_total_time_bar.png")

ggsave(png_loss, p_loss, width = 10, height = 5, dpi = 300)
ggsave(png_iter, p_iter, width = 10, height = 5, dpi = 300)
ggsave(png_total, p_total, width = 8, height = 5, dpi = 300)

# Build PDF report.
pdf_path <- file.path(out_dir, "simulation_time_comparison_case_norm_mode_tol_1e-04_max_iter_10000.pdf")
pdf(pdf_path, width = 11, height = 8.5, onefile = TRUE)
print(p_loss)
print(p_iter)
print(p_total)

plot.new()
title("Total Simulation Time Table", cex.main = 1.2)
text_total <- capture.output(print(total_df, row.names = FALSE))
text(0.01, 0.95, paste(text_total, collapse = "\n"), adj = c(0, 1), family = "mono", cex = 0.9)

dev.off()

# Build HTML report.
html_path <- file.path(out_dir, "simulation_time_comparison_case_norm_mode_tol_1e-04_max_iter_10000.html")

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

meta_line <- sprintf(
  "Settings: case_norm_mode, tol_outer=%s, tol_inner_base=%s, tol_inner_floor=%s, max_iter_pg_target=%s",
  unique(per_path$tol_outer),
  unique(per_path$tol_inner_base),
  unique(per_path$tol_inner_floor),
  unique(total_df$max_iter_pg_target)
)

html_lines <- c(
  "<!DOCTYPE html>",
  "<html><head><meta charset='utf-8'><title>Simulation Time Comparison</title>",
  "<style>",
  "body{font-family:Arial,Helvetica,sans-serif; margin:24px; line-height:1.4;}",
  "h1,h2,h3{margin:0 0 10px 0;}",
  ".meta{margin-bottom:18px;color:#333;}",
  ".plot{margin:18px 0 30px 0;}",
  "img{max-width:1100px;width:100%;border:1px solid #ddd;}",
  "table{border-collapse:collapse; margin:10px 0 24px 0; font-size:13px;}",
  "th,td{border:1px solid #bbb; padding:6px 8px; text-align:left;}",
  "th{background:#f5f5f5;}",
  "</style></head><body>",
  "<h1>Simulation Comparison Report</h1>",
  sprintf("<div class='meta'>Generated at: %s</div>", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  sprintf("<div class='meta'>%s</div>", escape_html(meta_line)),
  "<h2>Requested Plots</h2>",
  "<div class='plot'><h3>Training Loss vs log(lambda)</h3><img src='simulation_training_loss_vs_loglambda.png' alt='Training loss vs log lambda'></div>",
  "<div class='plot'><h3>Outer Iterations vs log(lambda)</h3><img src='simulation_outer_iterations_vs_loglambda.png' alt='Outer iterations vs log lambda'></div>",
  "<div class='plot'><h3>Total Simulation Time</h3><img src='simulation_total_time_bar.png' alt='Total simulation time'></div>",
  "<h2>Total Simulation Time Table</h2>",
  build_html_table(total_df),
  "</body></html>"
)

writeLines(html_lines, con = html_path, useBytes = TRUE)

cat("Generated report files:\n")
cat(sprintf("- %s\n", html_path))
cat(sprintf("- %s\n", pdf_path))
