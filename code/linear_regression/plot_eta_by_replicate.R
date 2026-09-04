## code/linear_regression/plot_eta_by_replicate.R
##
## Plots the density of eta_hat (one point estimate per replicate, from
## sequential_K_for_eta) for each congruence level, using the checkpoint
## produced by eta_estimation_by_replicate.R. Can be re-run at any point
## while that script is still running (or resuming) to see progress so far
## -- it just reads whatever's in the checkpoint at the time.

library(dplyr)
library(ggplot2)

RESULTS_PATH <- "samples/linear_regression/eta_by_replicate.RData"
FIG_PATH     <- "figures/lm_eta_density.png"

stopifnot(file.exists(RESULTS_PATH))
results_df <- readRDS(RESULTS_PATH)

ok <- results_df %>% filter(is.na(error))
n_failed <- nrow(results_df) - nrow(ok)
if (n_failed > 0) {
  cat(sprintf("Note: %d / %d combos failed and are excluded from the plot (see results_df$error).\n",
              n_failed, nrow(results_df)))
}

CONGRUENCE_LEVELS <- c("High congruence", "Small congruence", "No congruence")
ok$congruence <- factor(ok$congruence, levels = CONGRUENCE_LEVELS)

cat("Replicates per congruence level going into the plot:\n")
print(table(ok$congruence))

medians <- ok %>%
  group_by(congruence) %>%
  summarise(median_eta = median(eta_hat), .groups = "drop")

# named so fill and vline color always match by congruence label, regardless
# of factor/group ordering -- avoids the manual eyeballed-order matching in
# the original sample.r plot (3 separate geom_vline calls with hardcoded
# colors, error-prone if the group order ever changed)
palette <- c(
  "High congruence"  = "#1b7837",
  "Small congruence" = "#e08214",
  "No congruence"    = "#762a83"
)

labels <- c(
  "High congruence"  = "No discrepancy",
  "Small congruence" = "Small discrepancy",
  "No congruence"    = "Large discrepancy"
)

p <- ggplot(ok, aes(x = eta_hat, fill = congruence)) +
  geom_density(alpha = 0.5) +
  geom_vline(data = medians, aes(xintercept = median_eta, color = congruence),
             linetype = "dashed", linewidth = 0.6, show.legend = FALSE) +
  scale_fill_manual(values = palette, labels = labels) +
  scale_color_manual(values = palette, labels = labels) +
  labs(x = expression(hat(eta)), y = "", fill = "Discrepancy level") +
  theme_bw()

print(p)

if (!dir.exists(dirname(FIG_PATH))) dir.create(dirname(FIG_PATH), recursive = TRUE)
ggsave(FIG_PATH, p, width = 6, height = 4)
cat(sprintf("Saved plot to %s\n", FIG_PATH))
