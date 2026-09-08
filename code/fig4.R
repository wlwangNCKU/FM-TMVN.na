library(ggplot2)
################################################################################
#
#   Filename: fig4.R
#   Purpose: Plot and compare BIC values of FM‐TMVN and GMM models across 
#            different numbers of clusters (g = 1 to 6).
#   Input data files: data/stone_result.RData
#   Output data files: results/fig4.eps 
#   R Version: R-4.4.1
#   Required R packages: ggplot2
#
################################################################################

load(paste(SPATH, "/data/stone_result.RData", sep=''))

g = 1:6
bic_fm_tmvn = round(sapply(g, function(g) total[[ paste0(g, "_1") ]]$BIC),3)
bic_gmm = round(sapply(g, function(g) total[[ paste0(g, "_2") ]]$BIC),3)
df = data.frame(
  g = rep(g, 2),
  BIC = c(bic_gmm, bic_fm_tmvn),
  Model = factor(rep(c("FM-TMVN", "GMM"), each = length(g)))
)
highlight_points = data.frame(
  g = c(2, 2),
  BIC = c(955.0744, 1029.274)
)

postscript(paste0(SPATH, "results/fig4.eps", width = 5, height = 5))
ggplot(df, aes(x = g, y = BIC, color = Model, shape = Model, linetype = Model)) +
  geom_line(size = 1.1) +
  geom_point(size = 3) +
  geom_point(data = highlight_points, aes(x = g, y = BIC),
             inherit.aes = FALSE, color = "red", shape = 1, size = 5, stroke = 1.5) +
  scale_color_manual(values = c("FM-TMVN" = "blue", "GMM" = "#00BFBF")) +
  scale_shape_manual(values = c("FM-TMVN" = 16, "GMM" = 15)) +
  scale_linetype_manual(values = c("FM-TMVN" = "solid", "GMM" = "solid")) +
  scale_x_continuous(breaks = 1:6) +
  labs(x = "Number of Clusters (g)", y = "BIC", color = NULL, shape = NULL, linetype = NULL) +
  theme_classic(base_size = 14) +
  annotate("text", x = 2.4, y = 955.0744 - 8, label = "955.07", color = "red", size = 3.3) +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7)) +
  theme(
    legend.position = c(0.182, 0.89),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.4)
  )
dev.off()