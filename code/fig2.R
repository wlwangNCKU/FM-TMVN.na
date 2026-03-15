library(ggplot2)
################################################################################
#
#   Filename: fig2.R
#   Purpose:  Generate RMSE plots of estimated mixture model parameters across
#             varying sample sizes (n = 300, 600, 1200, 1800) and missing rates
#             (0%, 10%, 20%, 30%) from simulation results.
#   Input data files: results/simulation/para/*.csv (estimated parameters)
#                     results/simulation/sd/*.csv   (estimated standard errors)
#   Output data files: results/fig2.eps 
#   R Version: R-4.4.1
#
################################################################################
diff_result = matrix(NA, 16, 11)
for (i in 1:16){
  position = switch(i,"(300_0)","(600_0)","(1200_0)","(1800_0)","(300_0.1)","(600_0.1)","(1200_0.1)","(1800_0.1)","(300_0.2)","(600_0.2)","(1200_0.2)","(1800_0.2)","(300_0.3)","(600_0.3)","(1200_0.3)","(1800_0.3)")
  paraest = as.matrix(read.csv(paste0(SPATH,"/results/simulation/para/",position,".csv"),header = F))
  true.para=c(0.4,1.2,1.8,1.4,0.9,0.8,3.8,3.2,0.7,-0.6,1.3)
  diff = t(t(paraest)-true.para)
  diff_average = t(colMeans(diff^2))
  diff_result[i,] = diff_average
}
data = data.frame(
  parameter = rep(c("pi_1","mu_11","mu_12","Sigma_111","Sigma_121","Sigma_122","mu_21","mu_22","Sigma_211","Sigma_221",
                    "Sigma_222"), each = 16),
  dropout = rep(c("0%","0%","0%","0%","10%","10%","10%","10%","20%","20%","20%","20%","30%","30%","30%","30%"), times = 11),
  N = rep(c(300,600,1200,1800), times = 44),
  MSE = sqrt(as.vector(diff_result))
)
data$dropout = factor(data$dropout, levels = c("0%", "10%", "20%", "30%"))
data$parameter = factor(data$parameter, 
                        levels = c("pi_1", "mu_11", "mu_12", "Sigma_111", "Sigma_121", "Sigma_122",
                                   "mu_21", "mu_22", "Sigma_211", "Sigma_221", "Sigma_222"),
                        labels = c(expression(pi[1]),
                                   expression(bold(mu)[11]), expression(bold(mu)[12]), 
                                   expression(bold(Sigma)[111]), expression(bold(Sigma)[121]), expression(bold(Sigma)[122]),
                                   expression(bold(mu)[21]), expression(bold(mu)[22]), 
                                   expression(bold(Sigma)[211]), expression(bold(Sigma)[221]), expression(bold(Sigma)[222])))
data$N = factor(data$N)

postscript(paste0(SPATH, "results/fig2.eps", width = 12, height = 12))
ggplot(data, aes(x = N, y = MSE, color = dropout, group = dropout, linetype = dropout, shape = dropout)) +
  geom_line() + 
  geom_point(size = 2.5) + 
  facet_wrap(~ parameter, scales = "free", ncol = 4, labeller = label_parsed) + 
  labs(x = "Sample size (n)", y = "Root Mean Squared Error (RMSE)", color = "Missing Rate") +
  theme(
    strip.text = element_text(size = 20), 
    legend.position = c(0.95, 0.055), 
    legend.justification = c(1, 0), 
    legend.direction = "vertical",
  ) +
  guides(
    color = guide_legend(title = "Missing Rate"),
    linetype = guide_legend(title = "Missing Rate"),
    shape = guide_legend(title = "Missing Rate")
  )
dev.off()
