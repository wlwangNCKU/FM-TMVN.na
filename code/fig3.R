library(ggplot2)
library(GGally)
library(rlang)
################################################################################
#
#   Filename: fig3.R
#   Purpose:  Generate pairwise scatter–histogram plots with marginal histograms  
#             and correlation annotations for the four fluorescent markers in  
#             the HSCT dataset.
#   Input data files: data/source/hsct.csv
#   Output data files: results/fig3.eps 
#   R Version: R-4.4.1
#   Required R packages: ggplot2, GGally, rlang
#
################################################################################

data = read.csv(paste(SPATH, "/data/source/hsct.csv", sep=''), header = T)[,-1]
postscript(paste0(SPATH, "results/fig3.eps", width = 10, height = 10))
row0 = apply(data, 1, function(row) any(row == 0))
data.m = data[row0,]
data.o = data[!row0,]
set.seed(123)  
data_half = data.frame()
cls_levels = unique(data.m$cls)
for (g in cls_levels) {
  subdata = data.m[data.m$cls == g, ] 
  n = nrow(subdata)
  idx = sample(seq_len(n), size = floor(n/80))
  subdata.o = data.o[data.o$cls == g, ]
  n.o = nrow(subdata.o)
  idx.o = sample(seq_len(n.o), size = floor(n.o/1.5))
  data_half = rbind(data_half, subdata[idx, ],subdata.o[idx.o, ])
}
data = data_half
data$cls = as.factor(data$cls)
vars = sapply(data, is.numeric)
data_plot = data[, vars]
cls = data$cls

# Diagonal: histograms
yy = function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_histogram(position = 'identity', alpha = 0.7,
                   aes(y = ..density.., fill = cls), bins = 30)+
    theme_void()
}

# Lower triangle: scatter plots
zz = function(data, mapping, ...) {
  x_var = as_name(mapping$x)
  y_var = as_name(mapping$y)
  p = ggplot(data, mapping = mapping) +
    geom_point(aes(colour = cls), size = 0.3)
  # Find observations where x == 0 or y == 0
  zero_points = data[data[[x_var]] == 0 | data[[y_var]] == 0, ]
  # If such points exist, plot them as black 'x'
  if (nrow(zero_points) > 0) {
    p = p + geom_point(data = zero_points,
                       aes(x = .data[[x_var]], y = .data[[y_var]]),
                       shape = 4, size = 0.7, color = "black",
                       inherit.aes = FALSE)}
  return(p)
}
ggpairs(data_plot,
        mapping = aes(colour = cls),
        upper = list(continuous = wrap("cor", size = 4,stars = F)),
        lower = list(continuous = zz),
        diag  = list(continuous = yy),
        axisLabels = 'none',
        legend = c(2, 2))
dev.off()
