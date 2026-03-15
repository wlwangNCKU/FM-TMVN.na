rm(list = ls())
SPATH = paste(getwd(),"/Data_and_code_FM_TMVN.na",sep="")

# Re-produce Figure 1
source(paste(SPATH, '/code/fig1.R', sep=''))

# Re-produce Figure 2
source(paste(SPATH, '/code/fig2.R', sep=''))

# Re-produce Figure 3
source(paste(SPATH, '/code/fig3.R', sep=''))

# Re-produce Figure 4
source(paste(SPATH, '/code/fig4.R', sep=''))

# Re-produce Figure 5
source(paste(SPATH, '/code/fig5.R', sep=''))

# Re-produce Figure F1
source(paste(SPATH, '/code/figF1.R', sep=''))

# Re-produce Table 1
source(paste(SPATH, '/code/Table1.R', sep=''))

# Re-produce Table 2
source(paste(SPATH, '/code/Table2.R', sep=''))

# Re-produce Table 3
source(paste(SPATH, '/code/Table3.R', sep=''))

# Re-produce Table F1
source(paste(SPATH, '/code/TabF1.R', sep=''))

# Re-produce Simulation
source(paste(SPATH, '/code/simulation.R', sep=''))

# Re-produce hsct_result.RData
source(paste(SPATH, '/code/fit_hsct.R', sep=''))

# Re-produce stone_result.RData
source(paste(SPATH, '/code/fit_stone.R', sep=''))
