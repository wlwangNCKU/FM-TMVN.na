# FM-TMVN.na
Supplement: Data and Code for "Maximum likelihood estimation for mixtures of truncated multivariate normal distributions with missing information" by Wan-Lun Wang, Victor H. Lachos, Pei-Ting Ho and Tsung-I Lin*

# Author responsible for the code #
For questions, comments or remarks about the code please contact responsible authors, Dr. Tsung-I Lin (tilin@nchu.edu.tw) and Dr. Wan-Lun Wang (wangwl@gs.ncku.edu.tw)

# Configurations #
The code was written/evaluated in R with the following software 
versions: R version 4.5.0 (2025-04-11 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26100)

Matrix products: default
locale:
[1] LC_COLLATE=Chinese (Traditional)_Taiwan.utf8  
    LC_CTYPE=Chinese (Traditional)_Taiwan.utf8   
[3] LC_MONETARY=Chinese (Traditional)_Taiwan.utf8 
    LC_NUMERIC=C                                 
[5] LC_TIME=Chinese (Traditional)_Taiwan.utf8    

attached base packages:  
[1] grid      stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
 [1] relliptical_1.3.0 MASS_7.3-65       imager_1.0.3      magrittr_2.0.3    VIM_6.2.2         colorspace_2.1-1  mnormt_2.1.1     
 [8] mclust_6.1.1      MomTrunc_6.1      mvtnorm_1.3-3    

loaded via a namespace (and not attached):
 [1] tiff_0.1-12            xml2_1.3.8             class_7.3-23           robustbase_0.99-4-1    jpeg_0.1-11           
 [6] stringi_1.8.7          lattice_0.22-6         Matrix_1.7-3           e1071_1.7-16           nnet_7.3-20           
[11] vcd_1.4-13             hypergeo_1.2-14        deSolve_1.40           Formula_1.2-5          Ryacas0_0.4.4         
[16] laeken_0.5.3           FuzzyNumbers.Ext.2_3.2 purrr_1.0.4            tlrmvnmvt_1.1.2        Rdpack_2.6.4          
[21] abind_1.4-8            cli_3.6.5              rlang_1.1.6            contfrac_1.1-12        rbibutils_2.3         
[26] tools_4.5.0            ranger_0.17.0          boot_1.3-31            vctrs_0.6.5            png_0.1-8             
[31] bmp_0.3                settings_0.2.7         zoo_1.8-14             proxy_0.4-27           lifecycle_1.0.4       
[36] stringr_1.5.1          car_3.1-3              pkgconfig_2.0.3        data.table_1.17.0      glue_1.8.0            
[41] Rcpp_1.0.14            FuzzyNumbers_0.4-7     DEoptimR_1.1-3-1       lmtest_0.9-40          readbitmap_0.1.5      
[46] igraph_2.1.4           carData_3.0-5          elliptic_1.4-0         matrixcalc_1.0-6       compiler_4.5.0  


# Descriptions of the codes # 
Please extract the file "Data_and_code_FM_TMVN.na.zip" to the current working directory of R. You can use the getwd() function to determine the current directory path.

Before running the codes **fig1.R**, **fig2.R**, **fig3.R**, **fig4.R**, **fig5.R**, **figF1.R**, **fit_hsct.R**, **fit_stone.R**, **simulation.R**, **Table1.R**, **Table2.R**,
**Table3.R** and **TableF1.R**, one needs to install the following R packages:

    install.packages("mvtnorm")      # Version: 1.3-3 
    install.packages("MASS")         # Version: 7.3-65 
    install.packages("relliptical")  # Version: 1.3.0 
    install.packages("MomTrunc")     # Version: 6.1  
    install.packages("ggplot2")      # Version: 3.5.2
    install.packages("GGally")       # Version: 2.2.1
    install.packages("rlang")        # Version: 1.1.6
    install.packages("mclust")       # Version: 6.1.1
    install.packages("mnormt")       # Version: 2.1.1
    install.packages("VIM")          # Version: 6.2.2
    install.packages("imager")       # Version: 1.0.3
    install.packages("fBasics")      # Version: 4041.97

R codes for the implementation of our methodology are provided.

# Data and Code for function #
## Subfolder: ./function ##
`./function` contains the program (function) of

- (1) **FMTMVN.na.EM.R** for calculating the parameter estimates via the ECM algorithm under the FM-TMVN model with missing information;
- (2) **gener_na.R** for randomly generating missing values;
- (3) **GMIX.na.EM.R** for calculating the parameter estimates via the ECM algorithm under the GMM model with missing information;
- (4) **SE.FMTMVN.na.R** Computes standard errors of the FM-TMVN parameter estimates using the observed information matrix;
- (5) **Sim_FMTMVN.na.EM.R** Simulates datasets with missing values and performs ECM‐based estimation under the FM-TMVN model;
- (6) **Sim_GMIX.na.EM.R** Simulates datasets with missing values and performs ECM‐based estimation under the GMM model for comparative evaluation.

# Data and Code for code #
## Subfolder: ./code ##
`./code` contains

-(1) **fig1.R** main script to generate 3D perspective views and contour plots of finite mixtures of truncated bivariate normal (FM-TMVN) distributions, varying the number of components (g = 1, 2, 3) and truncation schemes (none, right, left, and double);
-(2) **simulation.R**  main script for re-generating simulated datasets and estimating model parameters under different sample sizes and missing rates; 
-(3) **fig2.R**  generates RMSE line‐plots of all estimated parameters across sample sizes (n = 300, 600, 1200, 1800) and missing‐data rates (0%, 10%, 20%, 30%);
-(4) **Table1.R** computes and exports a summary table of simulation results—reporting STD and IMSE for parameter estimates under missing‐data rates r = 0%, 10%, 20%, 30% at sample sizes n = 300, 600, 1200, 1800;
-(5) **fit_hsct.R** main script for model fitting (GMM, FM-TMVN) to the hsct dataset; saved results are stored in `./data/`;
-(6) **Table2.R** main script for summary statistics for each of the four manually defined cell clusters in the HSCT dataset; 
-(7) **fig3.R** main script for pairwise scatter plots of the four fluorescent markers in the HSCT dataset; 
-(8) **Table3.R** main script for cluster-label cross-tabulation comparing the clustering results of GMM and FM-TMVN with the true class labels in the HSCT dataset; 
-(9) **fit_stone.R** main script for model fitting (GMM, FM-TMVN) with g=1~6 to the Stone Flakes dataset; saved results are stored in `./data/`;
-(10) **fig4.R** main script for BIC values for the GMM and FM-TMVN models with the number of clusters varying from g = 1 to 6; 
-(11) **fig5.R** main script for projected scatter plots overlaid with fitted GMM and FM-TMVN contours for 5 pairs of variables of the stone flakes data;
-(12) **TableF1.R** main script for Parameter estimates and their standard errors in parentheses for the GMM and FM-TMVN models fitted to the HSCT dataset; 
-(13)  **figF1.R** Scatter-histogram plots of one simulation case with 600 random samples and 10% missingness. The cross symbols represent kNN imputed missing values.

# Data and Code for data #
## Subfolder: ./data ##
`./data/` contains 

-(1) **hsct_result.RData** : These files contains precomputed fits of the GMM and FM-TMVN models with g = 4 on the HSCT dataset, enabling direct generation of all related tables and figures without re-running the estimation;
-(2) **stone_result.Rdata** : Contains precomputed fits of the GMM and FM-TMVN models for g = 1–6 on the Stone Flakes dataset, allowing you to generate all corresponding figures and tables without re-running the ECM algorithms;

`./data/source` subfolder contains
	
-(3) **hsct.csv** contains the HSCT flow cytometry data used for clustering and imputation analyses in Section 6;
-(4) **StoneFlakes.csv** contains the Stone Flakes morphological measurements used for clustering experiments in Section 6.

# Data and Code for results #
## Subfolder: ./results ##
`./results`  contains

-(1) **fig1.eps**: 3D perspective and contour plots for the FM-TMVN model with g = 1~3;
-(2) **fig2.eps**: RMSE of estimated parameter values across varying samples sizes and missing rates based on simulation study;
-(3) **Table1.csv**: Simulation results showing the STD and IMSE of the maximum likelihood estimates under missing rates of r = 0%, 10%, 20% and 30%, across larger sample sizes n = 300, 600, 1200, 1800;
-(4) **Table2.csv**: Summary statistics for each of the four manually defined cell clusters in the HSCT dataset;
-(5) **fig3.eps**: Pairwise scatter plots of the four fluorescent markers in the HSCT dataset;
-(6) **Table3.csv**: Cluster-label cross-tabulation comparing the clustering results of GMM and FM-TMVN with the true class labels in the HSCT dataset;
-(7) **fig4.eps**: BIC values for the GMM and FM-TMVN models with the number of clusters varying from g = 1 to 6 in the StoneFlakes dataset;
-(8) **fig5.eps**: Scatter plots overlaid with fitted GMM and FM-TMVN contours for 5 pairs of variables of the stone flakes data;
-(9) **TableF1.csv**: Parameter estimates and their standard errors in parentheses for the GMM and FM-TMVN models fitted to the HSCT dataset;
-(10) **figF1.eps**: Scatter-histogram plots of one simulation case with 600 random samples and 10% missingness from FM-TMVN model;

`./data/results` subfolder contains
	
-(11) `./simulation/`: Subfolder storing intermediate simulation results, including estimated parameters (`./para/`), standard errors (`./sd/`) and RDatas (`./rdata/`) for various (n, missing rate) settings.

## Additional Remark ##
- Note (1): One can directly run each "source(.)" described in **master.R** file in the seperate R session to obtain the results;
- Note (2): To draw the 3D perspective and contour plots as shown in Figure 1, please  run the **fig1.R** scripts in subfolder `./code/`. The results have been stored in `./results/`.
- Note (3): R code **simulation.R** generates the intermediate results of Tables 1 and Figure 2 in the manuscript. Because the code requires a large amount of computation time, we have pre-saved the intermediate results under the `./results/simulation/` subfolders, including `./para/` for parameter estimates, `./sd/` for standard errors and `./rdata/` for RDatas.
- Note (4): To reproduce the results presented in Figure 2, please source the estimation results stored in the `./results/simulation/para/` subfolders, and then run the script **fig2.R** in the `./code/` folder. The figure will be saved as **fig2.eps** in the `./results/` folder.
- Note (5): To reproduce the results presented in Table 1, please source the estimation results stored in the `./results/simulation/para/` and `./results/simulation/sd/` subfolders, and then run the script **SimTable1.R** in the `./code/` folder. The final evaluation table will be saved as **SimTable1.csv** in the `./results/` folder.
- Note (6): Because the estimation procedures in **fit_hsct.R** involve fitting two different models across high dimensions, the computations are time-consuming. Therefore, we have saved the fitted results in **hsct.result.Rdata** under the `./data/` folder. Users can run subsequent scripts to directly reproduce the table.
- Note (7): To draw Table 2 and Figure 3 in this paper, please run the **Tab2.R**,**fig3.R** scripts in subfolder `./code/`.
- Note (8): To draw Tables 3 in this paper, please load the **hsct_result.RData** file in subfolder `./data/`, and then run the **Table3.R** scripts in subfolder `./code/`. 
- Note (9): Because the estimation procedures in **fit_stone.R** involve fitting two different models across g=1 to 6, the computations are time-consuming. Therefore, we have saved the fitted results in **stone.result.Rdata** under the `./data/` folder.
- Note (10): To draw Figures 4 and 5 in this paper, please load the **stone.result.Rdata** file in subfolder `./data/`, and then run the **fig4.R** and **fig5.R** scripts in subfolder `./code/`.
- Note (11): To reproduce the results presented in Table F1, please source the scripts **SE.FMTMVN.na.R** from the `./function/` folder, load the **hsct_result.RData** file in subfolder `./data/`, and then run the script **figF1.R** in the `./code/` folder. The figure will be saved as **figF1.eps** in the `./results/` folder.
- Note (12): To reproduce the results presented in Figure F1, please source the scripts **gener.na.R** from the `./function/` folder, and then run the script **figF1.R** in the `./code/` folder. The figure will be saved as **figF1.eps** in the `./results/` folder.
  
