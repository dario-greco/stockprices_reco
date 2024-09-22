# Replication files for "Improving out-of-sample forecasts of stock price indexes with forecast reconciliation and clustering"

Here you can find the files to replicate the results of the paper "Improving out-of-sample forecasts of stock price indexes with forecast reconciliation and clustering", written by Raffaele Mattera, George Athanasopoulos and Rob Hyndman.

The data needed for replicating the results are in the main folder. In particular:

- "DataP0.xlsx" includes the bottom stock prices considered for the main application to Dow Jones Index;
- "DJIA - composition at 2020.xlsx" includes the information used for metadata clustering, that is, the industry sector and the Exchange where the stocks are traded;
- "DataP0_2015.xlsx" includes the bottom stock prices considered for robustness to Dow Jones Index in the period 2015-2017;
- "DowJones2015.xlsx" includes information used for metadata clustering of the 2015-2017 dataset;
- "sp500_companies.xlsx" includes the bottom stock prices considered for the robustness with S&P500 data;
- "sp500_stocks.csv" includes information used for metadata clustering of the S&P500 data.

The main results are obtained from the following R codes, included in the folder "R_scripts" (R version: 4.4.1):

- "Reconciliation_clustering" to replicate the clustering results used in the main application (Section 4);
- "Reconciliation_forecasting" to obtain base forecasts and reconcilied forecasts (Section 5);
- "Reconciliation_accuracy" to evaluate out of sample forecast accuracy of the methods (Sections 5.1 and 5.2);
- "Reconciliation_timing" to replicate the results of the Section 5.3.

The sub-folder called "Robustness analysis" ("R_scripts/Robustness analysis") contains the codes to replicate the results of the Section 5.4 by means of the two codes:

- "DJI2015appl" considering the data of the Dow Jones Index between 2015-2017 and 
- "SP500appl" for the stocks in the S&P500 Index.

We use the R version 4.4.1 with the following packages:

hts_6.0.2, forecast_8.23.0, TSclust_1.3.1, cluster_2.1.6, foreach_1.5.2, rio_1.1.1, BatchGetSymbols_2.6.4,
dplyr_1.1.4, sandwich_3.1-0, lmtest_0.9-40, doParallel_1.0.17, ggplot2_3.5.1, here_1.0.1, iterators_1.0.14, reshape2_1.4.4 

The folder called "Figures" is empty, but it is needed to save figures obtained with the R codes.

To reproduce the results correctly, we suggest to unzip the file and running the R codes directly from the resulting folder. For any enquiry, please write to Raffaele Mattera (email: raffaele.mattera@uniroma1.it)
