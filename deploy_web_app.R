
rm(list=ls())
library(rsconnect)

rsconnect::setAccountInfo(name='tisnesesteban', token='77768124592EB05F595F0A2CBE44636C', 
                          secret='cjVawps2PWok3PNgSeqefViRF+dCN0OgSe7Jd8TT')

rsconnect::deployApp('/Users/estebantisnes/Documents/GitHub/ARIMA_Simulation_LAB/', 
                     account = "tisnesesteban")

