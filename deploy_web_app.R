rm(list=ls())
library(rsconnect)

rsconnect::setAccountInfo(name='', token='', 
                          secret='')

rsconnect::deployApp('/Users/estebantisnes/Documents/GitHub/ARIMA_Simulation_LAB/', 
                     account = "tisnesesteban")


