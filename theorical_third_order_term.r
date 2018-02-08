# 
set.seed(123)
#
library(dplyr)
library(data.table)
library(foreach)
library(doParallel)
library(doSNOW)
###############################################
rho = 0.8
mu = 3.0
sigma = 0.030
R=0.18
DT=30
m=200
h.max = 0.016
h=0.010
#
# pcf grid
n.grid= 500
grid=seq(0,R,length.out = n.grid)
gridt=R/n.grid
######   apprximate the Q1&Q2 ###############
# calculus grid
n_grid_c = 8000
bin_grid_c   = DT/n_grid_c 
grid_c= seq(bin_grid_c,DT,by = bin_grid_c)
calculus_bin_area = bin_grid_c^2