clear all
clc

%==============
% For 1000 run

rng('default')

nRun = 1000;
dim = 100;

w = rand(nRun, dim);
x = rand(nRun, dim);

addpath_toolbox

tic
val = Compute_GST_EXP1_vec_V3(x, w)
toc

