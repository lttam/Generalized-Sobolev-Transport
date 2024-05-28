clear all
close all
clc

%============
% For 1 run

rng('default')
addpath_toolbox

% number of points.
N = 20;
dim = 5;

x = rand(dim, N);
y = rand(dim, N);

wx = rand(N, 1);
wx = wx/sum(wx);

wy = rand(N, 1);
wy = wy/sum(wy);

mu = wx;
nu = wy;

EPS = 1e-10;
c = sqdistance(x, y);
c = (c + c')/2;
c(c < EPS) = 0;

epsilon = 0.1;

options.niter = 1000;
options.tau = 0;
options.verb = 0;
tol = 1e-6;

phi = @(X) exp(X) - X - 1;
invphi = @(y) -Lambert_W(-exp(-y-1), -1)-y-1;

tic
val = OrliczWasserstein(phi, invphi, mu, nu, c, epsilon, tol)
toc





