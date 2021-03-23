clc
clear all
close all
addpath('CTISfunction')
addpath('data0318');
addpath('FFWTensor-master');

load('pandadata.mat');
load('Mu.mat');
load('Mu1.mat');
load('calipar_0620bird.mat');

%% load Mask
Mu2 = repmat(Mu,[1,1,par.k]);
par.Mu1 = Mu1;
par.Mu2 = Mu2;

%% input
par.data = g;

par.iter = 3;
fctEM = CodeEMf1(par);

par.iter = 3;
fctCFW = CodeLRf(par);


par.iter = 3;
[fctCPat]=CodePatf(par);


