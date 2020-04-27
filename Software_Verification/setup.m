%% Description
%This script sets up path and folders.
%
% Author: Tamas Kiss <kiss.t@wigner.hu>

%% Parameters
clear par

par.ROOTDIR = '/home/umat/bognor/kNNSS';
par.FunctionDir = fullfile(par.ROOTDIR,'Function_Library'); %Our functions live in this directory.
par.TopDir = fullfile(par.ROOTDIR, 'Example_Data', 'RawData'); %This is the top folder for raw data.
par.DataDir = fullfile(par.TopDir, 'EDF'); %This is the folder in which your EDF files live. All of them will be processed.
par.IntDir = fullfile(par.ROOTDIR, 'Example_Data', 'IntermRes'); %Training sets (feature tables in .mat files) will be saved in this folder.

%% Set the path
restoredefaultpath;
rmpath(userpath)
clear RESTOREDEFAULTPATH_EXECUTED
addpath(par.FunctionDir)
rehash

