%% Description
%
% Before running this script, please execute setup.m first to set up path
% and folder parameters. 
%
% Once the feature space of the data is generated and partitioned into a
% training and a test set, SecondStep.m will train and save a kNN
% classifier using the training set.
%
% Author: Tamas Kiss <kiss.t@wigner.hu>

%% Parameters                                                              \par\

%IO
par.OneModelFile = fullfile(par.IntDir, 'ExampleSingleModel.mat'); %Matlab file to store kNN classifier
par.MultipleModelFile = fullfile(par.IntDir, 'ExampleMultipleModels.mat'); %Matlab file to store kNN classifier
par.TrainPF = '_train.mat'; % This will be appended to the partition put aside for training the classifier

%Model fitting
par.CombineStates = {... %This variable combines some manual scores into one score
    'NR' 'NA'
    'A' 'D';...
    'A' '-';...
    'W' 'WA';...
    'R' 'RA'};
par.ValMeth = 4; %Sets model validation method. Integer N will do N-fold cross-validation; 'resub' performs resubstitution validation.
par.CostM = 'rempref'; %Cost matrix can be 'offdiag' for equal cost or 'rempref' for increased cost for REM misclassification
par.REMBoost = 2; %REM misclassification cost
par.WakeStart = 0; %minutes -- beginning of initial wake period used for normalization
par.WakeLength = 30; %minutes -- end of initial wake period

%% Train a common model merging training data from all animals
files = dir([par.IntDir filesep '*' par.TrainPF]);
train_one_model(par.IntDir, {files.name}, par.OneModelFile, ...
    'ValMeth', par.ValMeth, 'cost_m', par.CostM, 'Remboost', par.REMBoost, ...
    'CombineStates', par.CombineStates, 'RemoveStates', {}, ...
    'PerLen', par.WakeLength, 'PerStart', par.WakeStart);
clear files

%% Train one classifier for each animal
files = dir([par.IntDir filesep '*' par.TrainPF]);
train_many_models(par.IntDir, {files.name}, par.MultipleModelFile, ...
    'ValMeth', par.ValMeth, 'postfn', '', 'RemFromFn', par.TrainPF);
clear files