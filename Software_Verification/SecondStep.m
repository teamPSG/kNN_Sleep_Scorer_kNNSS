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
par.DataFilt = 'A*'; % Filter to select which set of data to work on
par.OneModelFile = fullfile(par.IntDir, 'ExampleRatSingleModel.mat'); %Matlab file to store kNN classifier
par.MultipleModelFile = fullfile(par.IntDir, 'ExampleRatMultipleModels.mat'); %Matlab file to store kNN classifier
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

% The goal of this section here is to collect information from all animals
% and thus increase the size of our training set. This model will be able
% to generalize more than the next version, however, our experience shows
% that in some cases it will not be as accurate as training a model for
% each animal.

files = dir([par.IntDir filesep par.DataFilt par.TrainPF]);
train_one_model(par.IntDir, {files.name}, par.OneModelFile, ...
    'ValMeth', par.ValMeth, 'cost_m', par.CostM, 'Remboost', par.REMBoost, ...
    'CombineStates', par.CombineStates, 'RemoveStates', {}, ...
    'PerLen', par.WakeLength, 'PerStart', par.WakeStart);
clear files

%% Train one classifier for each animal

% In this version each animal has its own model. This can be more accurate
% than a general model but since it contains less data points it is not as
% good in generalization. This means on one hand that it won't work that
% well on data from other animals or on the same animal in very different
% conditions.

files = dir([par.IntDir filesep par.DataFilt par.TrainPF]);
train_many_models(par.IntDir, {files.name}, par.MultipleModelFile, ...
    'ValMeth', par.ValMeth, 'postfn', '', 'RemFromFn', par.TrainPF);
clear files ans
