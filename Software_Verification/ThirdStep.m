%% Description
%
% Before running this script, please execute setup.m first to set up path
% and folder parameters. 
%
% Using the kNN Classifier trained in Step 2, Step 3 is the prediction of
% labels for the unseen, testing part of the data.
%
% It will load the feature data and preprocess it. Then will load the model
% as well and run the prediction algorithm.
%
% Author: Tamas Kiss <kiss.t@wigner.hu>

%% Parameters

%IO
par.DataFilt = 'A*'; % Filter to select which set of data to work on
par.OneModelFile = fullfile(par.IntDir, 'ExampleRatSingleModel.mat'); %Matlab file to store kNN classifier
par.MultipleModelFile = fullfile(par.IntDir, 'ExampleRatMultipleModels.mat'); %Matlab file to store kNN classifier
par.PerfLogFile = fullfile(par.IntDir, 'PerformanceLogsRat.txt');
par.TestPF = '_test.mat'; % This will be appended to the partition put aside for testing

%Testing
par.TrainLabels = {'W', 'NR', 'R'}; %Manual labels to train
par.ROCNameMap = {... %Used to plot short labels on plot presenting model acuracy
    'A1_VEH_021815' 'A1';...
    'A2_VEH_022415' 'A2';...
    'A3_VEH_022715' 'A3';...
    'A4_VEH_051915' 'A4';...
    'A5_VEH_052615' 'A5';...
    'A6_VEH_052915' 'A6';...
    'A7_VEH_052615' 'A7'};
% par.ROCNameMap = {... %Used to plot short labels on plot presenting model acuracy
%     'tko116_01292016' '116';...
%     'tko118_01182016' '118';...
%     'tko120_01112016' '120';...
%     'tko122_01072016' '122';...
%     'tko124_01052016' '124';...
%     'tko128_01182016' '128'};
par.DeflateVal = 0.9; %90% of W and NR epochs are removed to compensate for too few REM epochs

%% Open log file
diary(par.PerfLogFile)

%% Load pre-calculated feature data
clear T
files = dir([par.IntDir filesep par.DataFilt par.TestPF]);
for expidx = 1:length(files)
    ifn = [par.IntDir filesep files(expidx).name];
    fprintf('Loading %s.\n', ifn)
    tmp = load(ifn);
    fn = canonize_fieldname(strrep(files(expidx).name, par.TestPF, ''));
    T.(fn) = tmp.T;
    epdur.(fn) = tmp.epdur;
    exps{expidx} = fn; %#ok<SAGROW>
    OriginalLabels.(fn) = T.(fn){:, 'Scores'};
end
fprintf('\n')
clear expidx tmp ifn fn files

%% Feature preprocessing
TNormd = preprocess_features(T, 'EpDur', epdur, 'IsTrainingSet', true);
TNonNormd = preprocess_features(T, 'DoNormalize', false, ...
    'EpDur', epdur, 'IsTrainingSet', true);

%% Load pre-trained single model
fprintf('Loading model from %s...', par.OneModelFile)
Single = load(par.OneModelFile);
fprintf(' done.\n\n')

%% Predict labels using single model
%Structure SingleLabels will store predicted labels when single model is
%used. Each animal is described by a field in this struch with a cell of
%strings array, same length as OriginalLabels.
fprintf('Predicting labels...')
clear labels
for expidx = 1:length(exps)
    SingleLabels.(exps{expidx}) = Single.mdl.only.predict(TNormd.(exps{expidx}){:, [false Single.seld_features]});
end
fprintf(' done.\n\n')
clear expidx

%% Plot and print a quick evaluation of estimation goodness with single model
mg = evaluate_model_goodness(SingleLabels, OriginalLabels, ...
    'PlotROC', par.TrainLabels, 'ROCNameMap', par.ROCNameMap, ...
    'ROCLayout', [1,3]);
set(gcf, 'Name', 'Prediction using single model')

%Print performance metrics
fprintf('True positive rates and false positive rates achieved using a single model for all animals:\n')
PM = NaN(length(exps), 2, length(par.TrainLabels));
for expidx = 1:length(exps)
    fprintf('%s:\n', exps{expidx})
    for lidx = 1:length(par.TrainLabels)
        tmpSEN = mg.(exps{expidx}).PerfMetr{par.TrainLabels{lidx}, 'SEN'};
        tmpFPR = mg.(exps{expidx}).PerfMetr{par.TrainLabels{lidx}, 'FPR'};
        PM(expidx, :, lidx) = [tmpSEN tmpFPR];
        fprintf('\t%s -- TPR: %0.2f; FPR: %0.2f\n', ...
            par.TrainLabels{lidx}, tmpSEN, tmpFPR);
    end
    fprintf('\n')
end
for lidx = 1:length(par.TrainLabels)
    fprintf('Average +- std of %s [SEN, FPR] = [%0.2f, %0.2f] +- [%0.2f, %0.2f].\n',...
        par.TrainLabels{lidx}, mean(PM(:, :, lidx)), std(PM(:,:, lidx)));
end
fprintf('\n\n')
clear mg expidx lidx tmp* PM

%% Load pre-trained multiple models
fprintf('Loading model from %s...', par.OneModelFile)
Multiple = load(par.MultipleModelFile);
fprintf(' done.\n\n')

%% Predict labels using multiple models
fprintf('Predicting labels...')
clear labels
for expidx = 1:length(exps)
    MultipleLabels.(exps{expidx}) = Multiple.mdl_set.(exps{expidx}).predict(TNonNormd.(exps{expidx}){:, [false Multiple.seld_features.(exps{expidx})]});
end
fprintf(' done.\n\n')
clear expidx

%% Plot and print a quick evaluation of estimation goodness with multiple model
mg = evaluate_model_goodness(MultipleLabels, OriginalLabels, ...
    'PlotROC', par.TrainLabels, 'ROCNameMap', par.ROCNameMap, ...
    'ROCLayout', [1,3]);
set(gcf, 'Name', 'Prediction using multiple models')

%Print performance metrics
fprintf('True positive rates and false positive rates achieved using one model for each animal:\n')
PM = NaN(length(exps), 2, length(par.TrainLabels));
for expidx = 1:length(exps)
    fprintf('%s:\n', exps{expidx})
    for lidx = 1:length(par.TrainLabels)
        tmpSEN = mg.(exps{expidx}).PerfMetr{par.TrainLabels{lidx}, 'SEN'};
        tmpFPR = mg.(exps{expidx}).PerfMetr{par.TrainLabels{lidx}, 'FPR'};
        PM(expidx, :, lidx) = [tmpSEN tmpFPR];
        fprintf('\t%s -- TPR: %0.2f; FPR: %0.2f\n', ...
            par.TrainLabels{lidx}, tmpSEN, tmpFPR);
    end
    fprintf('\n')
end
for lidx = 1:length(par.TrainLabels)
    fprintf('Average +- std of %s [SEN, FPR] = [%0.2f, %0.2f] +- [%0.2f, %0.2f].\n',...
        par.TrainLabels{lidx}, mean(PM(:, :, lidx)), std(PM(:,:, lidx)));
end
fprintf('\n\n')
clear mg expidx lidx tmp* PM

%% Extra step: Deflate single model
%Identify NR and W epochs in model and remove DeflateVal fraction of them
Wloc = find(strcmp(Single.mdl.only.Y, 'W'));
NRloc = find(strcmp(Single.mdl.only.Y, 'NR'));
DelW = round(length(Wloc)*par.DeflateVal);
DelNR = round(length(NRloc)*par.DeflateVal);
X = Single.mdl.only.X;
Y = Single.mdl.only.Y;
if DelW > 0 && DelNR > 0 && DelW < length(Wloc) && DelNR < length(NRloc)
    DelWIdx = datasample(Wloc, DelW, 'Replace', false);
    DelNRIdx = datasample(NRloc, DelNR, 'Replace', false);
    X([DelWIdx; DelNRIdx], :) = [];
    Y([DelWIdx; DelNRIdx], :) = [];
end
fprintf('%i W and %i NR were deleted. Originally %i observations -> %i new observations.\n', ...
    DelW, DelNR, Single.mdl.only.NumObservations, length(Y))

%Put remaining values back into model
fprintf('Retraining classifier...')
NewMod = ClassificationKNN.fit(X, Y, 'NumNeighbors', Single.mdl.only.NumNeighbors, ...
    'Cost', Single.mdl.only.Cost, 'NSMethod', 'exhaustive', 'PredictorNames', Single.mdl.only.PredictorNames, ...
    'ResponseName', 'Score', 'DistanceWeight', Single.mdl.only.DistanceWeight);
fprintf(' done.\n')

% Predict labels using deflated model
fprintf('Predicting labels...')
for expidx = 1:length(exps)
    DeflatedLabels.(exps{expidx}) = NewMod.predict(TNormd.(exps{expidx}){:, [false Single.seld_features]});
end
fprintf(' done.\n\n')
clear Wloc NRloc DelW DelNR DelWIdx DelNRIdx expidx

%% Plot and print what happened to prediction goodness
mg = evaluate_model_goodness(DeflatedLabels, OriginalLabels, ...
    'PlotROC', par.TrainLabels, 'ROCNameMap', par.ROCNameMap, ...
    'ROCLayout', [1,3]);
set(gcf, 'Name', 'Prediction using deflated single model')

%Print performance metrics
fprintf('True positive rates and false positive rates achieved using deflated single model:\n')
PM = NaN(length(exps), 2, length(par.TrainLabels));
for expidx = 1:length(exps)
    fprintf('%s:\n', exps{expidx})
    for lidx = 1:length(par.TrainLabels)
        tmpSEN = mg.(exps{expidx}).PerfMetr{par.TrainLabels{lidx}, 'SEN'};
        tmpFPR = mg.(exps{expidx}).PerfMetr{par.TrainLabels{lidx}, 'FPR'};
        PM(expidx, :, lidx) = [tmpSEN tmpFPR];
        fprintf('\t%s -- TPR: %0.2f; FPR: %0.2f\n', ...
            par.TrainLabels{lidx}, tmpSEN, tmpFPR);
    end
    fprintf('\n')
end
for lidx = 1:length(par.TrainLabels)
    fprintf('Average +- std of %s [SEN, FPR] = [%0.2f, %0.2f] +- [%0.2f, %0.2f].\n',...
        par.TrainLabels{lidx}, mean(PM(:, :, lidx)), std(PM(:,:, lidx)));
end
clear mg expidx lidx tmp* PM

%% Close log file
diary off
