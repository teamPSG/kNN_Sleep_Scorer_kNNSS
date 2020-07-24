%% Description
%
% Before running this script, please execute setup.m first to set up path
% and folder parameters. 
%
% This script is the first step to take for automated sleep scoring with
% kNNSS. It will first set the path to include functions called during the 
% process then it calculates and saves the feature table for a set of data
% files.
%
% EDF files are assumed to be in the par.DataDir folder. All EDF files in
% this folder will be processed and feature tables saved in new .mat files
% as specified below.
%
% Files containing feature tables for each experiment are stored in
% par.IntDir with par.PostFn appended to the file name and extension
% changed from .edf to .mat.
%
% As a last step features for each experiment are partitioned into a
% training and a test set and saved.
%
% Author: Tamas Kiss <kiss.t@wigner.hu>

%% Parameters                                                              \par\

%IO
par.DataFilt = 'A*'; % Filter to select which set of data to work on
par.PostFn = '_trset.mat'; %This part will be appended to the name of the EDF file and features will be saved in this file for each experiment.
par.DataExt = '.edf'; %File name extension for raw data
par.TestPF = '_test.mat'; % This will be appended to the partition put aside for testing
par.TrainPF = '_train.mat'; % This will be appended to the partition put aside for training the classifier

%Partitioninig
par.CombineStates = {...
    'W'  'WA';...
    'NR' 'NA';...
    'R'  'RA';...
    'A'  'D';...
    'A'  '-'};
par.TrainLabels = {'W', 'NR', 'R'}; %Manual labels to train

%% Calculate features
files = dir([par.DataDir filesep par.DataFilt par.DataExt]);
for file = files'
    filename = strrep(file.name, par.DataExt, '');
    fprintf('Processing %s.\n', filename);
    generate_statespace(par.TopDir, filename, 'Description', filename, ...
        'IsTrainingSet', true, ...
        'OutDir', par.IntDir, 'WriteTmpFile', [filename par.PostFn]);
end
clear file files filename ans

%% Load pre-calculated features                                            \ED.(exp); T.(exp); exps{eidx}; fns{eidx}\
clear T
files = dir([par.IntDir filesep par.DataFilt par.PostFn]);
for expidx = 1:length(files)
    ifn = [par.IntDir filesep files(expidx).name];
    fprintf('Loading %s.\n', ifn)
    tmp = load(ifn);
    expname = canonize_fieldname(strrep(files(expidx).name, par.PostFn, ''));
    ED.(expname) = tmp.epdur;
    T.(expname) = tmp.T;
    exps{expidx} = expname; %#ok<SAGROW>
    fns{expidx} = strrep(files(expidx).name, par.PostFn, ''); %#ok<SAGROW>
end
fprintf('\n')
clear expidx tmp expname ifn files 

%% Combine score labels

% Since the standard scoring time window is 10 or 30 seconds long, a number
% of states might be present in a single window. Scoring experts often
% times mark these by scoring an epoch as an artifactual one, like NA here
% stands for a non-REM sleep epoch that is not really characteristic to the
% non-REM state. These could be ignored altogether, however, to increase
% sample size, in this example we merge marked epochs with clean ones.

fprintf('Combining state labels...')
for expidx = 1:length(exps)    
    [idx, loc] = ismember(T.(exps{expidx}){:, 'Scores'}, par.CombineStates(:,2));    
    T.(exps{expidx})(idx, 'Scores') = par.CombineStates(loc(idx), 1);    
end
fprintf(' done.\n')
clear expidx idx loc

%% Partition data into training and test sets                              \T_train.(exp); T_test.(exp)\
clear T_test T_train
for expidx = 1:length(exps)
    fprintf('Partitioning %s...\n', exps{expidx})
    fprintf('\tSlecting training labels:')
    train_index = [];
    for lidx = 1:length(par.TrainLabels)
        fprintf(' %s', par.TrainLabels{lidx})
        tmp = find(ismember(T.(exps{expidx}){:, 'Scores'}, par.TrainLabels{lidx}));
        train_index = union(train_index, randsample(tmp, round(length(tmp)/2)));
    end
    test_index = setdiff(1:height(T.(exps{expidx})), train_index);
    
    T_train.(exps{expidx}) = T.(exps{expidx})(train_index, :);
    T_test.(exps{expidx}) = T.(exps{expidx})(test_index, :);
    fprintf('\nDone.\n\n')
end
fprintf('\n')
clear expidx C train_index test_index lidx tmp

%% Save training and test sets
clear T
for expidx = 1:length(exps)
    ofn = [par.IntDir filesep fns{expidx} par.TrainPF];
    fprintf('Saving data for %s\n', exps{expidx})
    % Training set
    fprintf('\t%s\n', ofn)
    T = T_train.(exps{expidx}); %#ok<NASGU>
    epdur = ED.(exps{expidx}); %#ok<NASGU>
    save(ofn, 'T', 'epdur')
    % Test set
    ofn = [par.IntDir filesep fns{expidx} par.TestPF];
    fprintf('\t%s\n', ofn)
    T = T_test.(exps{expidx}); %#ok<NASGU>
    save(ofn, 'T', 'epdur')
end
fprintf('\n')
clear expidx ofn T epdur T_test T_train ED fns exps
