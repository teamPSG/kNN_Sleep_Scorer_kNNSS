function [ mdl_set, seld_features, normfact ]  = train_many_models(drf, ts, mf, varargin)

%This function fits many models individually.
%
%Usage:
%  [Model SelectedFeatures, NormalizingFactor]= train_one_model(DataRootFolder, TrainigSet, ModelFile, ...)
%
%Input arguments:
%  DataRootFolder: string, absolute path for data location. Uses the format
%    inherited from the first SRI dataset, the function will look for a
%    sub-folder called FFT and one called EDF under this folder to get
%    scoring, and raw data, respectively.
%  TrainingSet: 1xN cellstr, a list of file names describing which
%    experiment to use for training the machine. All files (including one
%    in the training set) are used for testing.
%  ModelFile: string, full path to a Matlab file to store the KNN models
%    and its associated selected features.
%
%Output:
%  Model: struct of ClassificationKNN, the models created for the training
%    sets
%  SelectedFeatures: 1xM logical array, where M is the number of features
%    in calculated. Value is true for features used in building the model,
%    false othrewise.
%  NormalizingFactor: structure containing normalizing factors for each
%    features and for all experiments processed.
%
%Optional input arguments:
%  'NN': double, number of neighbors used in the kNN classifier. Default is
%    11.
%  'DistanceWeight': string, method used by kNN fitting to take distance of
%    points into account. See Matlab documentation for ClassificationKNN
%    for details. Default is 'squaredinverse'.
%  'FSDir': string, specifies direction of sequential feature selection.
%    Can be 'forward' and 'backward'. Default is 'forward'.
%  'cost_m'. string, specifies the type of cost matrix used in classifier.
%    Can be 'offdiag' or 'rempref'. Offdiag poses same penalty for
%    misclassifying any states, repref increases cost of REM
%    misclassification. This is useful since REM is highly underrepresented
%    in typical datasets. Default is 'offdiag'.
%  'REMboost': double, the value of REM misclassification cost. Default is
%    2.
%  'FeatFile': string, specifies a file that specifies features to be used
%    by the classifier. Default is empty.%
%  'FSKeep': string or cell of string, specifies what features to keep
%    during the feature selection procedure. If set to 'All', every
%    features will be used in the model. If 'Free' sequential feature
%    selection is performed allowing inclusion/exclusion of any features.
%    If 'Suggested', overlap of feature distributions in defferent states
%    will be calculated and most dissimilar features will be included plus
%    any additional features selected by sequential feature selection. If a
%    cell of string is specified, named features will be used plus any
%    other features selected by sequential feature selection. Default is
%    'Suggested'.
%  'UseFeat': cell of strings, explicitely specifies which features to use
%    in the model, omits sequential feature selection. Default is {}.
%  'OmitFeat': cell of strings, a list of features not to use in the model.
%    Default is {}.
%  'ValMeth': double of string, if a number N is given, N-fold cross
%    validation is performed in sequential feature selection.  If the
%    string 'resub' is given, resubstitution validation is performed.
%    Default is 5 to perform 5-fold cross validation.
%  'Verbose': boolean, if true the forward selection procedure displays its
%    history. Default is true.
%  'PerLen': double, the length of the period used for calculating feature
%    normalization factors. Given in minutes. Default is 30.
%  'PerStart': double, start time of the period used for normalization,
%    given in minutes. Default is the beginning of the experiment: 0.
%  'DoNormalize': boolean, turns on/off normalization to selected period.
%    Default is TRUE.
%  'NormFact': struct, stores normalizing factors for each experiment in
%    TrainingSet. Default is empty structure.
%  'NormMap': Nx2 cell of strings, specifies a correspondance between names
%    in TrainingSet and names in NormFact. Default is empty cell.
%  'Norm1st': boolean, if true normalization preceeds log transformation of
%    features, if false order is swapped. Default is true.
%  'NormMeth': string, specifies normalization method. Can be 'time': the
%    interval specified by 'PerStart' and 'PerLength' is used to calculate
%    normalization factors; 'state': if manual labels are also available
%    'NormStates' labels will be used to calculate normalizing factors;
%    'zscore': full length of data is used for zscoring. Default is 'time'.
%  'NormStates': cell of strings, specifies what state labels to include if
%    NormMeth is 'state'. Default is 'time'.
%  'NormAvg': string, the name of the function to be used for
%    normalization. Can be 'mean' or 'median'. Default is 'median'.
%  'postfn': string, the value of this argument is appended to the file
%    names specified in TrainingSet. Default value is empty string.
%  'ArtifMap': structure, containing a boolean for each experiment. Epochs
%    with TRUE ArtifMap values will be discarded from training. Default is
%    empty structure.
%  'AFMapUse': string, specifies how to combine artifacts detected by 
%    different methods. 'union': all epochs marked as artifact by automated
%    outlier detection and manual input will be marked as artifact;
%    'intersect': only epochs marked as artifact by both methods will be
%    marked as artifact; 'alone': only manually identified artifactual
%    epochs are marked as artifacts. Default is 'union'.
%  'RemoveStates': cell of strings. Epochs labelled with labels specified
%    here will be removed from the training set. Default is {'D' 'SS' 'U'}
%  'OutlierVar': cell of strings, specifies which variables in the feature
%    set are used for outlier detection. Default is {'emg_RMS' 'eeg_RMS'}.
%  'MisSigFeat': string, used to specify which signal channel to look at
%    to find epochs with missing signal in them. Default is 'eeg_RMS'.
%  'CombineStates', Nx2 cell of strings, specifies what states to combine.
%    First column is the label that will be finally assigned to labels
%    originally labelled by labels in the second column. Default is: 
%    {'W' 'WA'; 'NR' 'NA'; 'R' 'RA'}.
%  'RemFromFn': string, this part of the file name used for training data
%    will be removed from the field name of the model structure.
%See also postprocess_features, estimate_feature_goodness, select_features,
%train_classifier, and evaluate_model_goodness
%
%Author: Tamas Kiss <kiss.t@wigner.hu>

%% Parameters
p = inputParser;
%IO
addRequired(p, 'drf', @isstr)
addRequired(p, 'ts', @iscell)
addRequired(p, 'mf', @isstr)
addParamValue(p, 'postfn', '_trset.mat', @isstr); %#ok<*NVREPL> %File name for intermediate data
%Feature post-processing
addParamValue(p, 'ArtifMap', struct([]), @isstruct);
addParamValue(p, 'AFMapUse', 'union', @isstr); %union, intersect, alone
addParamValue(p, 'RemoveStates', {'D' 'SS' 'U'}, @iscell);
addParamValue(p, 'OutlierVar', {'emg_RMS' 'eeg_RMS'}, @iscell);
addParamValue(p, 'MisSigFeat', 'eeg_RMS', @isstr); %This feature is used to determine missing signal locations
addParamValue(p, 'CombineStates', {...
    'W' 'WA';...
    'NR' 'NA';...
    'R' 'RA'}, @iscell);
%Feature selection
addParamValue(p, 'FeatFile', '', @isstr);
addParamValue(p, 'NN', 3, @isnumeric);
addParamValue(p, 'DistanceWeight', 'squaredinverse', @isstr);
addParamValue(p, 'FSDir', 'forward', @isstr);
addParamValue(p, 'cost_m', 'OffDiag', @isstr); %can be 'offdiag' or 'REMpref'
addParamValue(p, 'REMboost', 2, @isnumeric); %If REMpref this is the number used instead of 1 for REM penalty
addParamValue(p, 'FSKeep', 'Suggested', @isstr);%'All';%{'eeg_RMS' 'emg_RMS' 'eeg_theta' 'eeg_delta' 'eeg_gamma'};
addParamValue(p, 'UseFeat', {}, @iscell); %{'eeg_delta_high' 'eeg_theta' 'eeg_gamma' 'eeg_RMS' 'emg_gamma' 'emg_RMS'};
addParamValue(p, 'OmitFeat', {}, @iscell); %{'eeg_delta','eeg_delta_low','eeg_delta_high','eeg_theta','eeg_alpha','eeg_sigma','eeg_beta','eeg_gamma'};
addParamValue(p, 'ValMeth', 5, @(x)validateattributes(x,{'char', 'numeric'}, {'nonempty'}));
%Etc
addParamValue(p, 'Verbose', false, @islogical);
addParamValue(p, 'PerLen', 30, @isnumeric); %This is given in minutes
addParamValue(p, 'PerStart', 0, @isnumeric); %Also in minutes
addParamValue(p, 'DoNormalize', false, @islogical);
addParamValue(p, 'NormFact', struct([]), @isstruct);
addParamValue(p, 'NormMap', {}, @iscell);
addParamValue(p, 'Norm1st', true, @islogical);
addParamValue(p, 'NormMeth', 'time', @isstr);
addParamValue(p, 'NormStates', {}, @iscell);
addParamValue(p, 'NormAvg', 'median', @isstr);
addParamValue(p, 'RemFromFn', '', @isstr);

parse(p, drf, ts, mf, varargin{:});

%% Load pre-calculated features
clear T
for fidx = 1:length(ts)
    filename = [drf filesep ts{fidx} p.Results.postfn];
    fprintf('Loading feature set from %s.\n', filename)
    tmp = load(filename);
    fn = canonize_fieldname(strrep(ts{fidx}, p.Results.RemFromFn, ''));
    T.(fn) = tmp.T;
    epdur.(fn) = tmp.epdur;
end
exps = fieldnames(T);
fprintf('\n')
clear filename tmp fn

%% Preprocess features
[T, rmidx, ~, normfact] = preprocess_features(T, 'IsTrainingSet', true, 'EpDur', epdur, ...
    'RemoveStates', p.Results.RemoveStates, 'CombineStates', p.Results.CombineStates, ...
    'PerLen', p.Results.PerLen, 'PerStart', p.Results.PerStart, 'DoNormalize', p.Results.DoNormalize, ...
    'ArtifMap', p.Results.ArtifMap, 'AFMapUse', p.Results.AFMapUse, ...
    'NormFact', p.Results.NormFact, 'NormMap', p.Results.NormMap, 'Norm1st', p.Results.Norm1st, ...
    'NormMeth', p.Results.NormMeth, 'NormStates',p.Results.NormStates, 'NormAvg', p.Results.NormAvg);
for expidx = 1:length(exps)
    T.(exps{expidx}) = T.(exps{expidx})(~rmidx.(exps{expidx}),:);
end
clear expidx rmidx

%% Create cost matrices
fprintf('Starting cost matrix generation...')
clear cost_m
for expidx = 1:length(exps)
    tmp.states = unique(T.(exps{expidx}){:,'Scores'});
    tmp.nostates = length(tmp.states);
    tmp.REMloc = find(strcmp(tmp.states, 'R'));
    if strcmpi(p.Results.cost_m, 'offdiag')
        cost_m.(exps{expidx}) = ...
            ones(tmp.nostates) - eye(tmp.nostates);
    elseif strcmpi(p.Results.cost_m, 'rempref')
        if ~isempty(tmp.REMloc)
            tmp.cost_m = ones(tmp.nostates) - eye(tmp.nostates);
            tmp.cost_m(tmp.REMloc, :) = p.Results.REMboost;
            tmp.cost_m(:, tmp.REMloc) = p.Results.REMboost;
            tmp.cost_m(tmp.REMloc, tmp.REMloc) = 0;
            cost_m.(exps{expidx}) = tmp.cost_m;
        else
            fprintf('\nREM was not seen for experiment %s. Using OffDiag.\n', ...
                exps{expidx});
            cost_m.(exps{expidx}) = ...
                ones(tmp.nostates) - eye(tmp.nostates);
        end
    else
        fprintf('\nUnknown cost matrix name %s. Using OffDiag.\n', ...
            p.Results.cost_m);
        cost_m.(exps{expidx}) = ...
            ones(tmp.nostates) - eye(tmp.nostates);
    end
end
fprintf(' done\n\n')
clear expidx tmp

%% Select features based on training data for each animal
fprintf('Starting feature selection...\n')
clear seld_features
if isempty(p.Results.FeatFile)
    if isempty(p.Results.UseFeat) && strcmpi(p.Results.FSKeep, 'Suggested')
%         for expidx = 1:length(exps)
%             tmpT.(exps{expidx}) = T.(exps{expidx});
%         end
        [~, best_features] = estimate_feature_goodness(T);
    end
    for expidx = 1:length(exps)
        fprintf('********* %s:\n', exps{expidx})
        if isempty(p.Results.UseFeat)
            if strcmpi(p.Results.FSKeep, 'All')
                seld_features.(exps{expidx}) = ...
                    true(1,width(T.(exps{expidx}))-1);
            elseif strcmpi(p.Results.FSKeep, 'Suggested')
                fprintf('estimate_feature_goodness selected:\n')
                disp(best_features.(exps{expidx}))
                seld_features.(exps{expidx}) = select_features(...
                    T.(exps{expidx}), cost_m.(exps{expidx}),...
                    'NNeighbors', p.Results.NN, 'Verbose', true, 'Summary', true, 'SfSDir', p.Results.FSDir, ...
                    'KeepIn', best_features.(exps{expidx}), 'KeepOut', p.Results.OmitFeat, ...
                    'ValMeth', p.Results.ValMeth);
            elseif strcmpi(p.Results.FSKeep, 'Free')
                seld_features.(exps{expidx}) = select_features(...
                    T.(exps{expidx}), cost_m.(exps{expidx}), ...
                    'NNeighbors', p.Results.NN, 'Verbose', p.Results.Verbose, 'Summary', true, ...
                    'SfSDir', p.Results.FSDir, 'ValMeth', p.Results.ValMeth);
            else
                seld_features.(exps{expidx}) = select_features(...
                    T.(exps{expidx}), cost_m.(exps{expidx}),...
                    'NNeighbors', p.Results.NN, 'Verbose', true, 'Summary', true, 'SfSDir', p.Results.FSDir, ...
                    'KeepIn', p.Results.FSKeep, 'KeepOut', p.Results.OmitFeat, ...
                    'ValMeth', p.Results.ValMeth);
            end
        else
            seld_features.(exps{expidx}) = ...
                ismember(T.(exps{expidx}).Properties.VariableNames(2:end), ...
                p.Results.UseFeat);
            if sum(seld_features.(exps{expidx})) ~= length(p.Results.UseFeat)
                fprintf('train_many_models: Warning: some prescribed features were missed (%i/%i)\n', ...
                    sum(seld_features.(exps{expidx})), length(p.Results.UseFeat))
            end
        end
        fprintf('\n\n')
    end
else
    fprintf('Loading file: %s\n', p.Results.FeatFile);
    load(p.Results.FeatFile)
end
disp(seld_features)
fprintf('done.\n\n')
clear expidx best_features

%% Train models
fprintf('Starting model training...')
clear mdl_set
for expidx = 1:length(exps)
    tmp.name = exps{expidx};
    tmp.mdl = train_classifier(T.(tmp.name)(:, [true seld_features.(tmp.name)]), ...
        cost_m.(tmp.name), 'NNeighbors', p.Results.NN, 'DistanceWeight', p.Results.DistanceWeight);
        mdl_set.(tmp.name) = tmp.mdl;
end
fprintf(' done.\n\n')
clear tmp expidx

%% Save model and feature selection
if ~isempty(mf)
    fprintf('Saving model and associated feature selection in %s.\n', mf);
    save(mf, 'mdl_set', 'seld_features', 'normfact')
end

end

