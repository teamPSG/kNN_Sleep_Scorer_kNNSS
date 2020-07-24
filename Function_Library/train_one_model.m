function [ mdl, seld_features, normfact ] = train_one_model(drf, ts, mf, varargin)
%This function fits one kNN model on data from multiple subjects
%
%Usage:
%  [Model SelectedFeatures, NormalizingFactor] = ...
%             train_one_model(DataRootFolder, TrainigSet, ModelFile, ...)
%
%Input arguments:
%  DataRootFolder: string, absolute path for data location.
%  TrainingSet: 1xN cellstr, a list of file names describing which
%    experiment to use for training the machine.
%  ModelFile: string, full path to a Matlab file containing the kNN model
%    and its associated selected features.
%
%Output:
%  Model: struct of ClassificationKNN, the models created for the training
%    sets.
%  SelectedFeatures: structure containing boolean arrays that are true at
%    feature indices that were selected to be used in the classifier, false
%    otherwise.
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
%  'DoPCA': string, specifies if features of their Principal Components are
%    used for training the kNN classifier. Possible values are: 'fisrt'
%    (calculate PCA for features of the first animal and then project
%    subsequent animals into this PCA space); 'each' (project each animal's
%    feature space into its own PCA space); 'no' (do not calculate PCA, go
%    with features themselves). Default is 'no'.
%  'NPCs': double, number of principal components to use for model fitting.
%    Default is 3.
%  'PlotExplained': boolean, if true and 'DoPCA' is set to calculate PCA,
%    the variance explained by principal components is plotted. Default is 
%    false for not plotting.
%  'FeatFile': string, specifies a file that specifies features to be used
%    by the classifier. Default is empty.
%
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
addParamValue(p, 'postfn', '', @isstr); %#ok<*NVREPL> %File name for intermediate data
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
addParamValue(p, 'DoPCA', 'no', @isstr);
addParamValue(p, 'NPCs', 3, @isnumeric);
addParamValue(p, 'PlotExplained', false, @islogical);
addParamValue(p, 'FeatFile', '', @isstr);
addParamValue(p, 'NN', 11, @isnumeric);
addParamValue(p, 'FSDir', 'forward', @isstr);
addParamValue(p, 'cost_m', 'offdiag', @isstr); %can be 'offdiag' or 'REMpref'
addParamValue(p, 'REMboost', 2, @isnumeric); %If REMpref this is the number used instead of 1 for REM penalty
addParamValue(p, 'FSKeep', 'Suggested', @(x)ischar(x)||iscell(x)); %'All'; %'Free' %{'eeg_RMS' 'emg_RMS' 'eeg_theta' 'eeg_delta' 'eeg_gamma'};
addParamValue(p, 'UseFeat', {}, @iscell); %{'eeg_delta_high' 'eeg_theta' 'eeg_gamma' 'eeg_RMS' 'emg_gamma' 'emg_RMS'};
addParamValue(p, 'OmitFeat', {}, @iscell); %{'eeg_delta','eeg_delta_low','eeg_delta_high','eeg_theta','eeg_alpha','eeg_sigma','eeg_beta','eeg_gamma'};
addParamValue(p, 'ValMeth', 5, @(x)validateattributes(x,{'char', 'numeric'}, {'nonempty'}));
%Etc
addParamValue(p, 'Verbose', true, @islogical);
addParamValue(p, 'PerLen', 30, @isnumeric); %This is given in minutes
addParamValue(p, 'PerStart', 0, @isnumeric); %Also in minutes
addParamValue(p, 'DoNormalize', true, @islogical);
addParamValue(p, 'NormFact', struct([]), @isstruct);
addParamValue(p, 'NormMap', {}, @iscell);
addParamValue(p, 'Norm1st', true, @islogical);
addParamValue(p, 'NormMeth', 'time', @isstr);
addParamValue(p, 'NormStates', {}, @iscell);
addParamValue(p, 'NormAvg', 'median', @isstr);
addParamValue(p, 'DistanceWeight', 'squaredinverse', @isstr);

parse(p, drf, ts, mf, varargin{:});

%% Load pre-calculated features
clear T
for fidx = 1:length(ts)
    filename = [drf filesep ts{fidx} p.Results.postfn];
    fprintf('Loading feature set from %s.\n', filename)
    tmp = load(filename);
    fn = canonize_fieldname(ts{fidx});
    T.(fn) = tmp.T;
    epdur.(fn) = tmp.epdur;
end
exps = fieldnames(T);
fprintf('\n')
clear filename tmp fn

%% Preprocess features

% This step is necessary to allow for the generation of a uniform set of
% features across animals. It includes removal of bad epochs as well as
% normalization to some standard signal section (like wake), and taking the
% logarithm of the signal to generate a more normal distribution of feature
% values. It will also merge states specified in p.Results.CombineStates.

[T, rmidx, ~, normfact] = preprocess_features(T, 'IsTrainingSet', true, 'EpDur', epdur, ...
    'RemoveStates', p.Results.RemoveStates, 'CombineStates', p.Results.CombineStates, ...
    'PerLen', p.Results.PerLen, 'PerStart', p.Results.PerStart, ...
    'ArtifMap', p.Results.ArtifMap, 'AFMapUse', p.Results.AFMapUse, ...
    'DoNormalize', p.Results.DoNormalize, ...
    'NormFact', p.Results.NormFact, 'NormMap', p.Results.NormMap, 'Norm1st', p.Results.Norm1st, ...
    'NormMeth', p.Results.NormMeth, 'NormStates',p.Results.NormStates, 'NormAvg', p.Results.NormAvg);
for expidx = 1:length(exps)
    T.(exps{expidx}) = T.(exps{expidx})(~rmidx.(exps{expidx}),:);
end
clear expidx rmidx

%% Generate unified training set

% The goal of this function is to use data from all possible source
% (animals) by merging epochs into one unified training set. This way the
% model can generalize better.

fprintf('Generating training set...')
tmpc = cellfun(@(fn) T.(fn), exps, 'UniformOutput', false);
train_set = vertcat(tmpc{:});
fprintf(' done\n\n')
clear tmpc

%% Create cost matrix

% This is a possible way to take care of imbalanced data. Or even in case
% of balanced data, if accurate classification of one state is more
% important than classification of other states, increasing penalty of
% misclassification can be useful.

fprintf('Cost matrix generation...')
clear cost_m
tmp.states = unique(train_set{:,'Scores'});
tmp.nostates = length(tmp.states);
tmp.REMloc = find(strcmp(tmp.states, 'R'));
if strcmpi(p.Results.cost_m, 'offdiag')
    cost_m = ones(tmp.nostates) - eye(tmp.nostates);
elseif strcmpi(p.Results.cost_m, 'rempref')
    if ~isempty(tmp.REMloc)
        cost_m = ones(tmp.nostates) - eye(tmp.nostates);
        cost_m(tmp.REMloc, :) = p.Results.REMboost;
        cost_m(:, tmp.REMloc) = p.Results.REMboost;
        cost_m(tmp.REMloc, tmp.REMloc) = 0;
    else
        fprintf('\nREM was not seen. Using OffDiag.\n');
        cost_m = ones(tmp.nostates) - eye(tmp.nostates);
    end
else
    fprintf('\nUnknown cost matrix name %s. Using OffDiag.\n', ...
        p.Results.cost_m);
    cost_m = ones(tmp.nostates) - eye(tmp.nostates);
end
fprintf(' done\n\n')
clear  tmp

%% Select features for the grouped training set
switch p.Results.DoPCA
    case 'all'
        fprintf('Calculating PCA for unified...\n')
        %In this version all data is used to calcualte PCA. Next version works
        %better because in that case zscoring is done on animals that
        %apparently helps positioning point clouds into same location in PCA
        %space
        [coef, pcs, ~, ~, explained] = pca(train_set{:,2:end}, 'VariableWeights', 'variance');
        ortocoef = diag(std(train_set{:, 2:end}))\coef;
        % To check if enough variance is explained
        if p.Results.PlotExplained
            figure
            pareto(explained)
            xlabel('Principal Component')
            ylabel('Variance Explained (%)')
            grid on
            vline(p.Results.NPCs, 'r')
        end
    case 'first'
        fprintf('Calculating PCA for first...\n')
        %In this version the first animal (quite ad hoc) data is used to
        %calculate PCA projection, then all animals are projected and PCs
        %concatenated into training set
        [coef, ~, ~, ~, explained] = pca(T.(exps{1}){:,2:end}, 'VariableWeights', 'variance');
        ortocoef = diag(std(T.(exps{1}){:, 2:end}))\coef;
        pcs = [];
        for expidx = 1:length(exps)
            pcs = [pcs; zscore(T.(exps{expidx}){:,2:end})*ortocoef]; %#ok<AGROW>
        end
        % To check if enough variance is explained
        if p.Results.PlotExplained
            figure
            pareto(explained)
            xlabel('Principal Component')
            ylabel('Variance Explained (%)')
            grid on
            vline(p.Results.NPCs, 'r')
        end
    case 'each'
        fprintf('Calculating PCA for each')
        pcs = [];
        coef = [];
        divider = ones(size(T.(exps{1}),2)-1, 1);
        for expidx = 1:length(exps)
            fprintf('.')
            [tmpcoef, tmp, ~, ~, ~] = pca(T.(exps{expidx}){:,2:end}, 'VariableWeights', 'variance');
            pcs = [pcs; tmp]; %#ok<AGROW>
            coef = [coef tmpcoef(:, 1:p.Results.NPCs) divider]; %#ok<AGROW>
        end
        fprintf('\n')
        if p.Results.PlotExplained
            disp(coef)
        end
    case 'no'
        
        % If dimensionality of the feature space is not decreased by a PCA
        % of some sort, features are selected by either comparing them to
        % each other (estimate_feature_goodness) or by trying the model 
        % performance on a subsample of epocs (select_features). The below
        % combines these two in different ways.
        
        fprintf('Starting feature selection...\n')
        clear seld_features
        if isempty(p.Results.FeatFile)
            if isempty(p.Results.UseFeat)
                if strcmpi(p.Results.FSKeep, 'All')
                    seld_features = true(1,width(train_set)-1);
                elseif strcmpi(p.Results.FSKeep, 'Suggested')
                    fprintf('estimate_feature_goodness selected:\n')
                    tmp.only = train_set;
                    [~, best_features] = estimate_feature_goodness(tmp);
                    disp(best_features.only)
                    seld_features = select_features(...
                        train_set, cost_m,...
                        'NNeighbors', p.Results.NN, 'Verbose', p.Results.Verbose, 'Summary', true, ...
                        'SfSDir', p.Results.FSDir, 'KeepIn', best_features.only, ...
                        'KeepOut', p.Results.OmitFeat, 'ValMeth', p.Results.ValMeth);
                elseif strcmpi(p.Results.FSKeep, 'Free')
                    seld_features = select_features(...
                        train_set, cost_m, ...
                        'NNeighbors', p.Results.NN, 'Verbose', p.Results.Verbose, 'Summary', true, ...
                        'SfSDir', p.Results.FSDir, 'ValMeth', p.Results.ValMeth);
                else
                    seld_features = select_features(train_set, cost_m, ...
                        'NNeighbors', p.Results.NN, 'Verbose', p.Results.Verbose, 'Summary', true, ...
                        'SfSDir', p.Results.FSDir, 'KeepIn', p.Results.FSKeep, ...
                        'KeepOut', p.Results.OmitFeat, 'ValMeth', p.Results.ValMeth);
                end
            else
                seld_features = ...
                    ismember(train_set.Properties.VariableNames(2:end), ...
                    p.Results.UseFeat);
                if sum(seld_features) ~= length(p.Results.UseFeat)
                    fprintf('Only features: ')
                    sfs = train_set.Properties.VariableNames([false seld_features]);
                    fprintf('%s ', sfs{:})
                    fprintf(' were found from the UseFeat argument.\n')
                else
                    fprintf('Using UseFeat to select ')
                    fprintf('%s ', p.Results.UseFeat{:})
                    fprintf('features for the model.\n')
                end
            end
        else
            fprintf('Loading file: %s\n', p.Results.FeatFile);
            load(p.Results.FeatFile)
        end
        disp(seld_features)
        fprintf('done.\n\n')
        clear trainidx best_features
    otherwise
        error('train_one_model:unknown DoPCA mode')
end

%% Train model

% In case of a kNN training is essentially labelling points in feature
% space using the manual labels.

fprintf('Starting model training...')
if strcmp(p.Results.DoPCA, 'no')
    mdl.only = train_classifier(train_set(:, [true seld_features]), ...
        cost_m, 'NNeighbors', p.Results.NN, 'DistanceWeight', p.Results.DistanceWeight);
else
    fprintf(' (using %i PCs)', p.Results.NPCs);
    varnames = {};
    for vnidx = 1:p.Results.NPCs
        varnames = [varnames {['PC' num2str(vnidx)]}]; %#ok<AGROW>
    end
    mdl.only = train_classifier(...
        [train_set(:, 'Scores') array2table(pcs(:, 1:p.Results.NPCs), 'VariableNames', varnames)], ...
        cost_m, 'NNeighbors', p.Results.NN, 'DistanceWeight', p.Results.DistanceWeight);
end
fprintf(' done.\n\n')
clear train_set cost_m

%% Save model and feature selection
if ~isempty(mf)
    switch p.Results.DoPCA
        case 'all'
            fprintf('Saving PCA-based model and associated orthonormal coefficients in %s.\n', mf);
            PCAMeth = 'all'; %#ok<NASGU>
            save(mf, 'mdl', 'ortocoef', 'PCAMeth')
            seld_features = ortocoef; %Because this is to be returned.
        case 'first'
            fprintf('Saving PCA-based model and associated orthonormal coefficients in %s.\n', mf);
            PCAMeth = 'first'; %#ok<NASGU>
            save(mf, 'mdl', 'ortocoef', 'PCAMeth')
            seld_features = ortocoef; %Because this is to be returned.
        case 'each'
            fprintf('Saving PCA-based model in %s.\n', mf);
            PCAMeth = 'each'; %#ok<NASGU>
            save(mf, 'mdl', 'PCAMeth')
            seld_features = NaN; %Because this is to be returned.
        case 'no'
            fprintf('Saving model and associated feature selection in %s.\n', mf);
            save(mf, 'mdl', 'seld_features', 'normfact')    
    end
end

end

