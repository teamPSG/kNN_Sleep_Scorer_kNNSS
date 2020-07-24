function [T, rmidx, oldlabels, normfact] = preprocess_features(T, varargin)

% This function performs a pre-processing of features enabling fitting of
% models and predicting labels. Pre-processing steps include removal of
% data with artifacts, missing data segments, removal of manually selected
% epochs, normalization, and log transformation.
%
%Usage:
%  [PreProFeat, RemoveIdx, OldLabels, NormalizingFactors] = ...
%             preprocess_features(Features, ...)
%
%Input arguments:
%  Features: struct of tabels with features. This can be a feature table
%    saved to file.
%
%Outputs:
%  PreProFeat: a structure of faeture tables, same as the input Features
%    follwing preprocessing
%  RemoveIdx: a structure of booleans with TRUE values at epochs that were
%    removed being outliers.
%  OldLabels: struct of cell of strings containing manually assigned labels
%    from Ffeatures.
%  NormalizingFactors: struct of double vectors, the normalization
%    coefficients for features.
%
%Optional input arguments:
%  'RemoveStates': cell of strings. Epochs labelled with labels specified
%    here will be removed from the training set. Default is {'D' 'SS' 'U'}
%  'misshistbinno': double, a helper variable defining the number of bins
%    for a histogram used for detecting epochs with missing signal. Default
%    is 100.
%  'MisSigFfeat': string, used to specify which signal channel to look at
%    to find epochs with missing signal in them. Default is 'eeg_RMS'.
%  'CombineStates', Nx2 cell of strings, specifies what states to combine.
%    First column is the label that will be finally assigned to labels
%    originally labelled by labels in the second column. Default is: 
%    {'W' 'WA'; 'NR' 'NA'; 'R' 'RA'}.
%  'IsTrainingSet': boolean, if true feature table contains manual labels
%    in its first column. Default is false.
%  'EpDur': struct of doubles, duration of an epoch for each experiment
%    given in seconds. Default is [], in which case 10 seconds is used with
%    a warning.
%  'PerLen': double, the length of the period used for calculating feature
%    normalization factors. Given in minutes. Default is 30.
%  'PerStart': double, start time of the period used for normalization,
%    given in minutes. Default is the beginning of the experiment: 0.
%  'WriteFlag': boolean, if set to true 'RemoveFlag' is used to replace
%    original manual label in feature table. Default is true.
%  'RemoveFlag': string, specifies label to be used in feature table to
%    mark epochs removed as outliers. Default is 'Removed'.
%  'DoNormalize': boolean, turns on/off normalization to selected period.
%    Default is TRUE.
%  'DoNormalize': boolean, turns on/off normalization to selected period.
%    Default is TRUE.
%  'DoOutlierDet': boolean, turns on/off outlier detection. Default is
%    true, outliers are detected
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
%  'NormSaveFile': string, if not empty, specifies name of file to store
%    normalization coefficients. Default is empty, do not save
%    coefficients.
%  'ArtifMap': structure, containing a boolean for each experiment. Epochs
%    with TRUE ArtifMap values will be discarded from training. Default is
%    empty structure.
%  'AFMapUse': string, specifies how to combine artifacts detected by 
%    different methods. 'union': all epochs marked as artifact by automated
%    outlier detection and manual input will be marked as artifact;
%    'intersect': only epochs marked as artifact by both methods will be
%    marked as artifact; 'alone': only manually identified artifactual
%    epochs are marked as artifacts. Default is 'union'.
%
% Author: Tamas Kiss <kiss.t@wigner.hu>

%% Parameters
p = inputParser;

addRequired(p, 'T', @isstruct);
addParamValue(p, 'RemoveStates', {'D' 'SS' 'U'}, @iscell); %#ok<*NVREPL>
addParamValue(p, 'misshistbinno', 100, @isnumeric);
addParamValue(p, 'MisSigFeat', 'eeg_RMS', @isstr); %This feature is used to determine missing signal locations
addParamValue(p, 'CombineStates', {...
    'W' 'WA';...
    'NR' 'NA';...
    'R' 'RA'}, @iscell);
addParamValue(p, 'IsTrainingSet', false, @islogical);
addParamValue(p, 'EpDur', [], @isstruct);
addParamValue(p, 'PerLen', 30, @isnumeric); %This is given in minutes
addParamValue(p, 'PerStart', 10, @isnumeric); %Also in minutes
addParamValue(p, 'WriteFlag', true, @islogical);
addParamValue(p, 'DoNormalize', true, @islogical);
addParamValue(p, 'DoLogaritmize', true, @islogical);
addParamValue(p, 'DoOutlierDet', true, @islogical);
addParamValue(p, 'RemoveFlag', 'Removed', @isstr);
addParamValue(p, 'NormFact', [], @isstruct);
addParamValue(p, 'NormMap', [], @iscell);
addParamValue(p, 'Norm1st', true, @islogical);
addParamValue(p, 'NormMeth', 'time', @isstr);
addParamValue(p, 'NormAvg', 'median', @isstr);
addParamValue(p, 'NormStates', {}, @iscell);
addParamValue(p, 'NormSaveFile', '', @isstr);
addParamValue(p, 'ArtifMap', [], @isstruct);
addParamValue(p, 'AFMapUse', 'union', @isstr); %union, intersect, alone

parse(p, T, varargin{:});

exps = fieldnames(T);

if isempty(p.Results.EpDur)
    fprintf('Warning::preprocess_features: Epoch duration not given. Using default 10 sec.\n')
    for expidx = 1:length(exps)
        epdur.(exps{expidx}) = 10;
    end
else
    epdur = p.Results.EpDur;
end

if ~isempty(p.Results.NormMap)
    for expidx = 1:size(p.Results.NormMap,1)
        NormMap{expidx,1} = canonize_fieldname(p.Results.NormMap{expidx, 1}); %#ok<AGROW>
        NormMap{expidx,2} = canonize_fieldname(p.Results.NormMap{expidx, 2}); %#ok<AGROW>
    end
    clear expidx
else
    if ~isempty(p.Results.NormFact)
        if length(exps) ~= length(fieldnames(p.Results.NormFact))
            error('preprocess_features:: cannot generate NormMap automatically.')
        else
            fprintf('Generating NormMap...\n')
            NormMap = [exps fieldnames(p.Results.NormFact)];
        end
    end
end

fprintf('\n')

%% Flag designated states

% There might be some state labels that we don't want to train, like
% artifactual states, or states when animals were dosed or when the
% experiment was disturbed in any ways. If corresponding epochs are marked
% they can be easily removed here.

if isempty(p.Results.ArtifMap) && ~strcmp(p.Results.AFMapUse, 'alone')
    if p.Results.IsTrainingSet
        if ~isempty(p.Results.RemoveStates)
            stxt = '';
            for sidx = 1:length(p.Results.RemoveStates)
                stxt = [stxt ' ' p.Results.RemoveStates{sidx}]; %#ok<AGROW>
            end
            fprintf('Flagging states:%s for removal...', stxt)
        else
            fprintf('No state type flagged for removal.\n\n')
        end
        
        for expidx = 1:length(exps)
            if ~isempty(p.Results.RemoveStates)
                rmidx.(exps{expidx}) = ismember(T.(exps{expidx}){:, 'Scores'}, p.Results.RemoveStates);
            else
                rmidx.(exps{expidx}) = false(height(T.(exps{expidx})), 1);
            end
        end
        
        if ~isempty(p.Results.RemoveStates), fprintf(' done.\n\n'); end
        
        clear stxt sidx expidx
    else
        for expidx = 1:length(exps)
            rmidx.(exps{expidx}) = false(height(T.(exps{expidx})), 1);
        end
    end
end

%% Flag outliers

% This removes epochs where the signal went crazy. Features to watch are
% set in the p.Results.OutlierVariables parameter, defaulting to eeg_RMS
% and emg_RMS. Ie. if there is too much (measured by deviation from mean)
% EEG or EMG power in a given epoch, it will be removed.

if p.Results.DoOutlierDet && isempty(p.Results.ArtifMap) && ~strcmp(p.Results.AFMapUse, 'alone')
    fprintf('Flagging outliers...');
    for expidx = 1:length(exps)
        rmidx.(exps{expidx})(flag_outliers(T.(exps{expidx}))) = true;
    end
    clear tmp testidx trainidx
    fprintf(' done.\n\n')
end

%% Flag segments where signal was missing
if isempty(p.Results.ArtifMap) && ~strcmp(p.Results.AFMapUse, 'alone')
    fprintf('Flagging segments w/ missing signal...')
    for expidx = 1:length(exps)
        [n, x] = hist(T.(exps{expidx}){:, p.Results.MisSigFeat}, p.Results.misshistbinno);
        mnl = find(n == max(n), 1, 'first');
        fzl = find(n==0, 1, 'first');
        if fzl < mnl %Zeros should be between the small value and the max value
            rmidx.(exps{expidx})(...
                T.(exps{expidx}){:, p.Results.MisSigFeat} <= x(fzl))...
                = true;
        end
        % The above heuristics did not work in one of CiTox's data, were 0
        % values not separated from others. Since different features behave
        % differently, I check all the below:
        if p.Results.IsTrainingSet
            rmidx.(exps{expidx})(any(isinf(T.(exps{expidx}){:, 2:end}), 2))...
                = true;
            rmidx.(exps{expidx})(any(isnan(T.(exps{expidx}){:, 2:end}), 2))...
                = true;
            rmidx.(exps{expidx})(any(T.(exps{expidx}){:, 2:end} == 0, 2))...
                = true;
        else
            rmidx.(exps{expidx})(any(isinf(T.(exps{expidx}){:, :}), 2))...
                = true;
            rmidx.(exps{expidx})(any(isnan(T.(exps{expidx}){:, :}), 2))...
                = true;
            rmidx.(exps{expidx})(any(T.(exps{expidx}){:, :} == 0, 2))...
                = true;
        end
    end
    fprintf(' done.\n\n')
    clear expidx n x mnl fzl
end

%% Mix manually calculated artifact map with rmidx

% One has the option to manually mark certain epochs as bad by specifying
% its index.

if ~isempty(p.Results.ArtifMap)
    fprintf('Adding artifact masks using ''%s''...', p.Results.AFMapUse)
    switch p.Results.AFMapUse
        case 'alone'
            rmidx = p.Results.ArtifMap;
        case 'union'
            for expidx = 1:length(exps)
                rmidx.(exps{expidx}) = and(...
                    rmidx.(exps{expidx}), p.Results.ArtifMap.(exps{expidx}));
            end
        case 'intersect'
            for expidx = 1:length(exps)
                rmidx.(exps{expidx}) = or(...
                    rmidx.(exps{expidx}), p.Results.ArtifMap.(exps{expidx}));
            end
        otherwise
            fprintf('\npreprocess_features:: Warning: unknown artifact use mode (%s). Using ''alone''.\n',...
                p.Results.AFMapUse)
            rmidx = p.Results.ArtifMap;
    end
    fprintf(' done.\n\n')
end

%% Normalize to W & do log transformation

% Normalization to the mean wake value of features and taking the logarithm
% is called here. Their order could matter, hence the below decision.

if p.Results.Norm1st
    norm2W
    logtrans
else
    logtrans
    norm2W
end

%% Write flag in table
if p.Results.WriteFlag
    fprintf('Removing bad signal...')
    for expidx = 1:length(exps)
        if p.Results.IsTrainingSet
            oldlabels.(exps{expidx}) = T.(exps{expidx}){:, 'Scores'};
            T.(exps{expidx}){rmidx.(exps{expidx}), 'Scores'} = ...
                cellstr(repmat(p.Results.RemoveFlag, sum(rmidx.(exps{expidx})), 1));
        else
            oldlabels.(exps{expidx}) = {};
        end
    end
    fprintf(' done.\n\n')
    clear expidx
else
    oldlabels = NaN;
end

%% Combine states
if p.Results.IsTrainingSet && ~isempty(p.Results.CombineStates)
    fprintf('Replacing state labels:\n')
    for sidx = 1:size(p.Results.CombineStates,1)
        fprintf('\t%s -> %s\n', p.Results.CombineStates{sidx,2}, p.Results.CombineStates{sidx,1})
    end
    fprintf('...')
    for expidx = 1:length(exps)
        [idx, loc] = ismember(table2cell(T.(exps{expidx})(:, 'Scores')), p.Results.CombineStates(:,2));
        T.(exps{expidx})(idx, 'Scores') = p.Results.CombineStates(loc(idx), 1);
    end
    clear expidx idx loc
    fprintf(' done.\n\n')
end

%% Log transformation -- nested function
    function logtrans
        if p.Results.DoLogaritmize
            fprintf('Performing log transformation...')
            for expidxN = 1:length(exps)
                if p.Results.IsTrainingSet
                    T.(exps{expidxN}){:, 2:end} = log(T.(exps{expidxN}){:, 2:end});
                else
                    T.(exps{expidxN}){:, :} = log(T.(exps{expidxN}){:, :});
                end
            end
            fprintf(' done.\n\n')
        end
    end

%% Normalize to wake -- nested function
    function norm2W
        if p.Results.DoNormalize
            if strcmp(p.Results.NormAvg, 'mean')
                navg = @mean;
            elseif strcmp(p.Results.NormAvg, 'median')
                navg = @median;
            else
                fprintf('Unknown method %s to determine first moment.\n', p.Results.NormAvg)
                error('Unknown method.')
            end
            if isempty(p.Results.NormFact)
                if strcmp(p.Results.NormMeth, 'time')
                    fprintf('Normalizing using given time period...')
                    for expidxN = 1:length(exps)
                        en = exps{expidxN};
                        periods = false(height(T.(en)), 1);
                        periods(round(p.Results.PerStart*60/epdur.(en))+1:round((p.Results.PerStart+p.Results.PerLen)*60/epdur.(en))) = true;
                        periods(rmidx.(en)) = false;
                        if p.Results.IsTrainingSet
                            normfact.(en) = navg(T.(en){periods, 2:end}, 1);
                        else
                            normfact.(en) = navg(T.(en){periods, :}, 1);
                        end
                    end
                elseif strcmp(p.Results.NormMeth, 'zscore')
                    fprintf('Normalizing using zscore...')
                    normfact = NaN;
                elseif strcmp(p.Results.NormMeth, 'state')
                    if p.Results.IsTrainingSet
                        fprintf('Normalizing using given states...')
                        for expidxN = 1:length(exps)
                            en = exps{expidxN};
                            periods = false(height(T.(en)), 1);
                            periods(ismember(T.(en){:, 'Scores'}, p.Results.NormStates)) = true;
                            normfact.(en) = navg(T.(en){periods, 2:end}, 1);
                        end
                    else
                        fprintf('Cannot normalize using state labels because dataset has no Scores field. Exiting.\n')
                        error('Not TrainingSet.')
                    end
                else
                    fprintf('%s unknown normalization method. Exiting.\n', par.Results.NormMeth)
                    error('Unknown method.')
                end
            else
                fprintf('Normalizing using given normalizing factors...\n')
                for nmidx = 1:size(NormMap, 1)
                    fprintf('\t%s -> %s\n', NormMap{nmidx,2}, NormMap{nmidx, 1})
                end
                for expidxN = 1:length(exps)
                    en = exps{expidxN};
                    normfact.(en) = p.Results.NormFact.(NormMap{ismember(NormMap(:,1), en), 2});
                end
            end
            if p.Results.IsTrainingSet
                if strcmp(p.Results.NormMeth, 'zscore')
                    for expidxN = 1:length(exps)
                        en = exps{expidxN};
                        T.(en){:, 2:end} = zscore(T.(en){:, 2:end});
                    end
                else
                    for expidxN = 1:length(exps)
                        en = exps{expidxN};
                        T.(en){:, 2:end} = bsxfun(@rdivide, T.(en){:, 2:end}, normfact.(en));
                    end
                end
            else
                if strcmp(p.Results.NormMeth, 'zscore')
                    for expidxN = 1:length(exps)
                        en = exps{expidxN};
                        T.(en){:, :} = zscore(T.(en){:, :});
                    end
                else
                    for expidxN = 1:length(exps)
                        en = exps{expidxN};
                        T.(en){:, :}     = bsxfun(@rdivide, T.(en){:,     :}, normfact.(en));
                    end
                end
            end
            fprintf(' done.\n\n')
            if ~isempty(p.Results.NormSaveFile)
                fprintf('Saving normalization coeffs in %s.\n\n', p.Results.NormSaveFile);
                save(p.Results.NormSaveFile, 'normfact')
            end
        else
            fprintf('Skipping normalization.\n\n');
            normfact = NaN;
        end
    end

end
