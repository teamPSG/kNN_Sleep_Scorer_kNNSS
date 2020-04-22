function [sel_feat, hist] = select_features(input_set, cost_m, varargin)
%Usage:
%  [SelectedFeatures, History] = select_features(InputSet, CostMatrix, ...)
%
%Description:
%  This function selects features that maximize goodness of prediction of
%  categories. It uses resubstitution (same dataset for training and
%  testing) or cross-validation, and forward or backward feature selection.
%
%Input variable:
%  InputSet: table, contains the scores and the metrics calculated from raw
%    data by generate_trainingset(). The first column of InputSet is called
%    'Scores' and contains strings describing the sleep stage. Subsequent
%    columns store the metrics (features) that describe a given time bin.
%    Each row describes a time bin.
%  CostMatrix: double, an NxN matrix describing the cost of
%    misclassification. N is the number of Score types in InputSet.
%
%Optional input arguments:
%  'Verbose': boolean, if true the forward selection procedure displays its
%    history. Default is true.
%  'NNeighbors': integer, the number of neighbors for the k nearest
%    neighbor classifier. Default 3.
%  'SfsDir': string, specifies the direction of sequential feature
%    selection. Possible values are 'forward', and 'backward'. Default is
%    backward.
%  'ValMeth': integer/string, specifies validation method. If an integer
%    value of N is given, N-fold cross validation will be performed. If the
%    string 'resub' is given, resubstitution validation is performed.
%    Default is 5 to perform 5-fold cross validation.
%  'Summary': logical, if true a summary of selected features is printed.
%    Default is true.
%  'KeepIn': cell of strings, the name of features to be kept during
%    feature selection. Default is none.
%  'KeepOut': cell of strings, the name of features to be omitted during
%    feature selection. Default is none.
%
%Output variables:
%  SelectedFeatures: boolean, a vector that has true values at each
%    selected features and false values at unselected features. The order
%    of features is specified in InputSet.
%  History: struct, returns information on which feature is chosen at each
%    step See doc sequentialfs for more help
%
%Dependencies:
%  This function uses the Statistics Toolbox of Matlab
%
%See also sequentialfs, ClassificationKNN
%
%Author: Tamas Kiss <kiss.t@wigner.hu>

%% Parse input parameters
p = inputParser;
addRequired(p, 'input_set', @istable);
addRequired(p, 'cost_m', @isnumeric);
addParamValue(p, 'Verbose', true, @islogical); %#ok<*NVREPL>
addParamValue(p, 'NNeighbors', 3, @(x)validateattributes(x,{'numeric'},{'integer','scalar','>',0}));
addParamValue(p, 'SfsDir', 'backward', @(x) any(validatestring(x, {'forward', 'backward'})));
addParamValue(p, 'Summary', true, @islogical);
addParamValue(p, 'ValMeth', 5, @(x)validateattributes(x,{'char', 'numeric'}, {'nonempty'}));
addParamValue(p, 'KeepIn', {}, @iscellstr);
addParamValue(p, 'KeepOut', {}, @iscellstr);
parse(p, input_set, cost_m, varargin{:});

tmpfeatl = length(unique(input_set{:, 'Scores'}));
if size(p.Results.cost_m, 1) ~= size(p.Results.cost_m, 2) || ...
        size(p.Results.cost_m, 1) ~= tmpfeatl
        fprintf(['train_classifier:: Warning: CostMatrix is of inappropriate size.'...
            'Using ones(%i)-eye(%i) instead.\n'], tmpfeatl, tmpfeatl)
        cost_m = ones(tmpfeatl)-eye(tmpfeatl);
end

if ~isempty(p.Results.KeepIn)
    keepin = ismember(input_set.Properties.VariableNames, p.Results.KeepIn);
    keepin = keepin(2:end);
else
    keepin = {};
end

if ~isempty(p.Results.KeepOut)
    keepout = ismember(input_set.Properties.VariableNames, p.Results.KeepOut);
    keepout = keepout(2:end);
else
    keepout = {};
end

%% Internal parameters
if p.Results.Verbose
    opts = statset('Display','iter');
else
    opts = statset('Display','off');
end

%% Assign numbers to state character representations
states = unique(input_set(:, 'Scores'));
scores = zeros(height(input_set),1);
for sidx = 1:numel(states)
    scores(strcmp(input_set{:, 'Scores'}, states{sidx,1})) = sidx;
end

%% Generate partitioning of data
if isnumeric(p.Results.ValMeth)
    c2 = cvpartition(scores, 'k', p.Results.ValMeth);
    opts = statset(opts, 'UseParallel', true);
elseif strcmp(p.Results.ValMeth, 'resub')
    c2 = cvpartition(numel(scores), 'resubstitution');
    opts = statset(opts, 'UseParallel', false);
else
    fprintf('select_features:: Warning: unknown validation method given for ''ValMeth'' optional argument. Using 5-fold cross validation.\n')
    c2 = cvpartition(scores, 'k', 5);
    opts = statset(opts, 'UseParallel', true);
end

%% Function handle for counting misclassification of objects in a KNN classifier
fun = @(Xtrain,Ytrain,Xtest,Ytest) ...
        sum(Ytest~= predict(ClassificationKNN.fit(Xtrain,Ytrain,'NumNeighbors',p.Results.NNeighbors,'Cost',cost_m,'NSMethod','exhaustive'), Xtest));

%% Feature selection
[sel_feat, hist] = sequentialfs(fun, input_set{:, 2:end}, scores, 'cv', c2,...
    'options',opts, 'direction',p.Results.SfsDir, 'KeepIn', keepin, 'KeepOut', keepout);

%% Some more useful output
if p.Results.Summary
    fprintf('Selected features:\n')
    disp(input_set.Properties.VariableNames([false sel_feat]))
end

end
