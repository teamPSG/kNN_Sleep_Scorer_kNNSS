function mdl = train_classifier(input_set, cost_m, varargin)
%Usage:
%  Model = train_classifier(InputSet, CostMatrix, ...)
%
%Description:
%  This function trains a KNN classifier using.
%
%Input variables:
%  InputSet: table, contains the scores and the metrics calculated from raw
%    data by generate_trainingset(). The first column of InputSet is called
%    'Scores' and contains strings describing the sleep stage. Subsequent
%    columns store the metrics (features) that describe a given time bin.
%    Each row describes a time bin.
%  CostMatrix: double, an NxN matrix describing the cost of
%    misclassification. N is the number of Score types in InputSet.
%
%Optional input variables:
%  'WriteTmpFile': string, name of a temporary file containing the training
%    set that will be written to disk.
%  'OutDir': string, path to output directory, defult is the working
%    directory (.)
%  'NNeighbors': integer, number of neighbors used by the KNN Classifier.
%    Default is 3.
%  'DistanceWeight': string, method used by kNN fitting to take distance of
%    points into account. See Matlab documentation for ClassificationKNN
%    for details. Default is 'squaredinverse'.
%
%Output variable:
%  Model: ClassificationKNN, the model that can be used to classify a test
%    set
%
%Dependencies:
%  This function uses the Statistics Toolbox of Matlab
%
%See also ClassificationKNN
%
%Author: Tamas Kiss <kiss.t@wigner.hu>

%% Parse input and set default parameters
p = inputParser;
addRequired(p, 'input_set', @istable);
addRequired(p, 'cost_m', @isnumeric);
addParamValue(p, 'WriteTmpFile', false, @isstr); %#ok<*NVREPL>
addParamValue(p, 'OutDir', '.', @isstr);
addParamValue(p, 'NNeighbors', 3, @isnumeric);
addParamValue(p, 'DistanceWeight', 'squaredinverse', @isstr);
parse(p, input_set, cost_m, varargin{:});

tmpfeatl = length(unique(input_set{:, 'Scores'}));
if size(p.Results.cost_m, 1) ~= size(p.Results.cost_m, 2) || ...
        size(p.Results.cost_m, 1) ~= tmpfeatl
        fprintf(['train_classifier:: Warning: CostMatrix is of inappropriate size. '...
            'Using ones(%i)-eye(%i) instead.\n'], tmpfeatl, tmpfeatl)
        cost_m = ones(tmpfeatl)-eye(tmpfeatl);
end

%% Train classifier
mdl = ClassificationKNN.fit(input_set{:, 2:end}, input_set{:, 'Scores'}, 'NumNeighbors', p.Results.NNeighbors, ...
        'Cost',cost_m, 'NSMethod','exhaustive', 'PredictorNames', input_set.Properties.VariableNames(:, 2:end), ...
        'ResponseName','Score', 'DistanceWeight', p.Results.DistanceWeight);

%% Write output to disk
if p.Results.WriteTmpFile
    save([p.Results.OutDir filesep p.Results.WriteTmpFile], 'mdl');
end
    
end

