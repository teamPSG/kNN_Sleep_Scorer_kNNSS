function [fgts, bigs] = estimate_feature_goodness(ift, varargin)

%This function measures the overlap between feature distributions of two
%states.
%
%Usage:
%  [FeatureGoodnessTableSet BestFeaturesSet] = ...
%      estimate_feature_goodness(InputFeatureTableSet, ...)
%
%Input arguments:
%  InputFeatureTableSet: struct of tables, contains the scores and the
%    metrics calculated from raw
%    data by generate_trainingset() for each experiment. The first column
%    of tables is called 'Scores' and contains strings describing the sleep
%    stage. Subsequent columns store the metrics (features) that describe a
%    given time bin. Each row describes a time bin.
%  label1, 2: strings, names of the sleep states to be considered for
%    calculating feature distribution overlap
%
%Optional input arguments:
%  'Verbose': logical, if set to true a plot with each feature will be
%    plotted to graphically represent distributions. Default is false.
%  'BigN': integer, specifying how many of the feature names for the
%    biggest goodness values to return. Default is 2.
%  'HistBinNo': integer, number of histogragram bins to be used for
%    calculating feature value distributions. Default is 100.
%  'PlotRatio': 1x2 integer, specifies the ratio of subplots in the
%    horizontal and vertical direction on the screen
%
%Output variables:
%  FeatureGoodnessTableSet: struct of struct of table, for each label pair
%    found in an experiment in ift a table is given with 1 row that
%    specifies the overlap between each given feature probability
%    distribution. If for the two states a certain feature shows the same 
%    distribution, the distributions totally overlap and the goodness value
%    will be 0. In case, however, when the two distributions are completely
%    disjunct the goodness value will be 1.
%  BestFeaturesSet: struct of cellstr, giving the BigN number of features
%    with the highest goodness values. This set is the union of best
%    features for all label pair calculations.
%
%
%Author: Tamas Kiss <kiss.t@wigner.hu>

%% Parse input and set default parameters
p = inputParser;
addRequired(p, 'ift', @isstruct);
addParamValue(p, 'Verbose', false, @islogical); %#ok<*NVREPL>
addParamValue(p, 'HistBinNo', 100, @isnumeric);
addParamValue(p, 'BigN', 2, @isnumeric);
addParamValue(p, 'PlotRatio', [2 3], @isnumeric);
parse(p, ift, varargin{:});

%% Main

% We want to see how much the distribution of any pair of features overlap.
% If they overlap then they contain the same information about states so
% their goodness value will be low. We'd like to use good features, ie.
% those that do not overlap too much.

exps = fieldnames(ift);
for expidx = 1:length(exps)
    fnames = ift.(exps{expidx}).Properties.VariableNames(2:end); %Names of features, 1 being Scores
    bigs.(exps{expidx}) = [];
    if p.Results.Verbose
        fprintf('Processing %s.\n', exps{expidx})
        tmp.ydim = ceil(p.Results.PlotRatio(1) * sqrt(length(fnames)/prod(p.Results.PlotRatio)));
        tmp.xdim = ceil(length(fnames)/tmp.ydim);
    end
    lnames = unique(ift.(exps{expidx}){:, 'Scores'});
    for l1idx = 1:length(lnames)-1
%fprintf('l1idx: %i\n', l1idx);
        l1 = lnames{l1idx};
        for l2idx = l1idx+1:length(lnames)
%fprintf('\tl2idx: %i\n', l2idx);
            l2 = lnames{l2idx};
            tmp.vals = NaN(1, length(fnames));
            if p.Results.Verbose
                tmp.fh = figure('Name', [exps{expidx} ', ' l1 '|' l2]);
                hold on
            end
            for fidx = 1:length(fnames)
%fprintf('%i ', fidx);
                tmp.l1 = ...
                    ift.(exps{expidx}){strcmp(ift.(exps{expidx}){:, 'Scores'}, l1), fnames{fidx}};
                tmp.l2 = ...
                    ift.(exps{expidx}){strcmp(ift.(exps{expidx}){:, 'Scores'}, l2), fnames{fidx}};
                tmp.xmin = min([min(tmp.l1) min(tmp.l2)]);
                tmp.xmax = max([max(tmp.l1) max(tmp.l2)]);
                tmp.xax = linspace(tmp.xmin, tmp.xmax, p.Results.HistBinNo);
                tmp.nl1 = hist(tmp.l1, tmp.xax);
                tmp.nl2 = hist(tmp.l2, tmp.xax);
                tmp.vals(fidx) = sum(abs(tmp.nl1/sum(tmp.nl1)-tmp.nl2/sum(tmp.nl2)))/2;
                if p.Results.Verbose
                    subplot(tmp.ydim, tmp.xdim, fidx)
                    plot(tmp.xax, tmp.nl1/sum(tmp.nl1))
                    hold on
                    plot(tmp.xax, tmp.nl2/sum(tmp.nl2))
                    grid on
                    title([fnames{fidx} ', d: ' num2str(tmp.vals(fidx))], 'Interpreter', 'None')
                    %legend(l1, l2)
                end
            end
%fprintf('\n')
            if p.Results.Verbose
                tightfit;
                set(tmp.fh, 'WindowStyle', 'Docked')
            end
            fgts.(exps{expidx}).(['l' l1 '_' l2]) = array2table(tmp.vals);
            fgts.(exps{expidx}).(['l' l1 '_' l2]).Properties.VariableNames = fnames;
            [~, i] = sort(tmp.vals, 'descend');
            bigs.(exps{expidx}) = union(bigs.(exps{expidx}), ...
                fgts.(exps{expidx}).(['l' l1 '_' l2]).Properties.VariableNames(i(1:p.Results.BigN)));
        end
    end
end
if p.Results.Verbose
    fprintf('Feature suggestions: \n')
    for expidx = 1:length(exps)
        fprintf('%s.\n', exps{expidx})
        disp(bigs.(exps{expidx}))
    end
end

end