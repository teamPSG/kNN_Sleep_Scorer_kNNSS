function descr = evaluate_model_goodness(labels, orig_labels, varargin)

%Usage:
%  PredictionGoodness = evaluate_model_goodness(Labels, OriginalLabels, ...)
%
%Description:
%  This function calculates confusion matrices for a set of models (the
%  ModelSet) and using this displays a set of measures to describe how well
%  the models do on the population level.
%
%Input variables:
%  Labels: structure, for each of the k experiments an N_kx1 cell of 
%    strings is given for the N epoch labels (scores).
%  OriginalLabels: structure, same as above. Field names of the structures
%    must match, also the number of corresponding N_ks to allow for
%    matching predicted and manual labels.
%
%Output variable:
%  PredictionGoodness: struct, contains the following fields for all model-test
%    data pair:
%    - ConfMatr: table, the confusion matrix C, where C(i, j) counts how
%      many times a test observation (state) "i" was classified as state
%      "j" by the model
%    - PerfMetr: table, contains a number of performance metrics for the
%      states the model was trained on. These metrics are described in the
%      help of confusion2PerformanceMetrics
%    - TrainTab: table, frequency table for the states in the training set
%    - TestTab: table, frequency table for the states in the test set
%    - RenTestTab: table, frequency table for states in the test set if
%      LabelRename was given. Confusion matrix and performance measures are
%      calculated on the renamed test set, with predicted scores also
%      renamed.
%    If the ROC figure was plotted Description also contains a field called
%    figha with the figure handle to it.
%  PredictedScore: a cell array containing the scores as predicted by the
%    classifier
%  RenamedTestLabels: cell array with renamed labels
%
%Optional Input variables:
%  'PlotROC': false/cell of strings: plots the ROC points for the model set
%    for states specified in the value string. No graphical output is
%    produced if set to false. Default is false.
%  'ROCLayout': 1x2 integer, how to organize subplots for ROC. If N is the
%    number of states to be plotted N must be equal to prod(ROCLayout).
%    Default is [2, ceil(length(N)/2)].
%  'ROCLabels': boolean, if set to true individual experiment names label
%    plot, no labels if false. Default is true
%  'PlotAvg': boolean, mean and standard deviation is shown on plot if set
%    to true. Default is true.
%  'ROCNameMap': Nx2 cell of strings, ROCLabels are taken from field names
%    of Label structure. 'ROCNameMap' renames labels. Default is {}, no
%    change in label.
%
%Dependencies:
%  This function calls kappa.m. Ref: Cardillo G. (2007) Cohen.s kappa:
%  compute the Cohen's kappa ratio on a 2x2 matrix. (well, it works for NxN
%  matrices)
%
%See also confusion2PerformanceMetrics, confusionmat, ClassificationKNN/predict
%
%Author: Tamas Kiss <kiss.t@wigner.mta.hu>

%% Parse input and set default parameters
p = inputParser;
addRequired(p, 'labels', @isstruct);
addRequired(p, 'orig_labels', @isstruct);
addParamValue(p, 'PlotROC', false, @iscellstr); %#ok<*NVREPL>
addParamValue(p, 'ROCLayout', [], @isnumeric);
addParamValue(p, 'ROCLabels', true, @islogical);
addParamValue(p, 'PlotAvg', true, @islogical);
addParamValue(p, 'PlotIndiv', true, @islogical);
addParamValue(p, 'ROCNameMap', {}, @iscellstr);
parse(p, labels, orig_labels, varargin{:});

%% Calculate confusion matrices and different performance measures
if length(fieldnames(labels)) ~= length(fieldnames(orig_labels))
    error('evaluate_model_goodness:: Error: Label and Original Label are of different length.')
end
if ~all(strcmp(fieldnames(labels), fieldnames(orig_labels)))
    error('evaluate_model_goodness:: Error: Label and Original Label contain different experiments.')
end
exps = fieldnames(labels);
allstates = {};
for expidx = 1:numel(exps)
    [confm, grpOrder] = confusionmat(orig_labels.(exps{expidx}), labels.(exps{expidx}));
    [sstate, sorder] = sort(grpOrder);
    % Get rid of crazy states by hand
    if ~isempty(find(strcmp(sstate, '-'),1))
        fprintf('Replacing state ''-'' with ''X1''\n')
        sstate{strcmp(sstate, '-')} = 'X1';
    end
    allstates = union(allstates, sstate);
    confm = confm(sorder, sorder);
    descr.(exps{expidx}).PerfMetr = confusion2PerformanceMetrics(confm, sstate, false);
    descr.(exps{expidx}).ConfMatr = array2table(confm, 'VariableNames', sstate, 'RowNames', sstate);
end

%% Generate Name list for ROC plot
RN = exps;
if ~isempty(p.Results.ROCNameMap)
    for expidx = 1:numel(exps)
        RN{expidx} = p.Results.ROCNameMap{strcmp(p.Results.ROCNameMap(:,1), exps{expidx}), 2};
    end
end

%% Plot ROC curves for required states
if iscell(p.Results.PlotROC)
    plotstates = ismember(p.Results.PlotROC, allstates);
    if ~isempty(plotstates) && sum(plotstates) > 0
        plotstatesi = find(plotstates);
        if isempty(p.Results.ROCLayout)
            suplxno = ceil(length(plotstatesi)/2);
            suplyno = 2;
        else
            if prod(p.Results.ROCLayout) < length(plotstatesi)
                fprintf('evaluate_model_goodness:: Error: number of states to plot is %i, but ROCLayout is [%i, %i].\nROC not plotted.\n', ...
                    length(plotstatesi), p.Results.ROCLayout(1), p.Results.ROCLayout(2))
                return
            end
            suplxno = p.Results.ROCLayout(2);
            suplyno = p.Results.ROCLayout(1);
        end
        numFields = numel(length(fieldnames(descr)));
        figha = figure;
        for stateidx = 1:length(plotstatesi)
            subplot(suplyno, suplxno, stateidx)
            plot([0 1], [0 1], 'r--')
            hold on
            %This is to check if the current state exists in all experiments
            if all(struct2array(structfun(@(x) sum(strcmp(x.PerfMetr.Properties.RowNames, ...
                    p.Results.PlotROC{plotstatesi(stateidx)})), descr, 'UniformOutput', false)))
                if p.Results.PlotIndiv
                    plot(structfun(@(x) x.PerfMetr{p.Results.PlotROC{plotstatesi(stateidx)}, 'FPR'}, descr), ...
                        structfun(@(x) x.PerfMetr{p.Results.PlotROC{plotstatesi(stateidx)}, 'SEN'}, descr), 'b*')
                end
                if p.Results.ROCLabels
                    text(structfun(@(x) x.PerfMetr{p.Results.PlotROC{plotstatesi(stateidx)}, 'FPR'}, descr), ...
                        structfun(@(x) x.PerfMetr{p.Results.PlotROC{plotstatesi(stateidx)}, 'SEN'}, descr), RN, 'Interpreter', 'None')
                end
                if p.Results.PlotAvg
                    fprAvg = nanmean(structfun(@(x) x.PerfMetr{p.Results.PlotROC{plotstatesi(stateidx)}, 'FPR'}, descr));
                    fprSem = nanstd(structfun(@(x) x.PerfMetr{p.Results.PlotROC{plotstatesi(stateidx)}, 'FPR'}, descr))/sqrt(numFields);
                    senAvg = nanmean(structfun(@(x) x.PerfMetr{p.Results.PlotROC{plotstatesi(stateidx)}, 'SEN'}, descr));
                    senSem = nanstd(structfun(@(x) x.PerfMetr{p.Results.PlotROC{plotstatesi(stateidx)}, 'SEN'}, descr))/sqrt(numFields);
                    plot([fprAvg-fprSem, fprAvg+fprSem], [senAvg, senAvg], 'r')
                    plot([fprAvg, fprAvg], [senAvg-senSem, senAvg+senSem], 'r')
                end
            else
                for expidx = 1:numel(exps)
                    if sum(strcmp(descr.(exps{expidx}).PerfMetr.Properties.RowNames, ...
                            p.Results.PlotROC{plotstatesi(stateidx)})) > 0
                        sen(expidx) = descr.(exps{expidx}).PerfMetr{p.Results.PlotROC{plotstatesi(stateidx)}, 'SEN'};  %#ok<AGROW>
                        fpr(expidx) = descr.(exps{expidx}).PerfMetr{p.Results.PlotROC{plotstatesi(stateidx)}, 'FPR'};  %#ok<AGROW>
                        if p.Results.PlotIndiv
                            plot(descr.(exps{expidx}).PerfMetr{p.Results.PlotROC{plotstatesi(stateidx)}, 'FPR'}, ...
                                descr.(exps{expidx}).PerfMetr{p.Results.PlotROC{plotstatesi(stateidx)}, 'SEN'}, 'b*')
                        end
                        if p.Results.ROCLabels
                            text(descr.(exps{expidx}).PerfMetr{p.Results.PlotROC{plotstatesi(stateidx)}, 'FPR'}, ...
                                descr.(exps{expidx}).PerfMetr{p.Results.PlotROC{plotstatesi(stateidx)}, 'SEN'}, RN{expidx}, 'Interpreter', 'None')
                        end
                    else
                        %This state does not exist for this experiment
                        sen(expidx) = NaN; %#ok<AGROW>
                        fpr(expidx) = NaN; %#ok<AGROW>
                    end
                end
                if p.Results.PlotAvg
                    fprAvg = nanmean(fpr);
                    fprSem = nanstd(fpr)/sqrt(length(fpr));
                    senAvg = nanmean(sen);
                    senSem = nanstd(sen)/sqrt(length(sen));
                    plot([fprAvg-fprSem, fprAvg+fprSem], [senAvg, senAvg], 'r')
                    plot([fprAvg, fprAvg], [senAvg-senSem, senAvg+senSem], 'r')
                end
            end
            title(['ROC for state ' p.Results.PlotROC{plotstatesi(stateidx)}])
            xlabel('False Positive Rate')
            ylabel('True Positive Rate')
            grid on
            xlim([0 1])
            ylim([0 1])
        end
        descr.figha = figha;
    end
end

end
