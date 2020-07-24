function rmidx = flag_outliers(ift, varargin)

%This function finds outliers in a feature table using median and mad.
%
%Usage:
%  FlaggedEpochs = flag_outliers(FeatureTable, ...)
%
%Input arguments
%  FeatureTable: table, contains features and optionally manual scores as
%    returned by generate_statespace().
%
%Output variable
%  FlaggedEpochs: logical, where true the epoch was flagged as outlier
%
%Optional input arguments
%  Method: can be 1 or 2. The first method uses a threshold ('MADAlpha') to
%    detect epochs where median of 'OutlierVariables' is over the
%    threshold. Method 2 calculates a histogram with 'HistBinNo' bins, then
%    identifies the location where 'OutlierVariable' goes below the value
%    of 'HistThresh' and values above this are flagged as outliers. Default
%    is 1.
%  'MADAlpha': 1xN doubles, where N is numel('OutlierVariables'). Default
%    is [15 25].
%  'OutlierVariables': cell of strings, names of features used to detect
%    outliers. Default is {'emg_RMS' 'eeg_RMS'}.
%  'iftName': string, name of data. If not empty, plots information about
%    outlier detection. Default is empty.
%  'HistBinNo': double, see above, default is 500.
%  'HistThresh': double, see above, default is 14.
%  'IsTrainingSet': boolean, if true feature table contains manual labels
%    in its first column. Default is false.
%
%Author: Tamas Kiss <kiss.t@wigner.hu>

%% Parameters
p = inputParser;
addRequired(p, 'ift', @istable)
addParamValue(p, 'Method', 1, @isnumeric); %#ok<*NVREPL>
addParamValue(p, 'MADAlpha', [15 25], @isnumeric);
addParamValue(p, 'OutlierVariables', {'emg_RMS' 'eeg_RMS'}, @iscell);
addParamValue(p, 'iftName', '', @isstr);
addParamValue(p, 'HistBinNo', 500, @isnumeric);
addParamValue(p, 'HistThresh', 14, @isnumeric);
addParamValue(p, 'IsTrainingSet', false, @islogical);
parse(p, ift, varargin{:});

ov = p.Results.OutlierVariables;

%% Main
ifarm = table2array(ift(:, ov));

switch p.Results.Method
    %This is the standard, classical way in which if the value of a
    %variable is over the median by MADAlpha then the epoch with that value
    %is treated as outlier
    case 1
        rmidx = sum(bsxfun(@gt, bsxfun(@rdivide, abs(bsxfun(@minus, ifarm, median(ifarm))), mad(ifarm,1)), p.Results.MADAlpha), 2) > 0;
    %This is another take on the same topic: a histogram is calculated and
    %values of the variable in histogram bins over HistThresh are taken as
    %outliers. This does not have any assumption of normality or anything.
    case 2
        rmidx = false(size(ifarm,1),1);
        thr = NaN(length(ov), 1);
        for ovidx = 1:length(ov)
            [n, x] = hist(ifarm(:, ovidx), p.Results.HistBinNo);
            maxloc = find(n == max(n), 1, 'first');
            thr(ovidx) = x(maxloc-1+find(n(maxloc:end) < p.Results.HistThresh, 1, 'first'));
            rmidx = rmidx | ifarm(:, ovidx) > thr(ovidx);
        end
end

%% Plot if needed
if ~isempty(p.Results.iftName)
    cols = ['r' 'g' 'b' 'm' 'c' 'y'];
    figure('WindowStyle', 'Docked')
    hold on
    if p.Results.Method == 1
        thr = NaN(length(ov), 1);
    end
    myl = NaN(length(ov), 1);
    for pidx = 1:length(ov)
        %if p.Results.Method == 1
            %if length(p.Results.MADAlpha) == 1
            %    ma = p.Results.MADAlpha;
            %else
                ma = p.Results.MADAlpha(pidx);
            %end
        %end
        plot(ift{:, ov{pidx}}, cols(pidx))
        axis tight
        %For Method #1
        if p.Results.Method == 1
            thr(length(ov)+pidx) = median(ift{:, ov{pidx}})+mad(ift{:, ov{pidx}}, 1)*ma;
            hline(thr(length(ov)+pidx), [cols(pidx) '--'])
        else
        %For Method #2
            hline(thr(pidx), cols(pidx))
        end
        myl(pidx) = min(ift{:, ov{pidx}});     
    end
    ylim([min(myl) 1.5*max(thr)])
    title(p.Results.iftName, 'Interpreter', 'None')
    legend(ov, 'Interpreter', 'None')
end

end
