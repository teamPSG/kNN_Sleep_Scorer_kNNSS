%% Description
% Since REM and Cataplexy prediction is somewhat crappy, try to boost it
% by removing some W and NR points from the model
%
% Author: Tamas Kiss <kiss.t@wigner.mta.hu>
% Date: 2020-02-19

%% Parameters
% Model realization specific parameters
%ATX3
par.ModelFile = '/home/umat/bognor/Takeda_Temporary/ATX3_50_model.mat';
par.TestExpFilter = '.*';
par.intdir = '/home/umat/bognor/Takeda_Temporary/50partitioned'; %Training and test sets from here
par.TestingSetExp = '*_test';
par.NormMap = {};
par.ROCNameMap = {...
    'No1_6p_11p_test' 'A1'; ...
    'No2_6p_11p_test' 'A2'; ...
    'No3_6p_11p_test' 'A3'; ...
    'No4_6p_11p_test' 'A4'; ...
    'No5_6p_11p_test' 'A5'; ...
    'No6_6p_11p_test' 'A6'; ...
    'No7_6p_11p_test' 'A7'; ...
    'No8_6p_11p_test' 'A8'; ...
    };

%DTA
% par.ModelFile = '/home/umat/bognor/Takeda_Temporary/DTA_50_model.mat';
% par.TestExpFilter = '.*';
% par.intdir = '/home/umat/bognor/Takeda_Temporary/50partitioned'; %Training and test sets from here
% par.TestingSetExp = '*_test';
% par.NormMap = {...
%     'DTA_Tg_No_10_6p_7p_test'  'DTA_Tg_No_10_6p_7p_train'; ...
% 	'DTA_Tg_No_2_6p_7p_test'   'DTA_Tg_No_2_6p_7p_train'; ...
% 	'DTA_Tg_No_6_6p_7p_test'   'DTA_Tg_No_6_6p_7p_train'; ...
% 	'DTA_Tg_No_10_7p_11p_test' 'DTA_Tg_No_10_7p_11p_train'; ...
% 	'DTA_Tg_No_2_7p_11p_test'  'DTA_Tg_No_2_7p_11p_train'; ...
% 	'DTA_Tg_No_6_7p_11p_test'  'DTA_Tg_No_6_7p_11p_train'; ...
% 	'DTA_Tg_No_11_6p_7p_test'  'DTA_Tg_No_11_6p_7p_train'; ...
% 	'DTA_Tg_No_3_6p_7p_test'   'DTA_Tg_No_3_6p_7p_train'; ...
% 	'DTA_Tg_No_7_6p_7p_test'   'DTA_Tg_No_7_6p_7p_train'; ...
% 	'DTA_Tg_No_11_7p_11p_test' 'DTA_Tg_No_11_7p_11p_train'; ...
% 	'DTA_Tg_No_3_7p_11p_test'  'DTA_Tg_No_3_7p_11p_train'; ...
% 	'DTA_Tg_No_7_7p_11p_test'  'DTA_Tg_No_7_7p_11p_train'; ...
% 	'DTA_Tg_No_12_6p_7p_test'  'DTA_Tg_No_12_6p_7p_train'; ...
% 	'DTA_Tg_No_4_6p_7p_test'   'DTA_Tg_No_4_6p_7p_train'; ...
% 	'DTA_Tg_No_8_6p_7p_test'   'DTA_Tg_No_8_6p_7p_train'; ...
% 	'DTA_Tg_No_12_7p_11p_test' 'DTA_Tg_No_12_7p_11p_train'; ...
% 	'DTA_Tg_No_4_7p_11p_test'  'DTA_Tg_No_4_7p_11p_train'; ...
% 	'DTA_Tg_No_8_7p_11p_test'  'DTA_Tg_No_8_7p_11p_train'; ...
% 	'DTA_Tg_No_1_6p_7p_test'   'DTA_Tg_No_1_6p_7p_train'; ...
% 	'DTA_Tg_No_5_6p_7p_test'   'DTA_Tg_No_5_6p_7p_train'; ...
% 	'DTA_Tg_No_9_6p_7p_test'   'DTA_Tg_No_9_6p_7p_train'; ...
% 	'DTA_Tg_No_1_7p_11p_test'  'DTA_Tg_No_1_7p_11p_train'; ...
% 	'DTA_Tg_No_5_7p_11p_test'  'DTA_Tg_No_5_7p_11p_train'; ...
% 	'DTA_Tg_No_9_7p_11p_test'  'DTA_Tg_No_9_7p_11p_train'};
% par.ROCNameMap = {...
%     'DTA_Tg_No_10_6p_7p_test'  'A10_6p-7p'; ...
% 	'DTA_Tg_No_2_6p_7p_test'   'A2_6p-7p'; ...
% 	'DTA_Tg_No_6_6p_7p_test'   'A6_6p-7p'; ...
% 	'DTA_Tg_No_10_7p_11p_test' 'A10_7p-11p'; ...
% 	'DTA_Tg_No_2_7p_11p_test'  'A2_7p-11p'; ...
% 	'DTA_Tg_No_6_7p_11p_test'  'A6_7p-11p'; ...
% 	'DTA_Tg_No_11_6p_7p_test'  'A11_6p-7p'; ...
% 	'DTA_Tg_No_3_6p_7p_test'   'A3_6p-7p'; ...
% 	'DTA_Tg_No_7_6p_7p_test'   'A7_6p-7p'; ...
% 	'DTA_Tg_No_11_7p_11p_test' 'A11_7p-11p'; ...
% 	'DTA_Tg_No_3_7p_11p_test'  'A3_7p-11p'; ...
% 	'DTA_Tg_No_7_7p_11p_test'  'A7_7p-11p'; ...
% 	'DTA_Tg_No_12_6p_7p_test'  'A12_6p-7p'; ...
% 	'DTA_Tg_No_4_6p_7p_test'   'A4_6p-7p'; ...
% 	'DTA_Tg_No_8_6p_7p_test'   'A8_6p-7p'; ...
% 	'DTA_Tg_No_12_7p_11p_test' 'A12_7p-11p'; ...
% 	'DTA_Tg_No_4_7p_11p_test'  'A4_7p-11p'; ...
% 	'DTA_Tg_No_8_7p_11p_test'  'A8_7p-11p'; ...
% 	'DTA_Tg_No_1_6p_7p_test'   'A1_6p-7p'; ...
% 	'DTA_Tg_No_5_6p_7p_test'   'A5_6p-7p'; ...
% 	'DTA_Tg_No_9_6p_7p_test'   'A9_6p-7p'; ...
% 	'DTA_Tg_No_1_7p_11p_test'  'A1_7p-11p'; ...
% 	'DTA_Tg_No_5_7p_11p_test'  'A5_7p-11p'; ...
% 	'DTA_Tg_No_9_7p_11p_test'  'A9_7p-11p'};

%IO
par.functiondir = '/home/umat/bognor/PfSS/Deployment/Functions/'; %This is the folder in which functions will live
par.postfn = '_trset.mat'; %This part will be appended to the name of the EDF file and features will be saved in this file for each experiment

%Testing
par.WakeStart = 0; %minutes
par.WakeLength = 30; %minutes
par.CombineStates = {...
    'W' 'Wake'
    'W' 'Wake ?'
    'NR' 'SW';...
    'R' 'REM';...
    'C' 'C1';...
    'C' 'C2'};
par.fullorder = {'W' 'NR', 'R', 'C'};
par.NormMeth = 'time';
par.NormStates = {'W'};
par.avgf = @nanmean;
par.df = @nanstd;

par.PlotLastROC = true;
par.PlotNormFactComp = false;
par.Movie = false;
par.Plot3D = false;

par.sigma = 0.25;
par.NewREM = 1000;
par.InflateREM = false;
par.ReduceNRW = true;
par.ModelParameters = 0.92*ones(1,50);

%% Set the path
restoredefaultpath;
rmpath('/home/umat/Documents/MATLAB')
clear RESTOREDEFAULTPATH_EXECUTED
addpath(par.functiondir)
rehash

%% Load model
fprintf('Loading base model from %s...', par.ModelFile)
load(par.ModelFile)
fprintf(' done.\n\n')

%% Load pre-calculated features for testing using par.TestingSetExp
warning('OFF', 'MATLAB:table:DimNameNotValidIdentifierBackCompat')
clear T
files = dir([par.intdir filesep par.TestingSetExp par.postfn]);
fns = {files.name};
fns = fns(not(cellfun(@isempty, regexp(fns, par.TestExpFilter))))';
for expidx = 1:length(fns)
    ifn = [par.intdir filesep fns{expidx}];
    fprintf('Loading %s.\n', ifn)
    tmp = load(ifn);
    expname = canonize_fieldname(strrep(files(expidx).name, par.postfn, ''));
    epdur.(expname) = tmp.epdur;
    T.(expname) = tmp.T;
end
expnames = fieldnames(T);
fprintf('\n')
clear expidx tmp expname files ifn
warning('ON', 'MATLAB:table:DimNameNotValidIdentifierBackCompat')

%% Feature preprocessing
[T, rmidx, ~, feat_nfact] = preprocess_features(T, 'IsTrainingSet', true, 'EpDur', epdur, ...
    'PerLen', par.WakeLength, 'PerStart', par.WakeStart, 'WriteFlag', false, ...
    'CombineStates', par.CombineStates, ...
    'NormMeth', par.NormMeth, 'NormStates', par.NormStates, ...
    'NormFact', normfact, 'NormMap', par.NormMap);
clear epdur

%% Model pre-conditioning
fprintf('Preconditioning model.\n\n')
OrigPrior = mdl.only.Prior;
%mdl.only.Prior = OrigPrior .* [1 1.5 1];
mdl.only.Cost = ones(4)-eye(4);
mdl.only.DistanceWeight = 'squaredinverse';
%mdl.only.BreakTies = 'nearest';
%mdl.only.IncludeTies = true;
mdl.only.NumNeighbors = 10;

%% Big for loop through models
pari = 1;
clear sen fpr mcc PT
for modpar = par.ModelParameters
    %% Classifier REM inflation
    % Inflate REM-labelled part of classifier by adding an extra point to
    % each existing REM epoch
    if par.InflateREM
        REMLoc = find(strcmp(mdl.only.Y, 'R'));
        NewREM = round(length(REMLoc)*modpar); %Use this when number of new points is the parameter to vary
        %NewREM = round(length(REMLoc)*par.NewREM); %Use this when number of new points is not the parameter to vary
        if NewREM > 0
            %jitter = repmat(normrnd(0, par.sigma, 1, size(mdl.only.X,2)), NewREM, 1); %Same shift for all observations
            jitter = normrnd(0, par.sigma, NewREM, size(mdl.only.X,2)); %Use this when sitter size is not the parameter to vary
            %jitter = normrnd(0, modpar, NewREM, size(mdl.only.X,2)); %Use this when jitter size is the parameter to vary
            InflIdx = datasample(REMLoc, NewREM);
            X = [mdl.only.X; mdl.only.X(InflIdx, :) + jitter];
            Y = [mdl.only.Y; mdl.only.Y(InflIdx)];
        else
            X = mdl.only.X;
            Y = mdl.only.Y;
        end
        fprintf('Classifier input data expanded by %i new REM points.\n', NewREM)
        clear REMLoc NewREM
    end
    
    %% Classifier NR and W reduction
    % Delete some NR and W points
    if par.ReduceNRW
        Wloc = find(strcmp(mdl.only.Y, 'W'));
        NRloc = find(strcmp(mdl.only.Y, 'NR'));
        DelW = round(length(Wloc)*modpar);
        DelNR = round(length(NRloc)*modpar);
        X = mdl.only.X;
        Y = mdl.only.Y;
        if DelW > 0 && DelNR > 0 && DelW < length(Wloc) && DelNR < length(NRloc)
            DelWIdx = datasample(Wloc, DelW, 'Replace', false);
            DelNRIdx = datasample(NRloc, DelNR, 'Replace', false);
            X([DelWIdx; DelNRIdx], :) = [];
            Y([DelWIdx; DelNRIdx], :) = [];
        end
        fprintf('%i W and %i NR were deleted. Originally %i observations -> %i new observations.\n', ...
            DelW, DelNR, mdl.only.NumObservations, length(Y))
        clear Wloc NRloc DelW DelNR DelWIdx DelNRIdx
    end
    
    %% Retrain classifier if needed
    if par.ReduceNRW || par.InflateREM
        fprintf('Retraining classifier...')
        NewMod = ClassificationKNN.fit(X, Y, 'NumNeighbors', mdl.only.NumNeighbors, ...
            'Cost', mdl.only.Cost, 'NSMethod', 'exhaustive', 'PredictorNames', mdl.only.PredictorNames, ...
            'ResponseName', 'Score', 'DistanceWeight', mdl.only.DistanceWeight);
        fprintf(' done.\n')
    else
        NewMod = mdl.only;
    end
    clear X Y
    
    %% Predict labels
    fprintf('Predicting labels...')
    clear labels
    for expidx = 1:length(expnames)
        labels.(expnames{expidx}) = NewMod.predict(T.(expnames{expidx}){:, [false seld_features]});
    end
    fprintf(' done.\n\n')
    clear expidx
    
    %% Calculate confusion matrix-derived measures
    for expidx = 1:length(expnames)
        [cm, go] = confusionmat(T.(expnames{expidx}){:, 'Scores'}, labels.([expnames{expidx}]), 'order', par.fullorder);
        tmp = confusion2PerformanceMetrics(cm, go, false);
        sen(:, expidx, pari) = tmp{:, 'SEN'}; %#ok<SAGROW>
        fpr(:, expidx, pari) = tmp{:, 'FPR'}; %#ok<SAGROW>
        mcc(expidx, pari) = m_corr_coeff(cm); %#ok<SAGROW>
    end
    if ~iscell(modpar)
        PT{pari} = [num2str(modpar) ': ' num2str(par.avgf(mcc(:, pari)), '%0.2f') '\pm' num2str(par.df(mcc(:, pari)), '%0.2f')]; %#ok<SAGROW>
    else
        PT{pari} = [cell2mat(modpar) ': ' num2str(par.avgf(mcc(:, pari)), '%0.2f') '\pm' num2str(par.df(mcc(:, pari)), '%0.2f', '%0.2f')]; %#ok<SAGROW>
    end
    pari = pari + 1;
    clear expidx cm go tmp mcc
    
end
clear modpar pari

%% Plot ROC
figure
p = panel();
xs = ceil(length(par.fullorder)/2);
ys = 2;
p.pack(ys, xs);
for sidx = 1:length(par.fullorder)
    p(floor((sidx-1)/xs)+1,mod(sidx-1, ys)+1).select();
    plot(squeeze(par.avgf(fpr(sidx, :, :), 2)), squeeze(par.avgf(sen(sidx, :, :), 2)), '*-')
    if sidx == 3
        %tpts = [1 30 40 47 50];
        tpts = [];
    else
        tpts = [];
    end
    text(squeeze(par.avgf(fpr(sidx, :, tpts), 2)), squeeze(par.avgf(sen(sidx, :, tpts), 2)), PT(tpts))
    title(par.fullorder{sidx})
    xlim([0 1])
    ylim([0 1])
    grid on
    xlabel('False Positive Rate')
    if sidx == 1 || sidx == 3
        ylabel('Sensitivity')
    end
    
    %For saving in PNG
    set(gcf, 'PaperPosition', [-4.1132    2.4579   16.7263    6.0842])
end
clear sidx p tpts xs ys

%% Individual ROC on last run
if par.PlotLastROC
    ol = structfun(@(x) x{:, 'Scores'}, T, 'UniformOutput', false);
    evaluate_model_goodness(labels, ol, ...
        'PlotROC', par.fullorder, 'ROCLayout', [2,2], ...
        'ROCLabels', true, 'ROCNameMap', par.ROCNameMap);%, 'PlotIndiv', false);

    %For saving in PNG
    set(gcf, 'PaperPosition', [-4.1132    2.4579   16.7263    6.0842])
    ax = get(gcf, 'Children');
    set(ax(4), 'OuterPosition', [0 0.5 0.5 0.5])
    set(ax(2), 'OuterPosition', [0 0 0.5 0.5])
    set(ax(3), 'OuterPosition', [0.5 0.5 0.5 0.5])
    set(ax(1), 'OuterPosition', [0.5 0 0.5 0.5])
    clear ol ax
end

%% Average confusion matrix for last run
for expidx = 1:length(expnames)
    cm(:, :, expidx) = confusionmat(T.(expnames{expidx}){:, 'Scores'}, labels.([expnames{expidx}]), 'order', par.fullorder); %#ok<SAGROW>
end
fprintf('Average confusion matrix:\n')
disp(array2table(mean(cm,3), 'VariableNames', par.fullorder, 'RowNames', par.fullorder))
fprintf('Std of confusion matrices:\n')
disp(array2table(std(cm,[], 3), 'VariableNames', par.fullorder, 'RowNames', par.fullorder))
%clear expidx cm

%% Compare normalizing factors
if par.PlotNormFactComp
    figure
    hold on
    mfns = fieldnames(normfact);
    ffns = fieldnames(feat_nfact);
    clrs = {'r' 'b' 'g' 'm'};
    for fidx = 1:length(expnames)
        rat = normfact.(mfns{fidx})./feat_nfact.(ffns{fidx});
        h(fidx) = plot(rat, clrs{fidx}); %#ok<SAGROW>
        plot(find(seld_features), rat(seld_features), '*', 'Color', clrs{fidx})
        fprintf('comparing %s with %s.\n', mfns{fidx}, ffns{fidx})
    end
    ylabel('NormFact ratio')
    legend(h, expnames, 'Interpreter', 'None')
    axis tight
    grid on
    hline(1, 'k')
    vline(find(seld_features))
    clear fidx mfns ffns
end

%% Compare original and new models with data to classify
if par.Plot3D
    animark = {'s', 'p', '^', '*', '>', '<', 'h', 'd'};
    [coeff, score, ~, ~, explained, mu] = pca(mdl.only.X);
    [Ncoeff, Nscore, ~, ~, Nexplained, Nmu] = pca(NewMod.X);
    pcs = (NewMod.X-mu)*coeff(:, 1:3); %#ok<NASGU>
    
    mdlREMi = strcmp(mdl.only.Y, 'R');
    mdlNRi = strcmp(mdl.only.Y, 'NR');
    mdlWi = strcmp(mdl.only.Y, 'W');
    mdlCi = strcmp(mdl.only.Y, 'C');
    
    NmdlREMi = strcmp(NewMod.Y, 'R');
    NmdlNRi = strcmp(NewMod.Y, 'NR');
    NmdlWi = strcmp(NewMod.Y, 'W');
    NmdlCi = strcmp(NewMod.Y, 'C');
    
    figure
    hold on
    %plot3(pcs(NmdlNRi,1), pcs(NmdlNRi,2), pcs(NmdlNRi,3), 'b.')
    plot3(score(mdlNRi,1), score(mdlNRi,2), score(mdlNRi,3), 'b.')
    %plot3(pcs(NmdlWi,1), pcs(NmdlWi,2), pcs(NmdlWi,3), 'g.')
    plot3(score(mdlWi,1), score(mdlWi,2), score(mdlWi,3), 'g.')
    %plot3(pcs(NmdlREMi,1), pcs(NmdlREMi,2), pcs(NmdlREMi,3), 'r.')
    plot3(score(mdlREMi,1), score(mdlREMi,2), score(mdlREMi,3), 'r.')
    %plot3(pcs(NmdlCi,1), pcs(NmdlCi,2), pcs(NmdlCi,3), 'm.')
    plot3(score(mdlCi,1), score(mdlCi,2), score(mdlCi,3), 'mo')
    xlabel('1st principal component')
    ylabel('2nd principal component')
    zlabel('3rd principal component')
    grid on
    
    lgtxt = {};
    for eidx = 1:length(expnames)
        datREMi = strcmp(T.(expnames{eidx}){:, 'Scores'}, 'C');
        if sum(datREMi) > 0
            pcs = (T.(expnames{eidx}){datREMi, [false seld_features]} - mu) * coeff(:, 1:3);
            plot3(pcs(:,1), pcs(:,2), pcs(:,3), 'LineStyle', 'none', ...
                'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'Marker', animark{eidx})
            lgtxt = [lgtxt expnames{eidx}]; %#ok<AGROW>
        end
    end
    xlim([-10 10])
    ylim([-10 10])
    legend([{'NR points'}, {'W points'}, {'REM points'} {'C points'} lgtxt], 'Interpreter', 'none')
    
    %Create movie
    if par.Movie
        OptionZ.FrameRate=15;OptionZ.Duration=10;OptionZ.Periodic=true;
        CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10],'/home/umat/bognor/CondBModel_PAC3D',OptionZ)
    end
    
    clear eidx animark *coeff *score *mu *explained pcs *mdl*i datREMi lgtxt
end
