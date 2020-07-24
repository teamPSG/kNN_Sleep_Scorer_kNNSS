function [T, eegchanno, emgchanno, epdur, epdat, Fs_eeg, Fs_emg, eeg_band, emg_band, parameters] = ...
    generate_statespace(top_folder, filename, varargin)

%Usage:
%  [T, EpochedData, EEGChannelNumber, EMGChannelNumber, EpochDuration] = ...
%     generate_statespce(top_folder, filename, ...)
%
%Description:
%  This function reads an EDF file from disk as well as the associated
%  "FFT" text file (if training set is to be created, see below) that
%  contains manual scoring. Calculates metrics for classification from EDF
%  data and concatenates it with manual scoring values (if specified, see
%  below). Metrics calculation is implemented in other functions
%  that should be called within genetare_trainingset. Presently called
%  metrics are:
%   - Power in frequency bands for EEG & EMG signals (relative & absolute)
%   - Hjorth parameters for EEG & EMG signals
%
%Notes:
%  - Only one EEG and one EMG channel is used. Trivial modifications are
%    needed to accommodate for more.
%  - Sampling rates on the two channels are assumed to be equal. Some
%    modifications, including checking readedf() is required to change
%    this.
%
%Inputs:
%  top_folder: string, specifies the absolute path to the data. Unless the
%    optional input arguments inEDF and inFFT are specified the function
%    assumes that EDF files are stored in a subfolder called EDF and FFT
%    text files in a subfolder called FFT
%  filename: string, the base name of the files that describe the
%    experiment. The '.edf', and '.txt' extensions are appended to the file
%    names later. EDF and FFT base file names on disk must match.
%
%Outputs:
%  T: table, contains manual score as string in its first column and
%    metrics data in subsequent columns. Table header gives metrics names.
%  EpochedData: NxMxL double array, where N is the number of channels, M is
%    number time points in a single epoch, and L is the number of epochs
%    used for calculating features. Note that EpochedData might not contain
%    the whole length of data, depending on features calculated.
%  EEGChannelNumber: double, the number of channel in EDF file used for EEG
%  EMGChannelNumber: double, the number of channel in EDF file used for EMG
%  EpochDuration: double, length of an epoch in seconds
%  If WriteTmpFile is set, a file is also generated with the features and
%    scores concatenated
%
%Optional input arguments as name-value pairs:
%IO-related
%  'WriteTmpFile': string, name of a temporary file containing the training
%    set that will be written to disk. If left unspecified temporary file
%    will not be written.
%  'SavePSpecFile': string, specifies name of file used to store epoched
%    power spectrum. Power spectrum is not saved if left unspecified
%    (default).
%  'OutDir': string, path to output directory, defult is the working
%    directory (.)
%  'EDFFolder': string, storage subfolder of .edf files under top_folder.
%    Defaults to 'EDF'
%  'FFTFolder': string, storage subfolder of .txt files under top_folder
%    containing manual scoring. Defaults to 'FFT'
%  'ManScoreVar': Nx1 cell of strings, manual scores. A long awaited
%    addition to the function allowing to pass a variable with scores
%    instead of reading a file. N is the number of epochs scored. Default
%    is {}. If ManScoreType == 'SOTE' scores are expected to be passed to
%    the function in a variable specified here.
%
%Data-related
%  'DataFileType': string, describes what file format the raw input data is
%    in. Originally only EDF was planned but life is nonlinear. Possible
%    values now are: 'EDF' (default) or 'D1MAT'.
%  'ManScoreType': string, can be SRI or other pre-defined value to
%    describe what type of manual score there available. Please see Section 
%    "Load FFT file and extract manual scoring and epoch duration" in the
%    code for available options and how to add further types. Default is
%    SRI.
%  'EpDur': double, the length of a manually scored epoch in sec. This
%    value will set the epoch duration and will be used for metric
%    calculation. The function will first try to read this value from the
%    FFT file, if fails will set it to 10 sec (default for EpDur).
%  'SlideTime': double, step size of epoch window in sec. In case of using
%    non-overlapping windows (typical for sleep scoring) set it to be equal
%    with EpDur. Default value is 10 sec.
%
%Analysis-related
%  'MTFreq': 1x2 double, frequency range of multitaper calculation given in
%    Hz. Defaults to [0.5 80].
%  'Tapers': 1x2 double, a numeric vector [TW K] where TW is the
%    time-bandwidth product and K is the number of tapers to be used (less
%    than or equal to 2TW-1). Default is [5 9].
%  'Description': string, a free text description of the trainingset.
%    Default: ''
%  'CalcFeatures': logical, if true the features are calculated, if false
%    features are not calculated only Scores are returned in T. This
%    operation makes sense if the epoched data is to be used with Scores
%    associated to each epoch. Default is true.
%  'EEGChannel2Use': integer, if more than one EEG channels are present in
%    the input data this will choose the one to use for training the
%    classifier. Default is 1.
%  'EMGCHannel2Use': integer, same as above with EMG.
%  'EEGChannelName': string, name of the EEG channel to be used for
%    calculations. If multiple channels have the same name use
%    EEGChannel2Use below to select the one to be used. Default is 'EEG'.
%  'EMGChannelName': string, same as above with EMG.
%  'EEGChannelNo': double, explicitely specifies the channel number of the
%    EEG channel to be used for feature calculations. If unspecified the
%    string in 'EEGChannelName' will be used to determine which EDF file
%    channel contains EEG information.
%  'EMGChannelNo': double, same as above for the EMG channel
%  'EEGMontage': 1x2 cell of strings, names of 2 EEG channels that are used
%    to calculate a difference betwen the two and the difference channel
%    will be used as EEGChannel
%  'FirstEpoch': integer, if scoring started later than start of file, this
%    parameter can take care of it. Default is 1.
%  'LastEpoch': integer, if scoring finished before end of recording, this
%    parameter can take care of it. Default is Inf, meaning to use up to
%    the last epoch.
%  'IsTrainingSet': logical, if set to true, the function reads the manual
%    scores associated with the EDF files and creates a unified table
%    containing features -- scores pairs. Default is false.
%
%Dependencies:
%  This function uses edfread by Brett Shoelson
%  This function uses hjorth by Alois Schloegl (part of BIOSIG-toolbox)
%
%Author: Tamas Kiss <kiss.t@wigner.hu>

%% Parse input and set default parameters
p = inputParser;
addRequired(p, 'top_folder', @isstr);
addRequired(p, 'filename', @isstr);
addParamValue(p, 'WriteTmpFile', false, @isstr); %#ok<*NVREPL>
addParamValue(p, 'SavePSpecFile', false, @isstr);
addParamValue(p, 'OutDir', '.', @isstr);
addParamValue(p, 'EDFFolder', 'EDF', @isstr);
addParamValue(p, 'FFTFolder', 'FFT', @isstr);
addParamValue(p, 'ManScoreVar', {}, @iscell);
addParamValue(p, 'DataFileType', 'EDF', @isstr);
addParamValue(p, 'ManScoreType', 'SRI', @isstr);
addParamValue(p, 'EpDur', 10, @isnumeric);
addParamValue(p, 'SlideTime', 10, @isnumeric);
addParamValue(p, 'MTFreq', [0.1 120], @(x)validateattributes(x, {'numeric'}, {'size', [1 2]}));
addParamValue(p, 'Tapers', [5 9], @(x)validateattributes(x, {'numeric'}, {'size', [1 2]}));
addParamValue(p, 'Description', '', @isstr);
addParamValue(p, 'CalcFeatures', true, @islogical);
addParamValue(p, 'EEGChannel2Use', 1, @isnumeric);
addParamValue(p, 'EMGChannel2Use', 1, @isnumeric);
addParamValue(p, 'EEGChannelName', 'EEG', @isstr);
addParamValue(p, 'EMGChannelName', 'EMG', @isstr);
addParamValue(p, 'EEGChannelNo', [], @isnumeric);
addParamValue(p, 'EMGChannelNo', [], @isnumeric);
addParamValue(p, 'EEGMontage', {}, @iscell);
addParamValue(p, 'FirstEpoch', 1, @isnumeric);
addParamValue(p, 'LastEpoch', Inf, @isnumeric);
addParamValue(p, 'IsTrainingSet', false, @islogical);


parse(p, top_folder, filename, varargin{:});

%% Further parameter setting used for metrics calculation [see notes inside]

%These settings are used later when features are calculated. If further
%metrics are introduced that require parameters to be set, please specify
%them here.

eeg_band.names = {'eeg_delta','eeg_delta_low','eeg_delta_high','eeg_theta','eeg_alpha','eeg_sigma','eeg_beta','eeg_gamma'};
eeg_band.ranges=[0.5,4; 0.5,2.5; 2.5,4.0; 4,8; 8,12; 12,16; 16,24; 24,50];
emg_band.names = {'emg_delta','emg_theta','emg_alpha','emg_sigma','emg_beta','emg_gamma','emg_high'};
emg_band.ranges=[0.5,4; 4,8; 8,12; 12,16; 16,24; 24,50; 80,120];
relat_eeg_band.names = strrep(eeg_band.names, 'eeg_', 'eeg_R_');
if length(eeg_band.names) ~= length(eeg_band.ranges) || length(emg_band.names) ~= length(emg_band.ranges)
    fprintf('generate_statespace:: Error: inappropriate parameters set.')
    T = NaN;
    return
end
eeg_band.names{end+1} = 'eeg_RMS';
emg_band.names{end+1} = 'emg_RMS';

if ~p.Results.IsTrainingSet
    epdur = p.Results.EpDur;
    mst = 'NoScore';
    noscore = true;
else
    mst = p.Results.ManScoreType;
    noscore = false;
end

%% Load FFT file and extract manual scoring and epoch duration [see notes inside]

%This section is very use-dependent. Different people and organizations
%save their manual scores in different file formats. If you have your own
%format, please specify how to read it here in a new 'case'.

switch mst
    case 'SRI'
        warning('off', 'MATLAB:table:ModifiedVarnames')
        T = readtable([top_folder filesep p.Results.FFTFolder filesep filename '.txt'], 'Delimiter', '\t', 'ReadVariableNames', false, 'HeaderLines', 10);
        warning('on', 'MATLAB:table:ModifiedVarnames')
        scores = T(:,2);
        scores.Properties.VariableNames = {'Scores'};
        fid = fopen([top_folder filesep p.Results.FFTFolder filesep filename '.txt']);
        fgetl(fid); fgetl(fid); fgetl(fid);
        tmp = fgetl(fid);
        fclose(fid);
        tmp = str2double(tmp(17:18));
        if (~isnan(tmp))
            epdur = tmp;
        else
            epdur = p.Results.EpDur;
            fprintf('generate_statespace:: Warning: could not read epoch duration interval from FFT file. Using %2.0f sec instead.\n', epdur);
        end
        clear T
    case 'CiTox'
        warning('off', 'MATLAB:table:ModifiedVarnames')
        T = readtable([top_folder filesep p.Results.FFTFolder filesep filename '.csv']);
        warning('on', 'MATLAB:table:ModifiedVarnames')
        T(1,:) = [];
        scores = T(:, 'LargeAnimalSleep');
        scores.Properties.VariableNames = {'Scores'};
        epdur = etime(datevec(T{1,2}), datevec(T{1,1}));
        clear T
    case 'Irina'
        warning('off', 'MATLAB:table:ModifiedVarnames')
        if strcmp(p.Results.DataFileType, 'D1MAT')
            filename = strrep(filename, '.mat', '');
        end
        T = readtable([top_folder filesep p.Results.FFTFolder filesep filename '.txt'], 'Delimiter', '\t', 'ReadVariableNames', false);
        warning('on', 'MATLAB:table:ModifiedVarnames')
        scores = T(3:end,2);
        scores.Properties.VariableNames = {'Scores'};
        epdur = p.Results.EpDur;
        recstart = datevec(T{1,2});
        scorestart = datevec(T{2,2});
        dataskip = scorestart - recstart;
        dataskip = sum(dataskip(4:6).*[3600 60 1]); % given in seconds since Fs is not yet known
        fprintf('\t%i scores will be used. Recording start @%s, scoring start @%s, skipping %i sec.\n', ...
            height(scores), T{1,2}{:}, T{2,2}{:}, dataskip)
        clear T
    case 'TSK'
        load([top_folder filesep p.Results.FFTFolder filesep filename '.mat'])
        scores = cell2table(scores, 'VariableNames', {'Scores'}); %#ok<NODEF>
        epdur = p.Results.EpDur;
    case 'SOTE'
        scores = cell2table(p.Results.ManScoreVar, 'VariableNames', {'Scores'});
        epdur = p.Results.EpDur;
    case 'SimpleCSV'
        epdur = p.Results.EpDur;
        msif = [top_folder filesep p.Results.FFTFolder filesep filename '.csv'];
        fprintf('Reading manual score from %s...', msif)
        T = readtable(msif, 'ReadVariableNames', false, 'Delimiter', 'tab');
        scores = T(:, 'Var2');
        fprintf(' read %i records (%0.2f min @ epdur = %0.2f sec).\nScores found:\n', height(scores), height(scores)*epdur/60, epdur)
        tabulate(scores{:,:})
        scores.Properties.VariableNames = {'Scores'};
        clear T
    case 'NoScore'
        % Do nothing
    otherwise
        fprintf('generate_statespace:: Error: unknown manual score type %s.\n', p.Results.ManScoreType)
        T= NaN;
        epdat = NaN;
        return;
end

%% Load data, get sampling rate

%In this section data is loaded and using labels in the header we try to
%guess which channel is the EEG and which is the EMG channel. If there are
%multiple EEG channels either one of them can be used for further
%calculations or a montage (one minus the other) can be used.

if strcmp(p.Results.DataFileType, 'EDF')
    edfif = [top_folder filesep p.Results.EDFFolder filesep filename '.edf'];
    fprintf('Reading EDF file: %s...', edfif)
    [header, data] = edfread(edfif);
    fprintf(' done.\n')
    if isempty(p.Results.EEGMontage)
        if isempty(p.Results.EEGChannelNo)
            eegchanno = find(strcmpi(header.label, p.Results.EEGChannelName));
        else
            eegchanno = p.Results.EEGChannelNo;
        end
    else
        if length(p.Results.EEGMontage) == 2
            eegchanno(1) = find(strcmp(header.label, p.Results.EEGMontage{1}));
            eegchanno(2) = find(strcmp(header.label, p.Results.EEGMontage{2}));
        else
           fprintf('generage_statespace:: Error: EEGMontage has to contain exactly 2 names.\n')
           T = NaN;
           return 
        end
    end
    if isempty(p.Results.EMGChannelNo)
        emgchanno = find(strcmpi(header.label, p.Results.EMGChannelName));
    else
        emgchanno = p.Results.EMGChannelNo;
    end
    %eegchanno = find(not(cellfun(@isempty, strfind(header.label, p.Results.EEGChannelName)))); %This approach would be better but not backward compatible
    %emgchanno = find(not(cellfun(@isempty, strfind(header.label, p.Results.EMGChannelName))));
    Fs_eeg = (header.samples(eegchanno)/header.duration);
    Fs_emg = (header.samples(emgchanno)/header.duration);
    if length(Fs_eeg) > 1
        if Fs_eeg(1) ~= Fs_eeg(2)
            fprintf('generage_statespace:: Error: EEG channels to be used in montage have different sampling rate.\n')
            T = NaN;
            return 
        else
            Fs_eeg = Fs_eeg(1);
        end
    end
elseif strcmp(p.Results.DataFileType, 'D1MAT')
    if isempty(p.Results.EEGMontage) || length(p.Results.EEGMontage) > 2
        fprintf('generage_statespace:: Error: EEGMontage cannot be empty or have more than 2 elements with DataFileType D1MAT.\n')
        T = NaN;
        return
    end
    load([top_folder filesep filename]) 
    emgchanno = find(ismember(description{ismember(description(:,1), 'Channel names'), 2}, p.Results.EMGChannelName)); %#ok<NODEF>
    eegchanno = find(ismember(description{ismember(description(:,1), 'Channel names'), 2}, p.Results.EEGMontage));
    Fs_eeg = description{ismember(description(:,1), 'Sampling rate'), 2};
    Fs_emg = Fs_eeg;
else
    fprintf('generate_statespace:: Warning: Unknown input data file type.\n')
end

%% Check if assumptions mentioned in notes are fulfilled
if length(eegchanno) > 1 && length(p.Results.EEGMontage) < 2
    fprintf('generate_statespace:: Warning: multiple EEG channels in %s. Using #%i\n', ...
        filename, p.Results.EEGChannel2Use)
    eegchanno = eegchanno(p.Results.EEGChannel2Use);
end
if length(emgchanno) > 1
    fprintf('generate_statespace:: Warning: multiple EMG channels in %s. Using #%i\n', ...
        filename, p.Results.EMGChannel2Use)
    emgchanno = emgchanno(p.Results.EMGChannel2Use);
end
if Fs_eeg ~= Fs_emg
    fprintf('generate_statespace:: Error: sampling rates on EEG and EMG channels are different.')
    T = NaN;
    return
else
    Fs = Fs_eeg;
end

fprintf('Read %i channels and %i data points (%0.2f min @ Fs = %0.2f).\n', size(data), size(data, 2)/Fs/60, Fs) 

if strcmp(p.Results.DataFileType, 'EDF')
    if length(p.Results.EEGMontage) == 2
        fprintf('Using channels #%i, and #%i (%s, %s) for EEG in montage,\nand channel #%i (%s) for EMG calculations.\n', ...
            eegchanno, header.label{eegchanno}, emgchanno, header.label{emgchanno})
    else
        fprintf('Using channel #%i (%s) for EEG and channel #%i (%s) for EMG calculations.\n', ...
            eegchanno, header.label{eegchanno}, emgchanno, header.label{emgchanno})
    end
elseif strcmp(p.Results.DataFileType, 'D1MAT')
    if length(p.Results.EEGMontage) == 1
        fprintf('Using channel #%i (%s) for EEG and channel #%i (%s) for EMG calculations.\n', ...
            eegchanno, description{ismember(description(:,1), 'Channel names'), 2}{eegchanno}, ...
            emgchanno, description{ismember(description(:,1), 'Channel names'), 2}{emgchanno})
    else
        fprintf('Using channels #%i, #%i (%s, %s) for EEG and channel #%i (%s) for EMG calculations.\n', ...
            eegchanno, description{ismember(description(:,1), 'Channel names'), 2}{eegchanno}, ...
            emgchanno, description{ismember(description(:,1), 'Channel names'), 2}{emgchanno})
    end
end

%% Number of samples in a given epoch and sliding time with above sampling rate
dN = floor(Fs*epdur);
dS = floor(Fs*p.Results.SlideTime);

%% Cut out EEG and EMG channels only to save memory for following computation
if length(p.Results.EEGMontage) == 2
    data = [data(eegchanno(1), :) - data(eegchanno(2), :); data(emgchanno, :)];
else
    data = data([eegchanno, emgchanno], :);
end

%% Cut out pieces used for training
if p.Results.IsTrainingSet
    if round(size(data,2)/Fs/60) ~= round(height(scores)*epdur/60)
        fprintf('Attention: data file is ~%i, score data is ~%i minutes long.\n', ...
            round(size(data,2)/Fs/60), round(height(scores)*epdur/60));
    end
    switch mst
        case {'SRI' 'CiTox' 'SOTE' 'SimpleCSV'}
            if isempty(p.Results.LastEpoch)
                LastEpoch = height(scores);
            else
                LastEpoch = p.Results.LastEpoch;
            end
            if isinf(LastEpoch)
                if ~noscore
                    scores = scores(p.Results.FirstEpoch:end,:);
                end
                data = data(:,(p.Results.FirstEpoch-1)*dN+1:end);
            else
                if ~noscore
                    scores = scores(p.Results.FirstEpoch:LastEpoch,:);
                end
                data = data(:, (p.Results.FirstEpoch-1)*dN+1:LastEpoch*dN);
            end
        case 'Irina'
            %If scoring is of the Irina type, all scores are used but scoring does
            %not start at the beginning of recording and does not last up to the
            %end of the recording so beginning of data will be skipped and not all
            %of data used to generate training set
            data = data(:, round(dataskip*Fs+1):round(dataskip*Fs+height(scores)*epdur*Fs));
        case 'TSK'
            if (length(data)/Fs)/epdur > height(scores)
                data = data(:, 1:end-dN);
            end
    end
end

%% Calculate spectrogram using Chronux's multitaper method and band metric
%Note: half of Fs_X is used below because of dummy downsampling

%Power in different frequency bands are calculated here. The bands defined
%here are classical EEG frequency bands that have activity correlated with
%specific behavioral or sleep stages. The frequency interval of a given
%band is specified in the e*g_band.ranges fields in Hz, and their
%corresponding names in the e*g_band.names.
if p.Results.CalcFeatures
    %EEG spectrogram:
    [eeg_band.values, pspec, faxis] = metric_band_power(eeg_band, data(1, 1:2:end), Fs/2, [epdur p.Results.SlideTime], p.Results.MTFreq, p.Results.Tapers); 
    pspec = pspec';%#ok<NASGU> %For backward compatibility
    %EMG spectrogram:
    emg_band.values = metric_band_power(emg_band, data(2, 1:2:end), Fs/2, [epdur p.Results.SlideTime], p.Results.MTFreq, p.Results.Tapers);
    %Putting together for output
    B = array2table([eeg_band.values emg_band.values]);
    B.Properties.VariableNames = [eeg_band.names emg_band.names];
end

%% Calculate Hjorth parameters

% Hjorth parameters are indicators of statistical properties used in signal
% processing. For details see:
% https://en.wikipedia.org/wiki/Hjorth_parameters

if p.Results.CalcFeatures
    %binno = floor(length(data)/dN);
    binno = floor((length(data)-dN)/dS)+1;
    eeg_hjo = NaN(binno,3);
    emg_hjo = NaN(binno,3);
    for bidx = 1:binno
        tidxs = dS*(bidx-1)+1:dS*(bidx-1)+dN;
        [a, m, c] = hjorth(data(:, tidxs)');
        eeg_hjo(bidx, :) = [a(1) m(1) c(1)];
        emg_hjo(bidx, :) = [a(2) m(2) c(2)];
    end
    H = array2table([eeg_hjo emg_hjo]);
    H.Properties.VariableNames = {'eeg_act' 'eeg_mob' 'eeg_cplx' 'emg_act' 'emg_mob' 'emg_cplx'};
end

%% Calculate relative EEG band powers

% Relative powers are also useful features: for example an overall increase
% of EEG signal power when transitioning from wake to sleep can suppress
% certain power peaks, like theta.

if p.Results.CalcFeatures
    RB = array2table(bsxfun(@rdivide, B{:,1:length(eeg_band.names)-1}, B{:,'eeg_RMS'}));
    RB.Properties.VariableNames = relat_eeg_band.names;
end

%% Output
%Some features might not be possible to calculate for all epochs (e.g. if
%window size is greater than an epoch's length as in band data by Dmitri).
%So let's assume values start with the first epoch and end wherever they
%end and so cut all features to the same length first.

if p.Results.CalcFeatures
    if ~noscore
        noep = min([height(scores) height(B), height(H)]);
        fprintf('Number of epochs: %i (S: %i, B: %i, H: %i).\n', noep, height(scores), height(B), height(H))
        T = [scores(1:noep, 1) B(1:noep,:) H(1:noep,:) RB(1:noep, :)];
        warning('off', 'MATLAB:table:DimNameNotValidIdentifierBackCompat')
        T.Properties.DimensionNames = {'Epoch', 'Scores_Features'};
        warning('on', 'MATLAB:table:DimNameNotValidIdentifierBackCompat')
    else
        noep = min([height(B), height(H)]);
        fprintf('Number of epochs: %i (B: %i, H: %i).\n', noep, height(B), height(H))
        T = [B(1:noep,:) H(1:noep,:) RB(1:noep, :)];
        T.Properties.DimensionNames = {'Epoch', 'Features'};
    end
    T.Properties.Description = p.Results.Description;
else
    if ~noscore
        noep = min([floor(length(data)/dN); height(scores)]);
        T = scores(1:noep, 1);
        T.Properties.Description = p.Results.Description;
        T.Properties.DimensionNames = {'Epoch', 'Score'};
    else
        noep = floor(length(data)/dN);
        T = table();
    end
end
if nargout > 4
    epdat = reshape(data(:, 1:dN*noep), size(data,1), dN, noep);
end
fprintf('\t%i records out.\n', noep)

%% Write output to disk
parameters = p.Results; 
if p.Results.WriteTmpFile
    ofn = [p.Results.OutDir filesep p.Results.WriteTmpFile];
    fprintf('Writing output to %s\n', ofn)
    save(ofn, 'T', 'eeg_band', 'emg_band', 'parameters', 'eegchanno', 'emgchanno', 'epdur');
end
if p.Results.SavePSpecFile
    ofn = [p.Results.OutDir filesep p.Results.SavePSpecFile];
    fprintf('Writing power spectrum to %s\n', ofn)
    save(ofn, 'pspec', 'faxis');
end
