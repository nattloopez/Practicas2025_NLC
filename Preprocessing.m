% Preprocessing of resting state EEG signals
% Authors: Natalia Lopez, Julia Reina, Guiomar Niso
% Cajal Institute (CSIC), Madrid, Spain 
% January, 2026 

clc; clear;

disp("=== My script has started.")

%% Set parameters

% Directory to store brainstorm database
BrainstormDbDir = fullfile(pwd, 'brainstorm_db'); 
ReportsDir = '/home/natalia/app-local/out_reports';
DataDir    = '/home/natalia/app-local/out_preproc';

% Protocol name (we use the same name as QC protocol because we will
% directly import it from database)
ProtocolName = 'Dataset1'; 

% PSD parameters
Win_length = 4;
Win_overlap = 50;

% NOTCH FILTER FREQUENCIES
% Manually set them or leave empty to use results from find peaks (excluding alpha)
Notch_freqs = [];

% Band pass filter frequencies
High_pass = 0.3;
Low_pass = 0;

% Re-referencing
REF_operation = 'AVERAGE';

% Peak to peak parameters 
DetectionThreshold = 0;          
RejectionThreshold = 100;        
WindowLength = 1.0;              

% Artifact detection parameters
ECG_name = '';
EOG_name = 'Fp1';
cardiac_name = 'cardiac';
blink_name = 'blink';
remove_time = 0.25;
sspEOG_select = 1;
sspECG_select = 1;

% Import EEG (epoch) parameters
epoch_length = 4;
ignore_short = 1;

% Head model parameters
headModel_sourcespace = 1;  % Cortex surface
headModel_meg = 3;  % Overlapping spheres
headModel_eeg = 2;  % OpenMEEG BEM
headModel_ecog = 2;  % OpenMEEG BEM
headModel_seeg = 2;  % OpenMEEG BEM
headModel_nirs = 1; 

% Noise covariance parameters 
noiseCov_sensortypes = 'EEG';
noiseCov_target = 1;
noiseCov_dcoffset = 1;  % Block by block, to avoid effects of slow shifts in data
noiseCov_identity = 1;
noiseCov_copycond = 0;
noiseCov_copysubj = 0;
noiseCov_copymatch = 0;
noiseCov_replacefile = 1; 

% Source computation parameters
sources_inverseMethod = 'minnorm';
sources_inverseMeasure = 'amplitude';
sources_sourceOrient = 'fixed';
sources_loose = 0.2;
sources_useDepth = 1;
sources_weightExp = 0.5;
sources_weightLimit = 10;
sources_noiseMethod = 'reg';
sources_noiseReg = 0.1;
sources_snrMethod = 'fixed';
sources_snrRms = 1e-06;
sources_snrFixed = 3;
sources_computeKernel = 1;
sources_dataTypes = 'EEG';

% Figure size (for report)
fig_width = 1000;
fig_height = 600;
fig_size = [200, 200, fig_width, fig_height];


%% Start Brainstorm and load protocol from database
disp('== Start Brainstorm (nogui)');

% Start Brainstorm in nogui mode if not already running
if ~brainstorm('status')
    brainstorm nogui
end

% Set the Brainstorm database directory
bst_set('BrainstormDbDir', BrainstormDbDir);

% Reset colormaps
bst_colormaps('RestoreDefaults', 'eeg');

% See Tutorial 1
disp(['- BrainstormDbDir:', bst_get('BrainstormDbDir')]);
disp(['- BrainstormUserDir:', bst_get('BrainstormUserDir')]); % HOME/.brainstorm (operating system)
disp(['- HOME env:', getenv('HOME')]);
disp(['- HOME java:', char(java.lang.System.getProperty('user.home'))]);

% Verify
disp(['BrainstormDbDir: ', bst_get('BrainstormDbDir')]);  

% Load the protocol by its name
disp(['== Loading protocol: ', ProtocolName]);
iProtocol = bst_get('Protocol', ProtocolName);
if isempty(iProtocol)
    error(['Protocol "', ProtocolName, '" not found in database at: ', BrainstormDbDir]);
end

% Set as current
gui_brainstorm('SetCurrentProtocol', iProtocol);

% Start report
bst_report('Start');


%% Store data in sFilesAll for normal use
disp('== Get all raw EEG data files in protocol');

% Select all continuous raw files (or whichever file type you want)
sFilesAll = bst_process('CallProcess', 'process_select_files_data', [], [], ...
    'subjectname',   '', ...     % empty = all subjects
    'condition',     '', ...     % empty = all conditions
    'tag',           '', ...     
    'includebad',    1, ...      
    'includeintra',  1, ...
    'includecommon', 0);

% Check if we have files
if isempty(sFilesAll)
    error('No data files found in this protocol');
else
    disp(['Found ', num2str(length(sFilesAll)), ' data files in protocol.']);
end

% Process: Ignore file names with tag: epoch
sFilesCont = bst_process('CallProcess', 'process_select_tag', sFilesAll, [], ...
    'tag',    'epoch', ...
    'search', 2, ...  % Search the file names
    'select', 2);  % Ignore the files with the tag


%% Loop over participants
for iSub = 1:length(sFilesCont)

    participant = sFilesCont(iSub).SubjectName; % Extract participant information
    disp(['=== Processing participant: ', num2str(participant)]);
    
    % Select participant
    sFilesRaw = sFilesCont(iSub);

    % Create structure for JSON report
    jsonData = struct();

    % Basic fields
    jsonData.participant = participant;
    jsonData.protocol = ProtocolName;
    jsonData.date = datestr(now, 'yyyy-mm-dd HH:MM:SS');

    % Reset notch frequencies
    Notch_freqs = [];
    

    %% Preprocessing
    
    % Compute frequencies for notch filter if not specified above
    if isempty(Notch_freqs)
        % Load json report from QC step and verify
        jsonPath = fullfile('out_reports', sprintf('QC-%s-%s.json', participant, ProtocolName)); 
        if ~exist(jsonPath, 'file')
            error('QC JSON file not found for participant %s. Path: %s', participant, jsonPath);
        end

        % Read file and store frequencies above 20Hz
        QC_data = jsondecode(fileread(jsonPath));
        allPeaks = QC_data.peakFrequencies;
        Notch_freqs = allPeaks(allPeaks > 20);

        % Validate
        if isempty(Notch_freqs)
            disp('Warning: No high-frequency peaks found in JSON. Using default 50Hz.');
            Notch_freqs = 50;
        else
            disp('=== Loaded notch frequencies from QC report. ');
        end
    end


    disp('=== Apply notch and band pass filter')
    % Process: Notch filter: 50Hz 60Hz 125Hz 150Hz 187Hz
    sFilesNotch = bst_process('CallProcess', 'process_notch', sFilesRaw, [], ...
        'sensortypes', 'EEG', ...
        'freqlist',    Notch_freqs, ...
        'cutoffW',     2, ...
        'useold',      0, ...
        'read_all',    1);
    
    % Process: band-pass filter
    sFilesBand = bst_process('CallProcess', 'process_bandpass', sFilesNotch, [], ...
        'sensortypes', 'EEG', ...
        'highpass',    High_pass, ...
        'lowpass',     Low_pass, ...
        'tranband',    0, ...
        'attenuation', 'strict', ...  
        'ver',         '2019', ... 
        'mirror',      0, ...
        'read_all',    1);
    
    disp('=== Apply montage and re-reference')
    % Process: Apply montage
    sFilesEEG = bst_process('CallProcess', 'process_montage_apply', sFilesBand, [], ...
        'montage',    'Average reference', ...
        'usectfcomp', 1, ...
        'usessp',     1);
    
    % Process: Re-reference EEG
    sFilesREF = bst_process('CallProcess', 'process_eegref', sFilesEEG, [], ...
        'eegref',      REF_operation, ...
        'sensortypes', 'EEG');
    
 
    disp('=== Detect artifacts and apply SSP')

    if ~isempty(ECG_name)
        % Process: Detect heartbeats
        bst_process('CallProcess', 'process_evt_detect_ecg', sFilesREF, [], ...
            'channelname', ECG_name, ...
            'timewindow',  [], ...
            'eventname',   cardiac_name);
    end
    
    if ~isempty(EOG_name)
        % Process: Detect eye blinks
        bst_process('CallProcess', 'process_evt_detect_eog', sFilesREF, [], ...
            'channelname', EOG_name, ...
            'timewindow',  [], ...
            'eventname',   blink_name);
    end

    if ~isempty(ECG_name) && ~isempty(EOG_name)
        % Process: Remove simultaneous
        bst_process('CallProcess', 'process_evt_remove_simult', sFilesREF, [], ...
            'remove', cardiac_name, ...
            'target', blink_name, ...
            'dt',     remove_time, ...
            'rename', 0);
    end

    if ~isempty(ECG_name)
        % Process: SSP ECG: cardiac
        ECGprojectors = bst_process('CallProcess', 'process_ssp_eog', sFilesREF, [], ...
            'eventname',   cardiac_name, ...
            'sensortypes', 'EEG', ...
            'usessp',      1, ...
            'select',      sspECG_select);
    
        % Process: Snapshot: SSP projectors
        bst_process('CallProcess', 'process_snapshot', ECGprojectors, [], ...
            'type',           'ssp', ...  % SSP projectors
            'modality',       4, ...  % EEG
            'orient',         1, ...  % left
            'time',           0, ...
            'contact_time',   [0, 0.1], ...
            'contact_nimage', 12, ...
            'threshold',      30, ...
            'surfsmooth',     30, ...
            'freq',           0, ...
            'rowname',        '', ...
            'mni',            [0, 0, 0], ...
            'Comment',        '');
    end

    % Add figure size modification when changing to show the SSP contact sheet! 
    
    if ~isempty(EOG_name)
        % Process: SSP EOG: blink
        EOGprojectors = bst_process('CallProcess', 'process_ssp_eog', sFilesREF, [], ...
            'eventname',   blink_name, ...
            'sensortypes', 'EEG', ...
            'usessp',      1, ...
            'select',      sspEOG_select);
    
        % Process: Snapshot: SSP projectors
        bst_process('CallProcess', 'process_snapshot', EOGprojectors, [], ...
            'type',           'ssp', ...  % SSP projectors
            'modality',       4, ...  % EEG
            'orient',         1, ...  % left
            'time',           0, ...
            'contact_time',   [0, 0.1], ...
            'contact_nimage', 12, ...
            'threshold',      30, ...
            'surfsmooth',     30, ...
            'freq',           0, ...
            'rowname',        '', ...
            'mni',            [0, 0, 0], ...
            'Comment',        '');
    end


    disp('== Detect bad channels and segments')
    % Process: Detect bad segments/trials: Peak-to-peak  EEG(0-100)
    sFilesPtP = bst_process('CallProcess', 'process_detectbad', sFilesREF, [], ...
        'timewindow', [], ...
        'meggrad',    [0, 0], ...
        'megmag',     [0, 0], ...
        'eeg',        [DetectionThreshold, RejectionThreshold], ...
        'ieeg',       [0, 0], ...
        'eog',        [0, 0], ...
        'ecg',        [0, 0], ...
        'win_length', WindowLength, ...
        'rejectmode', 2);  % Reject the entire segments/trials (all channels)
    
    
    %% Compute parameters for json report: bad channels, segments, blinks and cardiacs
    % Locate links to raw data
    rawFiles = bst_process('CallProcess','process_select_files_data', [], [], ...
        'subjectname',   participant, ...
        'condition',     '', ...
        'tag',           '', ...
        'includebad',    1, ...
        'includeintra',  0, ...
        'includecommon', 0);
    
    % Find the average raw file
    sRaw = [];
    for f = 1:length(rawFiles)
        if contains(rawFiles(f).FileName, 'Average_reference')
            sRaw = rawFiles(f);
            break;
        end
    end
    
    % Verify file was found
    if isempty(sRaw)
        warning('No raw link file found for subject %s', participant);
    else
        % Load the raw data structure from the database
        sRawData = in_bst_data(sRaw.FileName, 'F'); 
        
        % Bad channels are contained in the channelflag
        ChannelFlags = sRawData.F.channelflag;  % 1=good, -1=bad
        numBadChannels = 0;
        BadChannels = [];
        
        % Count and store indexes of bad channels
        for i = 1:length(ChannelFlags)
            if ChannelFlags(i) == -1
                numBadChannels = numBadChannels + 1;
                BadChannels(end+1) = i;
            end
        end
        % Display results and transfer to json structure
        fprintf("Found %i bad channels.\n", numBadChannels);
        jsonData.numBadChannels = numBadChannels;
        jsonData.badChannelIndexes = BadChannels;  % Conectar con tsv para poner nombres de canales en vez de numeros?
    end
    
    % Count number of events (bad segments, blinks and cardiacs)
    % Initialize
    numBadSegments = 0;
    BadSegments     = {};   
    numBlinks       = 0;
    Blinks          = {};   
    numCardiac      = 0;
    Cardiac         = {};   
    
    % Iterate through events
    for i = 1:length(sRawData.F.events)
        evt      = sRawData.F.events(i);
        evtTimes = evt.times;   
        evtLabel = evt.label;
    
        if contains(evtLabel, 'BAD', 'IgnoreCase', true)  %
            numBadSegments = numBadSegments + length(evtTimes);
            BadSegments{end+1} = evtTimes;  
    
        elseif contains(evtLabel, 'blink', 'IgnoreCase', true)
            numBlinks = numBlinks + length(evtTimes);
            Blinks{end+1} = evtTimes;      
    
        elseif contains(evtLabel, 'cardiac', 'IgnoreCase', true)
            numCardiac = numCardiac + length(evtTimes);
            Cardiac{end+1} = evtTimes;      
        end
    end
    
    % Add to JSON structure
    jsonData.numBadSegments = numBadSegments;
    jsonData.numBlinks      = numBlinks;
    jsonData.numCardiac     = numCardiac;
    jsonData.badSegments    = BadSegments;
    jsonData.blinkTimes     = Blinks;
    jsonData.cardiacTimes   = Cardiac;


    %% PSD post processing
    
    disp('=== Compute PSD post')
    % Process: Power spectrum density (Welch)
    sFilesPSDpost = bst_process('CallProcess', 'process_psd', sFilesPtP, [], ...
        'timewindow',  [], ...
        'win_length',  Win_length, ...
        'win_overlap', Win_overlap, ...
        'units',       'physical', ...  % Physical: U2/Hz
        'sensortypes', '', ...
        'win_std',     0, ...
        'edit',        struct(...
             'Comment',         'Power', ...
             'TimeBands',       [], ...
             'Freqs',           [], ...
             'ClusterFuncTime', 'none', ...
             'Measure',         'power', ...
             'Output',          'all', ...
             'SaveKernel',      0));
    
    % Process: Snapshot: Frequency spectrum
    % bst_process('CallProcess', 'process_snapshot', sFilesPSDpost, [], ...
    %     'type',           'spectrum', ...  % Frequency spectrum
    %     'modality',       4, ...  % EEG
    %     'orient',         1, ...  % left
    %     'time',           0, ...
    %     'contact_time',   [0, 0.1], ...
    %     'contact_nimage', 12, ...
    %     'threshold',      30, ...
    %     'surfsmooth',     30, ...
    %     'freq',           0, ...
    %     'rowname',        '', ...
    %     'mni',            [0, 0, 0], ...
    %     'Comment',        '');

    % View PSD manually to change size
    hFigPSD = view_spectrum(sFilesPSDpost.FileName, 'Spectrum');
    set(hFigPSD, 'Position', fig_size); 
    bst_report('Snapshot', hFigPSD, sFilesPSDpost.FileName, 'PSD Spectrum', fig_size);
    close(hFigPSD);
   

    %% Import recordings: convert to epochs for mean and std computation
    disp('=== Import recordings: convert to epochs')

    % Process: Import MEG/EEG: Time
    sFilesEpoch = bst_process('CallProcess', 'process_import_data_time', sFilesPtP, [], ...
        'subjectname',   participant, ...
        'condition',     'epoch', ...
        'timewindow',    [], ...
        'split',         epoch_length, ...
        'ignoreshort',   ignore_short, ...
        'usectfcomp',    0, ...
        'usessp',        1, ...
        'freq',          [], ...
        'baseline',      [], ...
        'blsensortypes', 'MEG, EEG');


    % Process: Add tag: epoch
    bst_process('CallProcess', 'process_add_tag', sFilesEpoch, [], ...
        'tag',            'epoch', ...
        'output',         'name');  % Add to file name
    

    %% Compute average and standard deviation for each sensor after processing
    disp('=== Compute average and stdev for each sensor')

    % Extract number of channels from one epoch
    sData = in_bst_data(sFilesEpoch(1).FileName, 'F');
    [nChan, nTime] = size(sData.F);
    allData = zeros(nChan, nTime * length(sFilesEpoch)); % eg: 19 channels x (2000 times x 198 epochs)
    
    % Go through all epochs to obtain a single matrix with all data
    for i = 1:length(sFilesEpoch)
        sEpoch = in_bst_data(sFilesEpoch(i).FileName, 'F');
        
        % Determine the indices for this epoch in the general matrix
        startId = (i-1) * nTime + 1;
        endId   = i * nTime;
        
        % Store in the matrix
        allData(:, startId:endId) = sEpoch.F;
    end
    
    % Compute average and standard deviation across the concatenated matrix
    avgData = mean(allData, 2);
    stdev   = std(allData, 0, 2);

    % Transfer to JSON structure
    jsonData.avgPerChannelPost = avgData;
    jsonData.stdPerChannelPost = stdev;

    disp('=== Finished computing average and stdev for each sensor')


    %% Generate 2D plot to visualize mean and std in each sensor
    disp('=== Show mean and std for each sensor')

    % view_topography requires a Brainstorm structure, so we use an
    % existing one but actually pass avgData and stdev to plot

    % First: mean figure
    hFigMean = view_topography(sFilesEpoch(1).FileName, 'EEG', '2DDisc', avgData, 1, "NewFigure");
    bst_report('Snapshot', hFigMean, sFilesEpoch(1).FileName, 'Mean per Channel (Post-processing)', fig_size);
    close(hFigMean);
    
    % Second: standard deviation
    hFigStd = view_topography(sFilesEpoch(1).FileName, 'EEG', '2DDisc', stdev, 1, "NewFigure");
    bst_report('Snapshot', hFigStd, sFilesEpoch(1).FileName, 'Standard deviation per Channel (Post-processing)', fig_size);
    close(hFigStd);

    disp('=== Finished displaying mean and std for each sensor')


    %% Compute sources
    
    disp('=== Compute head model') 
    % Process: Compute head model
    sFilesHead = bst_process('CallProcess', 'process_headmodel', sFilesREF, [], ...
        'Comment',     '', ...
        'sourcespace', headModel_sourcespace, ...  
        'meg',         headModel_meg, ...  
        'eeg',         headModel_eeg, ...  
        'ecog',        headModel_ecog, ...  
        'seeg',        headModel_seeg, ...  
        'nirs',        headModel_nirs, ...  
        'openmeeg',    struct(...
             'BemSelect',    [1, 1, 1], ...
             'BemCond',      [1, 0.0125, 1], ...
             'BemNames',     {{'Scalp', 'Skull', 'Brain'}}, ...
             'BemFiles',     {{}}, ...
             'isAdjoint',    0, ...
             'isAdaptative', 1, ...
             'isSplit',      0, ...
             'SplitLength',  4000), ...
        'nirstorm',    struct(...
             'FluenceFolder',    'https://neuroimage.usc.edu/resources/nst_data/fluence/', ...
             'smoothing_method', 'geodesic_dist', ...
             'smoothing_fwhm',   10), ...
        'channelfile', '');
    
    disp('=== Compute noise covariance matrix')
    % Process: Compute covariance (noise or data)
    sFilesNoise = bst_process('CallProcess', 'process_noisecov', sFilesREF, [], ...
        'baseline',       [], ...
        'datatimewindow', [], ...
        'sensortypes',    noiseCov_sensortypes, ...
        'target',         noiseCov_target, ...  
        'dcoffset',       noiseCov_dcoffset, ...  
        'identity',       noiseCov_identity, ...
        'copycond',       noiseCov_copycond, ...
        'copysubj',       noiseCov_copysubj, ...
        'copymatch',      noiseCov_copymatch, ...
        'replacefile',    noiseCov_replacefile);  
    
    % Process: Snapshot: Noise covariance
    % sFiles = bst_process('CallProcess', 'process_snapshot', sFilesNoise, [], ...
    %     'type',           'noiscov', ...  % Noise covariance
    %     'modality',       4, ...  % EEG
    %     'orient',         1, ...  % left
    %     'time',           0, ...
    %     'contact_time',   [0, 0.1], ...
    %     'contact_nimage', 12, ...
    %     'rowname',        '', ...
    %     'Comment',        '');

    % View noise covariance manually to change the size
    [sStudy, iStudy] = bst_get('AnyFile', sFilesNoise.FileName);
    NoiseCovFile = sStudy.NoiseCov.FileName;
    hFigNoise = view_noisecov(NoiseCovFile);
    set(hFigNoise, 'Position', fig_size); 
    bst_report('Snapshot', hFigNoise, NoiseCovFile, 'Noise covariance', fig_size);
    close(hFigNoise);

    
    disp('=== Compute sources')
    % Process: Compute sources [2018]
    sFilesSources = bst_process('CallProcess', 'process_inverse_2018', sFilesREF, [], ...
        'output',  1, ...  % Kernel only: shared
        'inverse', struct(...
             'Comment',        'MN: EEG', ...
             'InverseMethod',  sources_inverseMethod, ...
             'InverseMeasure', sources_inverseMeasure, ...
             'SourceOrient',   {{sources_sourceOrient}}, ...
             'Loose',          sources_loose, ...
             'UseDepth',       sources_useDepth, ...
             'WeightExp',      sources_weightExp, ...
             'WeightLimit',    sources_weightLimit, ...
             'NoiseMethod',    sources_noiseMethod, ...
             'NoiseReg',       sources_noiseReg, ...
             'SnrMethod',      sources_snrMethod, ...
             'SnrRms',         sources_snrRms, ...
             'SnrFixed',       sources_snrFixed, ...
             'ComputeKernel',  sources_computeKernel, ...
             'DataTypes',      {{sources_dataTypes}}));
    
    
    disp('=== Finished computing sources')
    

    %% Save report and end loop
    
    % Save report
    disp('=== Save report');
    % Desired filename
    outputName = fullfile(ReportsDir, sprintf('Preprocessing-%s-%s.html', participant, ProtocolName));
    
    % Save and then export to the custom name
    ReportFile = bst_report('Save', []);
    bst_report('Export', ReportFile, outputName);
    
    % Create JSON report
    jsonFile = fullfile(ReportsDir, sprintf('Preproc-%s-%s.json', participant, ProtocolName));
    fid = fopen(jsonFile, 'w');
    if fid == -1
        error("Cannot open JSON file.")
    end 
    fprintf(fid, '%s', jsonencode(jsonData, PrettyPrint=true));
    fclose(fid);
    
    disp(['Reports exported as: ', outputName, jsonFile]);

end

% Copy brainstorm_db data
copyfile([BrainstormDbDir,'/',ProtocolName], DataDir);


%% DONE
disp('=== Done!');
