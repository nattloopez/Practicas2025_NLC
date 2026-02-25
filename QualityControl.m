% Quality control of resting state EEG signals
% Authors: Natalia Lopez, Julia Reina, Guiomar Niso
% Cajal Institute (CSIC), Madrid, Spain 
% January, 2026 

clc; clear;

disp("=== My script has started.")

%% Set parameters

% Directory to store brainstorm database
BrainstormDbDir = fullfile(pwd, 'brainstorm_db'); 
ReportsDir = '/home/natalia/app-local/out_reports';
DataDir    = '/home/natalia/app-local/out_QC';

% BIDS directory with all participants
BIDS_dir = '/home/natalia/app-local/EEG_BIDS_restingdata_output';

% Protocol name and parameters
ProtocolName = 'DatasetAD'; 
UseDefaultAnat = 1; 
UseDefaultChannel = 0;

% Import BIDS parameters
nVertices = 15000; 

% Import EEG (epoch) parameters
epoch_length = 4;
ignore_short = 1;

% EEG positions
% ds004504:
% EEGcap = 'ICBM152: 10-20 19';
% ds003775
% EEGcap = 'Colin27: BioSemi 64 10-10';
% Home acquisition
EEGcap = 'ICBM152: BrainProducts EasyCap 64';

% PSD parameters
Win_length = 4;
Win_overlap = 50;
peaksThreshold = 1e-13;

% Figure size (for report)
fig_width = 1000;
fig_height = 600;
fig_size = [200 200 fig_width fig_height];


%% START BRAINSTORM
disp('== Start Brainstorm defaults');

% Start Brainstorm (no GUI)
brainstorm nogui

% Set Brainstorm database directory
bst_set('BrainstormDbDir',BrainstormDbDir) 
% Reset colormaps
bst_colormaps('RestoreDefaults', 'eeg');

% See Tutorial 1
disp(['- BrainstormDbDir:', bst_get('BrainstormDbDir')]);
disp(['- BrainstormUserDir:', bst_get('BrainstormUserDir')]); % HOME/.brainstorm (operating system)
disp(['- HOME env:', getenv('HOME')]);
disp(['- HOME java:', char(java.lang.System.getProperty('user.home'))]);


%% CREATE PROTOCOL 
disp('== Create protocol');

% Create new protocol
gui_brainstorm('CreateProtocol', ProtocolName, UseDefaultAnat, UseDefaultChannel);
disp('- New protocol created');

% Start report
bst_report('Start');

%% Import BIDS data
disp('=== Import BIDS dataset')

% Process: BIDS import
bst_process('CallProcess', 'process_import_bids', [], [], ...
    'bidsdir',       {BIDS_dir, 'BIDS'}, ...
    'nvertices',     nVertices, ...            % number of vertices for cortex surface reconstruction (if FS surfaces available)
    'mni',           'maff8', ...          % method for MNI normalization, if needed
    'anatregister',  'spm12', ...          % method for anatomy registration (if you want registration)
    'groupsessions', 1, ...                % whether to group multiple sessions into one subject
    'channelalign',  1);                   % align sensor positions (if head-shape or fiducials available)

% Store files that were imported
sFilesAll = bst_process('CallProcess', 'process_select_files_data', [], [], ...
    'tag',         '', ...
    'subjectname', 'All');

% Verify
if isempty(sFilesAll)
    error("Error: Database query returned no files. Check if BIDS import actually worked.");
else
    disp(["== Success: Found " num2str(length(sFilesAll)) " files in database"]);
end


%% Loop over participants
for iSub = 1:length(sFilesAll)

    % Extract participant information
    participant = sFilesAll(iSub).SubjectName;
    disp(['=== Processing participant: ', num2str(participant)]);
    
    % Select participant
    sFilesEEG = sFilesAll(iSub);

    % Create structure for JSON report
    jsonData = struct();

    % Basic fields
    jsonData.participant = participant;
    jsonData.protocol = ProtocolName;
    jsonData.date = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    
    
    %% Import recordings: convert to epochs for later use
    disp('=== Import recordings: convert to epochs')

    % Process: Import MEG/EEG: Time
    sFilesImport = bst_process('CallProcess', 'process_import_data_time', sFilesEEG, [], ...
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
    bst_process('CallProcess', 'process_add_tag', sFilesImport, [], ...
        'tag',            'epoch', ...
        'output',         'name');  % Add to file name


    %% Quality Control: check sensors 
    disp('=== Add EEG positions')
    
    % Process: Add EEG positions
    sFilesEEG = bst_process('CallProcess', 'process_channel_addloc', sFilesEEG, [], ...
        'channelfile', {'', ''}, ...
        'usedefault',  EEGcap, ... 
        'fixunits',    1, ...
        'vox2ras',     1, ...
        'mrifile',     {'', ''}, ...
        'fiducials',   []);

    disp('=== Check correct placement of sensors.')

    % Get cortex file and display it
    sSubject = bst_get('Subject', sFilesEEG.SubjectName);
    CortexFile = sSubject.Surface(sSubject.iCortex).FileName;
    hFig = view_surface(CortexFile);
    
    % Overlay the sensors and adjust transparency
    view_channels(sFilesEEG.ChannelFile, 'EEG', 1, 1, hFig, 0, []);
    panel_surface('SetSurfaceTransparency', hFig, 1, 0.5); % 50% transparency
    
    % Set orientation left (1)
    figure_3d('SetStandardView', hFig, 'left'); 
    bst_report('Snapshot', hFig, sFilesEEG.FileName, 'Sensor Layout: Left View', fig_size);

    % Set orientation right (3)
    figure_3d('SetStandardView', hFig, 'right'); 
    bst_report('Snapshot', hFig, sFilesEEG.FileName, 'Sensor Layout: Right View', fig_size);
    
    % Set orientation top (5) 
    figure_3d('SetStandardView', hFig, 'top'); 
    bst_report('Snapshot', hFig, sFilesEEG.FileName, 'Sensor Layout: Top View', fig_size);
    
    % Close figure
    close(hFig);

    disp('=== Finished checking sensors')


    %% Quality control: check the time series
    disp('=== Show raw recordings in column and butterfly mode')

    % Open figure, set it in column format and make snapshot
    hFigColumn = view_timeseries(sFilesEEG.FileName);
    pause(0.5);  % allow for complete rendering
    panel_record('SetDisplayMode', hFigColumn, 'column'); 
    bst_report('Snapshot', hFigColumn, sFilesEEG.FileName, 'Raw signals – Column view', fig_size);
    close(hFigColumn);

    % Open figure, set it in butterfly format and make snapshot
    hFigButterfly = view_timeseries(sFilesEEG.FileName);
    pause(0.5);  % allow for complete rendering
    panel_record('SetDisplayMode', hFigButterfly, 'butterfly'); 
    bst_report('Snapshot', hFigButterfly, sFilesEEG.FileName, 'Raw signals – Butterfly view', fig_size);
    close(hFigButterfly);
    
    disp('=== Finished showing raw recordings')


    %% Quality control: compute average and standard deviation for each sensor
    disp('=== Compute average and stdev for each sensor')

    % Extract number of channels from one epoch
    sData = in_bst_data(sFilesImport(1).FileName, 'F');
    [nChan, nTime] = size(sData.F);
    allData = zeros(nChan, nTime * length(sFilesImport)); % eg: 19 channels x (2000 times x 198 epochs)
    
    % Go through all epochs to obtain a single matrix with all data
    for i = 1:length(sFilesImport)
        sEpoch = in_bst_data(sFilesImport(i).FileName, 'F');
        
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
    jsonData.avgPerChannelPre = avgData;
    jsonData.stdPerChannelPre = stdev;

    disp('=== Finished computing average and stdev for each sensor')


    %% Quality control: generate 2D plot to visualize mean and std in each sensor
    disp('=== Show mean and std for each sensor')

    % view_topography requires a Brainstorm structure, so we use an
    % existing one but actually pass avgData and stdev to plot

    % First: mean figure
    hFigMean = view_topography(sFilesImport(1).FileName, 'EEG', '2DDisc', avgData, 1, "NewFigure");
    bst_report('Snapshot', hFigMean, sFilesImport(1).FileName, 'Mean per Channel (Pre-processing)', fig_size);
    close(hFigMean);
    
    % Second: standard deviation
    hFigStd = view_topography(sFilesImport(1).FileName, 'EEG', '2DDisc', stdev, 1, "NewFigure");
    bst_report('Snapshot', hFigStd, sFilesImport(1).FileName, 'Standard deviation per Channel (Pre-processing)', fig_size);
    close(hFigStd);

    disp('=== Finished displaying mean and std for each sensor')


    %% Quality control: check PSD pre-processing
    
    disp('=== Compute PSD prior') 
    % Process: Power spectrum density (Welch)
    sFilesPSDpre = bst_process('CallProcess', 'process_psd', sFilesEEG, [], ...
        'timewindow',  [], ...
        'win_length',  Win_length, ...
        'win_overlap', Win_overlap, ...
        'units',       'physical', ...  
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
    
    % PSD can be visualized with the snapshot or manually, the second
    % option allows changing figure size, snapshot doesn't 
    % Process: Snapshot: Frequency spectrum
    % bst_process('CallProcess', 'process_snapshot', sFilesPSDpre, [], ...
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

    % View PSD manually 
    hFigPSD = view_spectrum(sFilesPSDpre.FileName, 'Spectrum');
    set(hFigPSD, 'Position', fig_size); 
    bst_report('Snapshot', hFigPSD, sFilesPSDpre.FileName, 'PSD Spectrum', fig_size);
    close(hFigPSD);

    disp('=== Finished PSD computation')


    %% PSD peaks computation with find peaks function 
    disp('=== Compute peaks in the PSD')
    
    % Extract frequencies and PSD values
    sPSD = in_bst_timefreq(sFilesPSDpre.FileName);
    freqs = sPSD.Freqs;
    
    psdValues = sPSD.TF;      % [Nsignals x 1 x Nfreq]
    psdSqueezed = squeeze(psdValues);    % [Nsignals x Nfreq]
    
    % Average across channels 
    meanPSD = mean(psdSqueezed, 1);    % [1 x Nfreq]
    
    % Confirm this is a vector (required for findpeaks)
    if ~isvector(meanPSD)
        error('meanPSD is not a vector!');
    end

    % findpeaks requires Signal Processing Toolbox
    [peaks, peakFreqs] = findpeaks(meanPSD, freqs, 'Threshold', peaksThreshold);

    % The peak frequencies could then be used to determine frequencies to
    % filter. But only those greater than ~20Hz, otherwise we eliminate
    % alpha peak!
    
    % Save results in JSON structure
    jsonData.peakFrequencies = peakFreqs;
    jsonData.peaks = peaks; 
    
    disp('=== Finished computation of PSD peaks')
    disp('=== Display PSD peaks in Brainstorm report')

    % Create figure in the background
    hFigPeaks = figure('Visible', 'off'); 
    scatter(peakFreqs, 10*log10(peaks), 60, 'MarkerFaceColor', [0 .5 .5], 'MarkerEdgeColor', [0 .3 .3]);
    grid on;
    xlabel('Frequency (Hz)');
    xlim([0,250]);
    ylabel('Log Power (dB)');
    title(['PSD Peaks Detection - ', participant]);
    set(hFigPeaks, 'Position', fig_size); 
    
    % Make snapshot for report and close
    bst_report('Snapshot', hFigPeaks, sFilesPSDpre.FileName, 'Detected PSD Peaks Scatter', fig_size);
    close(hFigPeaks);


    disp('=== Quality Control finished')
    
    
    %% Save report and end loop
    disp('=== Save report');

    % Desired filename
    outputName = fullfile(ReportsDir, sprintf('QualityControl-%s-%s.html', participant, ProtocolName));
    
    % Save and then export to the custom name
    ReportFile = bst_report('Save', []);
    bst_report('Export', ReportFile, outputName);
    
    % Create JSON report
    jsonFile = fullfile(ReportsDir, sprintf('QC-%s-%s.json', participant, ProtocolName));
    fid = fopen(jsonFile, 'w');
    if fid == -1
        error("Cannot open JSON file.")
    end 
    fprintf(fid, '%s', jsonencode(jsonData, PrettyPrint=true));
    fclose(fid);
    
    disp(['=== Reports exported as: ', outputName, jsonFile]);


end

% Copy brainstorm_db data
copyfile([BrainstormDbDir,'/',ProtocolName], DataDir);


%% DONE
disp('=== Done!');
