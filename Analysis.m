% Analysis of resting state EEG signals in sensors and sources
% Authors: Natalia Lopez, Julia Reina, Guiomar Niso
% Cajal Institute (CSIC), Madrid, Spain 
% March, 2026 

clc; clear;

disp("=== My script has started.")

%% Set parameters

% Directory to store brainstorm database
BrainstormDbDir = fullfile(pwd, 'brainstorm_db'); 
ReportsDir = '/home/natalia/app_local/out_reports';
DataDir    = '/home/natalia/app_local/out_analysis';

% Protocol name (we use the same name as in other codes because we will
% directly import it from database)
ProtocolName = 'DatasetAD'; 
UseDefaultAnat = 1;

% PSD parameters
Win_length = 4;
Win_overlap = 50;

% Spatial smoothing parameters
FWHM = 3;

% Figure size (for report)
fig_width = 1000;
fig_height = 600;
fig_size = [200, 200, fig_width, fig_height];
fig_small = [200, 200, fig_width/2, fig_height/2];

% Frequency contact sheet parameters
FreqBands = [
   2   4;  
   5   7;   
   8  12;   
  15  29;   
  30  59;   
  60  90    
];
freqNames = {'delta', 'theta', 'alpha', 'beta', 'gamma1', 'gamma2'};
fontSize = 20;



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

% Select all continuous raw files
sFilesAll = bst_process('CallProcess', 'process_select_files_data', [], [], ...
    'subjectname',   '', ...     % empty = all subjects
    'condition',     '', ...     % empty = all conditions
    'tag',           '', ...     
    'includebad',    1, ...      
    'includeintra',  1, ...
    'includecommon', 0);

% Process: Ignore file names with tag: epoch
sFilesCont = bst_process('CallProcess', 'process_select_tag', sFilesAll, [], ...
    'tag',    'epoch', ...
    'search', 2, ...  % Search the file names
    'select', 2);  % Ignore the files with the tag

% Select only last output from preprocessing (files containing average)
sFilesFiltered = struct([]);

for i = 1:length(sFilesCont)
    if contains(sFilesCont(i).FileName, 'Average')
        sFilesFiltered = [sFilesFiltered; sFilesCont(i)];  % append struct
    end
end


% Check if we have files
if isempty(sFilesFiltered)
    error('No data files found in this protocol');
else
    disp(['Found ', num2str(length(sFilesFiltered)), ' data files in protocol.']);
end


%% Loop over participants
for iSub = 1:length(sFilesFiltered)

    participant = sFilesFiltered(iSub).SubjectName; % Extract participant information
    disp(['=== Processing participant: ', num2str(participant)]);
    
    % Select participant
    sFilesRaw = sFilesFiltered(iSub);

    % Create structure for JSON report
    jsonData = struct();
    % Basic fields
    jsonData.participant = participant;
    jsonData.protocol = ProtocolName;
    jsonData.date = datestr(now, 'yyyy-mm-dd HH:MM:SS');

    % Obtain condition for report naming (required if more than one
    % recording per subject)
    conditionName = sFilesRaw.Condition;
    
    % Then the rest of the parameters to include will be computed


    %% Analysis of resting state data on sensors
    
    disp('=== Separate in frequency bands')
    % Process: Power spectrum density (Welch), divide into frequency bands
    sFilesSensorBands = bst_process('CallProcess', 'process_psd', sFilesRaw, [], ...
        'timewindow',  [], ...      
        'win_length',  4, ...
        'win_overlap', 50, ...
        'units',       'physical', ...  % Physical: U2/Hz
        'sensortypes', '', ...
        'win_std',     0, ...
        'edit',        struct(...
             'Comment',         'Power,FreqBands', ...
             'TimeBands',       [], ...
             'Freqs',           {{'delta', '2, 4', 'mean'; 'theta', '5, 7', 'mean'; 'alpha', '8, 12', 'mean'; 'beta', '15, 29', 'mean'; 'gamma1', '30, 59', 'mean'; 'gamma2', '60, 90', 'mean'}}, ...
             'ClusterFuncTime', 'none', ...
             'Measure',         'power', ...
             'Output',          'all', ...
             'SaveKernel',      0));
    
    
    disp('=== Normalize spectrum')
    % Process: Spectrum normalization
    sFilesSensorNorm = bst_process('CallProcess', 'process_tf_norm', sFilesSensorBands, [], ...
        'normalize', 'relative2020', ...  % Relative power (divide by total power)
        'overwrite', 0);
    
    
    disp('=== Finished analysis of resting state data on sensors')
    
    
    %% Analysis of resting state data on sources
    
    % Select source files we want to process
    disp('=== Selecting source files');
    
    sFilesSources = bst_process('CallProcess', 'process_select_files_results', sFilesRaw, [], ...
        'subjectname',   participant, ...    % subject of interest
        'result',        '', ...             % empty = all source results
        'tag',           '', ...             % optionally filter on comment
        'includebad',    1, ...              % include bad results
        'includeintra',  0, ...
        'includecommon', 0);
    
    % Check that we got source files
    if isempty(sFilesSources)
        error(['No source files found for subject: ', participant]);
    else
        disp(['Found ', num2str(numel(sFilesSources)), ' source files for ', participant]);
    end
    
    disp('=== Separate in frequency bands')
    % Process: Power spectrum density (Welch), divide into frequency bands
    sFilesSourceBands = bst_process('CallProcess', 'process_psd', sFilesSources, [], ...
        'timewindow',  [], ...      
        'win_length',  4, ...
        'win_overlap', 50, ...
        'units',       'physical', ...  % Physical: U2/Hz
        'sensortypes', '', ...
        'win_std',     0, ...
        'edit',        struct(...
             'Comment',         'Power,FreqBands', ...
             'TimeBands',       [], ...
             'Freqs',           {{'delta', '2, 4', 'mean'; 'theta', '5, 7', 'mean'; 'alpha', '8, 12', 'mean'; 'beta', '15, 29', 'mean'; 'gamma1', '30, 59', 'mean'; 'gamma2', '60, 90', 'mean'}}, ...
             'ClusterFuncTime', 'none', ...
             'Measure',         'power', ...
             'Output',          'all', ...
             'SaveKernel',      0));
    
    
    disp('=== Normalize spectrum')
    % Process: Spectrum normalization
    sFilesSourceNorm = bst_process('CallProcess', 'process_tf_norm', sFilesSourceBands, [], ...
        'normalize', 'relative2020', ...  % Relative power (divide by total power)
        'overwrite', 0);
    
    disp('=== Perform spatial smoothing')
    % Process: Spatial smoothing (3.00 mm)
    sFilesSourceSmooth = bst_process('CallProcess', 'process_ssmooth_surfstat', sFilesSourceNorm, [], ...
        'fwhm',      FWHM, ...
        'method',    'fixed_fwhm', ...  % Fixed FWHM for all surfaces
        'overwrite', 0);
    
    
    %% Save snapshot of frequency contact sheet
    disp('=== Save frequency contact sheet')
    
    % Parameters for contact sheet
    nFreqs = length(FreqBands);
    centerFreqBands = mean(FreqBands, 2);

    % Load image and set label size and background color
    hFig = view_surface_data([], sFilesSourceSmooth.FileName, [], 'NewFigure');
    set(hFig, 'Color', [1 1 1]);
    set(hFig, 'Position', fig_small);

    % Create contact sheet and set parameters
    hFigContact = view_contactsheet(hFig, 'freq', 'fig', [], nFreqs, centerFreqBands(1:nFreqs));
    set(hFigContact, 'Color', [1 1 1]);
    hContactAxes = findobj(hFigContact, 'Type', 'axes');
    set(hContactAxes, 'FontSize', fontSize, 'FontWeight', 'bold');

    % Save image
    bst_report('Snapshot', hFigContact, sFilesSourceSmooth.FileName, 'Sources Frequency Contact Sheet', fig_size);
    close(hFig);
    close(hFigContact);
    

    %% Project to default template (if individual anatomy)
    
    % Process: project on default anatomy (surface)
    if UseDefaultAnat == 0
        bst_process('CallProcess', 'process_project_sources', sFilesSourceSmooth, [], ...
             'headmodeltype', 'surface');
    end


    %% Scouts 
    % Process: Scout time series: [68 scouts]
    sFilesScouts = bst_process('CallProcess', 'process_extract_scout', sFilesSourceSmooth, [], ...
        'timewindow',     [], ...
        'scouts',         {'Desikan-Killiany', {'bankssts L', 'bankssts R', 'caudalanteriorcingulate L', 'caudalanteriorcingulate R', 'caudalmiddlefrontal L', 'caudalmiddlefrontal R', 'cuneus L', 'cuneus R', 'entorhinal L', 'entorhinal R', 'frontalpole L', 'frontalpole R', 'fusiform L', 'fusiform R', 'inferiorparietal L', 'inferiorparietal R', 'inferiortemporal L', 'inferiortemporal R', 'insula L', 'insula R', 'isthmuscingulate L', 'isthmuscingulate R', 'lateraloccipital L', 'lateraloccipital R', 'lateralorbitofrontal L', 'lateralorbitofrontal R', 'lingual L', 'lingual R', 'medialorbitofrontal L', 'medialorbitofrontal R', 'middletemporal L', 'middletemporal R', 'paracentral L', 'paracentral R', 'parahippocampal L', 'parahippocampal R', 'parsopercularis L', 'parsopercularis R', 'parsorbitalis L', 'parsorbitalis R', 'parstriangularis L', 'parstriangularis R', 'pericalcarine L', 'pericalcarine R', 'postcentral L', 'postcentral R', 'posteriorcingulate L', 'posteriorcingulate R', 'precentral L', 'precentral R', 'precuneus L', 'precuneus R', 'rostralanteriorcingulate L', 'rostralanteriorcingulate R', 'rostralmiddlefrontal L', 'rostralmiddlefrontal R', 'superiorfrontal L', 'superiorfrontal R', 'superiorparietal L', 'superiorparietal R', 'superiortemporal L', 'superiortemporal R', 'supramarginal L', 'supramarginal R', 'temporalpole L', 'temporalpole R', 'transversetemporal L', 'transversetemporal R'}}, ...
        'scoutfunc',      'mean', ...  % Mean
        'pcaedit',        struct(...
             'Method',         'pca', ...
             'Baseline',       [-0.1, 0], ...
             'DataTimeWindow', [0, 1], ...
             'RemoveDcOffset', 'file'), ...
        'isnorm',         0, ...
        'concatenate',    0, ...
        'save',           1, ...
        'addrowcomment',  1, ...
        'addfilecomment', []);


    disp('=== Finished analysis of resting state on sources')
    
    
    %% Save report and end loop
    
    % Save report
    disp('=== Save report');
    % Desired filename
    outputName = fullfile(ReportsDir, sprintf('Analysis-%s-%s-%s.html', participant, conditionName, ProtocolName));
    
    % Save and then export to the custom name
    ReportFile = bst_report('Save', []);
    bst_report('Export', ReportFile, outputName);
    
    % Create JSON report
    jsonFile = fullfile(ReportsDir, sprintf('Analysis-%s-%s-%s.json', participant, conditionName, ProtocolName));
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


