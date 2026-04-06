% Import BIDS dataset to Brainstorm for further processing
% Authors: Natalia Lopez, Julia Reina, Guiomar Niso
% Cajal Institute (CSIC), Madrid, Spain 
% April, 2026 


clc; clear;

disp("=== My script has started.")

%% Set parameters

% Directory to store brainstorm database
BrainstormDbDir = fullfile(pwd, 'brainstorm_db'); 
ReportsDir = '/home/natalia/app_local/out_reports';
DataDir    = '/home/natalia/app_local/out_import';

% BIDS directory with all participants
BIDS_dir = '/home/natalia/app_local/ds003775';

% Protocol name and parameters
ProtocolName = 'Dataset003775'; 
UseDefaultAnat = 1; 
UseDefaultChannel = 0;

% Import BIDS parameters
nVertices = 15000; 


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

% Hide event files (tsv and json) to prevent warnings
disp('=== Temporarily hiding event files to suppress import warnings');
eventFiles = dir(fullfile(BIDS_dir, '**', '*_events.*'));
for i = 1:length(eventFiles)
    oldName = fullfile(eventFiles(i).folder, eventFiles(i).name);
    newName = fullfile(eventFiles(i).folder, [eventFiles(i).name, '.hidden']);
    movefile(oldName, newName); % Rename file
end

% Process: BIDS import
bst_process('CallProcess', 'process_import_bids', [], [], ...
    'bidsdir',       {BIDS_dir, 'BIDS'}, ...
    'nvertices',     nVertices, ...            % number of vertices for cortex surface reconstruction (if FS surfaces available)
    'mni',           'maff8', ...          % method for MNI normalization, if needed
    'anatregister',  'spm12', ...          % method for anatomy registration (if you want registration)
    'groupsessions', 1, ...                % whether to group multiple sessions into one subject
    'channelalign',  1);                   % align sensor positions (if head-shape or fiducials available)

% Restore all hidden files to normal
disp('=== Restoring event files');
hiddenFiles = dir(fullfile(BIDS_dir, '**', '*_events.*.hidden'));
for i = 1:length(hiddenFiles)
    oldName = fullfile(hiddenFiles(i).folder, hiddenFiles(i).name);
    newName = strrep(oldName, '.hidden', ''); % Remove the .hidden extension
    movefile(oldName, newName); % Rename back
end

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


% Copy brainstorm_db data
copyfile([BrainstormDbDir,'/',ProtocolName], DataDir);


%% DONE
disp('=== Done!');
