%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pipeline for the Auditory ERF data during various levels of subject
% movement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start-up: Add the appropriate toolboxes and scripts to your path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fieldtripDir    = 'D:\scripts\fieldtrip-master';
script_dir      = 'D:\Github\analyse_OPMEG';
mocap_func      = 'D:\Github\optitrack';
atlas_dir       = 'D:\Github\analyse_OPMEG\atlas\HCPMMP'

% Add Fieldtrip to path
disp('Adding Fieldtrip and analyse_OPMEG to your MATLAB path');
addpath(fieldtripDir)
ft_defaults;

% Add analyse_OPMEG Scripts to path
addpath(genpath(script_dir));
addpath(mocap_func);

% BIDS data directory:
data_dir        = 'D:\data\auditory_moving_ERF_BIDS\';
cd(data_dir);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: Load and Preprocess the Optitrack & Sensor-Level OPM Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Subject 1
save_dir        = fullfile(data_dir,'results','001_results');
motive_data     = fullfile(data_dir,'sub-001','ses-001','meg','motion',...
    'sub-001_ses-001_task-aef_run-003_eul_world.csv');

auditoryERF_preprocess(data_dir,save_dir, 1,...
    '001',motive_data)

auditoryERF_preprocess(data_dir,save_dir, 2,...
    '001',motive_data)

auditoryERF_preprocess(data_dir,save_dir, 3,...
    '001',motive_data)

% Subject 2
save_dir        = 'D:\data\auditory_moving_ERF_BIDS\results\002_results';
motive_data     = fullfile(data_dir,'sub-002','ses-001','meg','motion',...
    'sub-002_ses-001_task-aef_run-003_eul_world.csv');

auditoryERF_preprocess(data_dir,save_dir, 1,...
    '002',motive_data)

auditoryERF_preprocess(data_dir,save_dir, 2,...
    '002',motive_data)

auditoryERF_preprocess(data_dir,save_dir, 3,...
    '002',motive_data)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: Further process and plot the opti-track data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Subject 1
save_dir        = fullfile(data_dir,'results','001_results');
analyse_optitrack_data(save_dir)

% Subject 2
save_dir        = fullfile(data_dir,'results','002_results');
analyse_optitrack_data(save_dir)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3: Source-Level Analysis (Beamforming)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Subject 1
save_dir        = fullfile(data_dir,'results','001_results');

% Load pre-computed headmodel and sourcemodel
scannercast_dir = fullfile(data_dir,'fieldtrip_sourcespace','sub-001')
load(fullfile(scannercast_dir,'sub-001_desc-headmodel.mat'));
load(fullfile(scannercast_dir,'sub-001_desc-sourcemodel_5mm.mat'));

optitrack_beamforming_ERF(save_dir, atlas_dir, 1, headmodel, sourcemodel)
optitrack_beamforming_ERF(save_dir, atlas_dir, 2, headmodel, sourcemodel)
optitrack_beamforming_ERF(save_dir, atlas_dir, 3, headmodel, sourcemodel)

% Subject 2
save_dir        = fullfile(data_dir,'results','002_results');

% Load pre-computed headmodel and sourcemodel
scannercast_dir = fullfile(data_dir,'fieldtrip_sourcespace','sub-002')
load(fullfile(scannercast_dir,'sub-002_desc-headmodel.mat'));
load(fullfile(scannercast_dir,'sub-002_desc-sourcemodel_5mm.mat'));

optitrack_beamforming_ERF(save_dir, atlas_dir, 1, headmodel, sourcemodel)
optitrack_beamforming_ERF(save_dir, atlas_dir, 2, headmodel, sourcemodel)
optitrack_beamforming_ERF(save_dir, atlas_dir, 3, headmodel, sourcemodel)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 4: Source-Level Analysis (MNE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Subject 1
save_dir        = fullfile(data_dir,'results','001_results');

% Load pre-computed headmodel and sourcemodel
scannercast_dir = fullfile(data_dir,'fieldtrip_sourcespace','sub-001')
load(fullfile(scannercast_dir,'sub-001_desc-headmodel.mat'));
load(fullfile(scannercast_dir,'sub-001_desc-sourcemodel_5mm.mat'));

optitrack_MNE_ERF(save_dir, atlas_dir, 1, headmodel, sourcemodel)
optitrack_MNE_ERF(save_dir, atlas_dir, 2, headmodel, sourcemodel)
optitrack_MNE_ERF(save_dir, atlas_dir, 3, headmodel, sourcemodel)

% Subject 2
save_dir        = fullfile(data_dir,'results','002_results');

scannercast_dir = fullfile(data_dir,'fieldtrip_sourcespace','sub-002')
load(fullfile(scannercast_dir,'sub-002_desc-headmodel.mat'));
load(fullfile(scannercast_dir,'sub-002_desc-sourcemodel_5mm.mat'));

optitrack_MNE_ERF(save_dir, atlas_dir, 1, headmodel, sourcemodel)
optitrack_MNE_ERF(save_dir, atlas_dir, 2, headmodel, sourcemodel)
optitrack_MNE_ERF(save_dir, atlas_dir, 3, headmodel, sourcemodel)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 5: Investigate the Beamformer By Changing Regularisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save_dir        = fullfile(data_dir,'results','002_results');
investigate_regularisation(save_dir, headmodel,sourcemodel)
