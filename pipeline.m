%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis pipeline for the article:
%
% 'Using OPMs to measure neural activity in standing, mobile participants'
% Seymour et al., (2021)
%
% MATLAB and Python scripts were written by 
% Dr. Robert Seymour, February 2021-May 2021
% For enquiries, please contact: rob.seymour@ucl.ac.uk
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start-up: Add the appropriate toolboxes and scripts to your path
%           These paths are specific to my PC - change accordingly.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fieldtripDir    = 'D:\scripts\fieldtrip-master';
script_dir      = 'D:\Github\analyse_OPMEG';
mocap_func      = 'D:\Github\optitrack';
atlas_dir       = 'D:\Github\analyse_OPMEG\atlas\HCPMMP';

% Add Fieldtrip to path
disp('Adding Fieldtrip and analyse_OPMEG to your MATLAB path');
addpath(fieldtripDir)
ft_defaults;

% Add analyse_OPMEG Scripts to path
addpath(genpath(script_dir));
addpath(mocap_func);

% BIDS data directory. This is specific to my PC - change accordingly.
data_dir        = 'D:\data\auditory_moving_ERF_BIDS\';
cd(data_dir);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: Load and Preprocess the Motion Capture & Sensor-Level OPM Data
%         This corresponds to Figure 4 in the article
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

% This function plots sensor-level t-stats from the 
% OPM channel with the best M100 response (for Supplementary Figure 4)
plot_sensor_level_ERF_tstat(save_dir,'N3-TAN')

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

% This function plots sensor-level t-stats from the 
% OPM channel with the best M100 response (for Supplementary Figure 4)
plot_sensor_level_ERF_tstat(save_dir,'1A-TAN')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: Further process and plot the motion capture data
%         This corresponds to Figures 2-3 in the article
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
%         Whole-brain maps (Figure 5) were plotted with FSL
%         Figure 6 is made using the optitrack_beamforming_ERF.m function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Subject 1
save_dir        = fullfile(data_dir,'results','001_results');

% Load pre-computed headmodel and sourcemodel
mri             = ft_read_mri(fullfile(data_dir,'sub-001','ses-001',...
    'anat','001.nii')); mri.coordsys = 'neuromag';
scannercast_dir = fullfile(data_dir,'fieldtrip_sourcespace','sub-001')
load(fullfile(scannercast_dir,'sub-001_desc-headmodel.mat'));
load(fullfile(scannercast_dir,'sub-001_desc-sourcemodel_5mm.mat'));

optitrack_beamforming_ERF(save_dir, atlas_dir, 1, headmodel, sourcemodel)
optitrack_beamforming_ERF(save_dir, atlas_dir, 2, headmodel, sourcemodel)
optitrack_beamforming_ERF(save_dir, atlas_dir, 3, headmodel, sourcemodel)

% This function will create VE of continuous data rather than epoched data 
PSD_after_beamformer(save_dir, headmodel,sourcemodel,mri)

% Subject 2
save_dir        = fullfile(data_dir,'results','002_results');

% Load MRI pre-computed headmodel and sourcemodel
mri             = ft_read_mri(fullfile(data_dir,'sub-002','ses-001',...
    'anat','002.nii')); mri.coordsys = 'neuromag';
scannercast_dir = fullfile(data_dir,'fieldtrip_sourcespace','sub-002')
load(fullfile(scannercast_dir,'sub-002_desc-headmodel.mat'));
load(fullfile(scannercast_dir,'sub-002_desc-sourcemodel_5mm.mat'));

optitrack_beamforming_ERF(save_dir, atlas_dir, 1, mri, headmodel, sourcemodel)
optitrack_beamforming_ERF(save_dir, atlas_dir, 2, mri, headmodel, sourcemodel)
optitrack_beamforming_ERF(save_dir, atlas_dir, 3, mri, headmodel, sourcemodel)

% This function will create VE of continuous data rather than epoched data 
PSD_after_beamformer(save_dir, headmodel,sourcemodel,mri);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 4: Plot PSD of raw data, preprocessed data and beamformed data
%         This corresponds to Figure 7 in the manuscript
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_psd_preprocessing(fullfile(data_dir,'results','001_results'),...
    fullfile(data_dir,'results','002_results'))

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 5: Investigate the Beamformer By Changing Regularisation
%         This corresponds to Supplementary Figures S6-7              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mri             = ft_read_mri(fullfile(data_dir,'sub-002','ses-001',...
    'anat','002.nii')); mri.coordsys = 'neuromag';
scannercast_dir = fullfile(data_dir,'fieldtrip_sourcespace','sub-002')
load(fullfile(scannercast_dir,'sub-002_desc-headmodel.mat'));
load(fullfile(scannercast_dir,'sub-002_desc-sourcemodel_5mm.mat'));

save_dir        = fullfile(data_dir,'results','002_results');

investigate_regularisation(save_dir, headmodel,sourcemodel)
