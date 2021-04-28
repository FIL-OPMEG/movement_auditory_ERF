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
mocap_func      = 'D:\scripts\motioncapture_functions';
atlas_dir       = 'D:\Github\analyse_OPMEG\atlas\HCPMMP'

% Add Fieldtrip to path
disp('Adding Fieldtrip and analyse_OPMEG to your MATLAB path');
addpath(fieldtripDir)
ft_defaults;

% Add analyse_OPMEG Scripts to path
addpath(genpath(script_dir));
addpath(mocap_func);
addpath(HMM_dir);

% Data and save directories:
data_dir        = 'D:\data\auditory_moving_ERF\';

% cd to save dir
cd(save_dir)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: Load and Preprocess the Optitrack & Sensor-Level OPM Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Subject 1
save_dir        = 'D:\data\auditory_moving_ERF\001_results';
motive_data     = '20210310_aef_run3_euler_world.csv'

auditoryERF_preprocess(data_dir,save_dir, 1,...
    '001',motive_data)

auditoryERF_preprocess(data_dir,save_dir, 2,...
    '001',motive_data)

auditoryERF_preprocess(data_dir,save_dir, 3,...
    '001',motive_data)

% Subject 2
save_dir        = 'D:\data\auditory_moving_ERF\002_results';
motive_data     = 'sub-RS_aef_run6_mocap_euler_world.csv'

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
save_dir        = 'D:\data\auditory_moving_ERF\001_results';
analyse_optitrack_data(save_dir)

% Subject 2
save_dir        = 'D:\data\auditory_moving_ERF\002_results';
analyse_optitrack_data(save_dir)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3: Source-Level Analysis (Beamforming)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Subject 1
save_dir        = 'D:\data\auditory_moving_ERF\001_results';
scannercast_dir = '';
optitrack_beamforming_ERF(save_dir, atlas_dir, scannercast_dir, 1)
optitrack_beamforming_ERF(save_dir, atlas_dir, scannercast_dir, 2)
optitrack_beamforming_ERF(save_dir, atlas_dir, scannercast_dir, 3)

% Subject 2
save_dir        = 'D:\data\auditory_moving_ERF\002_results';
scannercast_dir = '';
optitrack_beamforming_ERF(save_dir, atlas_dir, scannercast_dir, 1)
optitrack_beamforming_ERF(save_dir, atlas_dir, scannercast_dir, 2)
optitrack_beamforming_ERF(save_dir, atlas_dir, scannercast_dir, 3)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 4: Source-Level Analysis (MNE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Subject 1
save_dir        = 'D:\data\auditory_moving_ERF\001_results';
scannercast_dir = '';
optitrack_MNE_ERF(save_dir, atlas_dir, scannercast_dir, 1)
optitrack_MNE_ERF(save_dir, atlas_dir, scannercast_dir, 2)
optitrack_MNE_ERF(save_dir, atlas_dir, scannercast_dir, 3)

% Subject 2
save_dir        = 'D:\data\auditory_moving_ERF\002_results';
scannercast_dir = '';
optitrack_MNE_ERF(save_dir, atlas_dir, scannercast_dir, 1)
optitrack_MNE_ERF(save_dir, atlas_dir, scannercast_dir, 2)
optitrack_MNE_ERF(save_dir, atlas_dir, scannercast_dir, 3)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 5: Investigate the Beamformer By Changing Regularisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save_dir        = 'D:\data\auditory_moving_ERF\002_results';
scannercast_dir = '';
investigate_regularisation(save_dir, scannercast_dir, 3)

