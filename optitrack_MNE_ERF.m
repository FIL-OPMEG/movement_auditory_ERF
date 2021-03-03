%% Paths (RS)
fieldtripDir    = 'D:\scripts\fieldtrip-master';
script_dir      = 'D:\Github\analyse_OPMEG';
data_dir        = 'D:\data\20201208_optitrack';
save_dir        = 'D:\data\20201208_optitrack';
mocap_func      = 'D:\scripts\motioncapture_functions';
HMM_dir         = 'D:\Github\HMM'

% Add Fieldtrip to path
disp('Adding Fieldtrip and analyse_OPMEG to your MATLAB path');
addpath(fieldtripDir)
ft_defaults;

% Add analyse_OPMEG Scripts to path
addpath(genpath(script_dir));
addpath(mocap_func);
addpath(HMM_dir);

% cd to save dir
cd(save_dir)

%%
run_num = 3;

%% Load data
load(['data_run' num2str(run_num) '.mat']);

%% Whole-brain 
% Prepare leadfield
cd('D:\Github\scannercast\examples\NA');
load('headmodel.mat');
clear sourcemodel
load('sourcemodel_5mm.mat');
mri = ft_read_mri('NA.nii');
mri.coordsys = 'neuromag';

% Prepare Leadfield
cfg                 = [];
cfg.method          = 'lcmv';
cfg.channel         = data.label;
cfg.grid            = sourcemodel;
cfg.grid.unit       = 'mm';
cfg.headmodel       = headmodel;
cfg.grad            = data.grad;
cfg.reducerank      = 2%(default = 3 for EEG, 2 for MEG)
cfg.normalize       = 'yes' ; %Normalise Leadfield: 'yes' for beamformer
cfg.normalizeparam  = 1;
lf                  = ft_prepare_leadfield(cfg);

%% Make leadfields symmetric across hemispheres
lf1 = reshape(lf.leadfield, lf.dim);
lf2 = flip(lf1,1);
for k = 1:numel(lf1)
    if ~isempty(lf1{k})&&~isempty(lf2{k})
        lf.leadfield{k} = [lf1{k} lf2{k}];
    else
        lf.leadfield{k} = [];
        lf.inside(k)    = false;
    end
end
clear lf1 lf2

% make a figure of the single subject{i} headmodel, and grid positions
figure; hold on;
ft_plot_vol(headmodel,  'facecolor', 'cortex', 'edgecolor', 'none');
alpha 0.5; camlight;
ft_plot_mesh(lf.pos(lf.inside,:),'vertexsize',1,'vertexcolor','r');
%ft_plot_sens(rawData_MEG.grad, 'style', 'r*'); view([0,0]);
ft_plot_sens(data.grad, 'style', 'r*'); view([0,0]);


%% Finally let's try some MNE
cfg         = [];
cfg.grad    = data.grad;   % sensor information
cfg.channel = data.label;  % the used channels
cfg.grid    = sourcemodel;   % source points
cfg.headmodel = headmodel;   % volume conduction model
cfg.singleshell.batchsize = 5000; % speeds up the computation
leadfield   = ft_prepare_leadfield(cfg);

%% Compute covariance matrix
cfg                  = [];
cfg.covariance       = 'yes';
cfg.vartrllength     = 2;
cfg.latency          = [0.08 0.12];
cfg.covariancewindow = [0 0.5];
avg                  = ft_timelockanalysis(cfg,data);
cfg.latency          = [-0.05 -0.01];
avg2                  = ft_timelockanalysis(cfg,data);

cfg               = [];
cfg.method        = 'mne';
cfg.grid          = leadfield;
cfg.headmodel     = headmodel;
cfg.mne.prewhiten = 'yes';
cfg.mne.lambda    = 3;
%cfg.mne.scalesourcecov = 'yes';
source1            = ft_sourceanalysis(cfg,avg);
source2            = ft_sourceanalysis(cfg,avg2);

cfg = [];
cfg.parameter = 'avg.pow'
cfg.operation = 'subtract';
sourceR = ft_math(cfg,source1,source2);

cfg = [];
%cfg.latency = 'all';
cfg.funparameter = 'pow';
figure;ft_sourceplot(cfg,sourceR);







