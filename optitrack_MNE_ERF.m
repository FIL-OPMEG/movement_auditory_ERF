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
run_num = 2;

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
cfg.latency          = 'all';
cfg.covariancewindow = [0 0.5];
avg                  = ft_timelockanalysis(cfg,data);

%% Let's MNE
cfg                     = [];
cfg.method              = 'mne';
cfg.sourcemodel         = leadfield;
cfg.headmodel           = headmodel;
cfg.mne.prewhiten       = 'yes';
cfg.mne.lambda          = 3;
cfg.mne.scalesourcecov  = 'yes';
sourceall               = ft_sourceanalysis(cfg,avg);

%% Project the dipole moment 
cfg                 = [];
cfg.projectmom      = 'yes';
sourceall           = ft_sourcedescriptives(cfg,sourceall);

%% Replace .pos field with template grid .pos
[t, r] = ft_version;

load(fullfile(r,'template/sourcemodel/standard_sourcemodel3d5mm.mat'));
template_grid = sourcemodel;
clear sourcemodel

sourceall.pos = template_grid.pos;

%%
addpath(genpath('D:\Github\MQ_MEG_Scripts'));
sourceR = get_source_pow(data,sourceall,[0.15 0.25]);

%% Interpolate
spm_brain = ft_read_mri('D:\scripts\fieldtrip-master\template\anatomy\single_subj_T1.nii');       
cfg              = [];
cfg.voxelcoord   = 'no';
cfg.parameter    = 'avg.pow';
cfg.interpmethod = 'nearest';
sourceI  = ft_sourceinterpolate(cfg, sourceR,spm_brain);

%%
% Change the colormap to RdBu
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
cmap = colormap(flipud(brewermap(64,'RdBu'))); % change the colormap

% Plot
cfg                 = [];
cfg.funparameter    = 'pow';
cfg.funcolormap        = cmap;
cfg.funcolorlim     = 'maxabs';
%cfg.maskparameter   = 'anat_mask';
ft_sourceplot(cfg,sourceI);

%% Export to nifti formt and use your favourite MRI software to visualise
cd(save_dir);
cfg = [];
cfg.filetype = 'nifti';
cfg.filename = ['MNE_M200_run' num2str(run_num)];
cfg.parameter = 'pow';
ft_sourcewrite(cfg,sourceI);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now the VE analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%% Load atlas
% Get the path and version of Fieldtrip
[~, r] = ft_version;

atlas_HCPMMP 	= ft_read_atlas(fullfile(HMM_dir,'HCPMMP',...
    'HCP-MMP1_combined_on_spm_brain.nii'));

fid = fopen(fullfile(HMM_dir,'HCPMMP','HCP-MMP1_combined_on_spm_brain123.txt'));
labels = textscan(fid,'%s');
fclose(fid);

atlas_HCPMMP.tissuelabel = labels{1};
atlas_HCPMMP.tissue            = atlas_HCPMMP.parcellation
atlas_HCPMMP                   = rmfield(atlas_HCPMMP,'parcellation');
atlas_HCPMMP.coordsys          = 'mni';

%%
cfg = [];
cfg.interpmethod = 'nearest';
cfg.parameter = 'tissue';
sourcemodel2 = ft_sourceinterpolate(cfg, atlas_HCPMMP, template_grid);

% Find atlas points
atlas_points = find(sourcemodel2.tissue==10);

VE = [];
for a = 1:length(atlas_points)
    VE = vertcat(VE, sourceall.avg.mom{atlas_points(a)});
end

VE = mean(VE,1);
figure; plot(sourceall.time,VE);
title(['MNE_' num2str(run_num)]);

% Save VE
save(['MNE_VE' num2str(run_num)],'VE');

%% Plot all on the same graph
figure;
cols = [0.4275    0.9804    0.3922;0.9804    0.5686    0.3843;
    0.2824    0.3137    0.9804];

for run_num = 1:3
    load(['MNE_VE' num2str(run_num) '.mat']);
    
    plot(sourceall.time,VE,'LineWidth',2,'Color',cols(run_num,:)); hold on;

end
xlim([-0.1 0.4]);
set(gca,'FontSize',18);
xlabel('Time (s)','FontSize',20);
ylabel('???','FontSize',20);
title('');
%legend({'Sitting';'Standing';'Standing + Moving'},'Location','SouthOutside');
print('run123_VE_ERF_MNE','-dpng','-r300');

%% Plot
figure;
cfg = [];
cfg.channel = VE_A1.label;
cfg.parameter = 'avg';
cfg.baseline = [-0.2 0];
cfg.showlegend    = 'no';
cfg.xlim    = [-0.1 0.4];
cfg.linecolor = [0.4275    0.9804    0.3922;0.9804    0.5686    0.3843;
    0.2824    0.3137    0.9804];
cfg.linewidth = 2;
cfg.ylim = [-2.8e-7 -2.8e-7];
figure; ft_singleplotER(cfg,avg_VE{1},avg_VE{2},avg_VE{3})
set(gca,'FontSize',18);
xlabel('Time (s)','FontSize',20);
ylabel('???','FontSize',20);
title('');
print('run123_VE_ERF','-dpng','-r300');
legend({'Sitting';'Standing';'Standing + Moving'},'Location','SouthOutside');





