function optitrack_MNE_ERF(save_dir,atlas_dir, run_num, headmodel, sourcemodel)
%% Hardcoded for now
scannercast_dir = 'D:\Github\scannercast\examples\NA';

%% Load data
cd(save_dir);
load(['data_run' num2str(run_num) '.mat']);

%% Finally let's try some MNE
cfg         = [];
cfg.grad    = data.grad;   % sensor information
cfg.channel = data.label;  % the used channels
cfg.grid    = sourcemodel;   % source points
cfg.headmodel = headmodel;   % volume conduction model
cfg.singleshell.batchsize = 5000; % speeds up the computation
lf          = ft_prepare_leadfield(cfg);

% make a figure of the single subject{i} headmodel, and grid positions
figure; hold on;
ft_plot_vol(headmodel,  'facecolor', 'cortex', 'edgecolor', 'none');
alpha 0.5; camlight;
ft_plot_mesh(lf.pos(lf.inside,:),'vertexsize',1,'vertexcolor','r');
%ft_plot_sens(rawData_MEG.grad, 'style', 'r*'); view([0,0]);
ft_plot_sens(data.grad, 'style', 'r*'); view([0,0]);

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
cfg.sourcemodel         = lf;
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

%% Get M100 power
source_pow_post = get_source_pow(data,sourceall,[0.08 0.12]);
source_pow_pre  = get_source_pow(data,sourceall,[-0.08 -0.04]);

cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'pow';
sourceR = ft_math(cfg,source_pow_post,source_pow_pre);
    
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
cfg.funcolormap     = cmap;
cfg.funcolorlim     = 'maxabs';
ft_sourceplot(cfg,sourceI);

%% Export to nifti formt and use your favourite MRI software to visualise
cd(save_dir);
cfg = [];
cfg.filetype = 'nifti';
cfg.filename = ['MNE_100_run' num2str(run_num)];
cfg.parameter = 'pow';
ft_sourcewrite(cfg,sourceI);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now the VE analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%% Load atlas
% Get the path and version of Fieldtrip
[~, r] = ft_version;

atlas_HCPMMP 	= ft_read_atlas(fullfile(atlas_dir,...
    'HCP-MMP1_combined_on_spm_brain.nii'));

fid = fopen(fullfile(atlas_dir,'HCP-MMP1_combined_on_spm_brain123.txt'));
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
if run_num == 3
    figure;
    cols = [0.4275    0.9804    0.3922;0.9804    0.5686    0.3843;
        0.2824    0.3137    0.9804];
    
    for r = 1:3
        load(['MNE_VE' num2str(r) '.mat']);
        plot(sourceall.time,VE,'LineWidth',2,'Color',cols(r,:)); hold on;
        t = t+1;
    end
    
    xlim([-0.1 0.4]);
    set(gca,'FontSize',18);
    xlabel('Time (s)','FontSize',20);
    ylabel('Dipole Moment (A.U.)','FontSize',20);
    title('');
    %legend({'Sitting';'Standing';'Standing + Moving'},'Location','SouthOutside');
    print('run123_VE_ERF_MNE','-dpng','-r300');
end

end





