function investigate_regularisation(save_dir, headmodel,sourcemodel)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Investigate various levels of beamformer regularisation on whole-brain
% M100 maps and source-level PSD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load data
cd(save_dir);

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
% 
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
cfg.covariancewindow = [0 0.5];
avg                  = ft_timelockanalysis(cfg,data);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START OF REGULARISATION LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

reg = {'0.01%','0.1%','1%','5%','10%','100%','1000%'};

for rr = 1:length(reg)
    
    cfg                    = [];
    cfg.channel            = data.label;
    cfg.grad               = data.grad;
    cfg.method             = 'lcmv';
    cfg.grid               = lf;
    cfg.headmodel          = headmodel;
    cfg.lcmv.keepfilter    = 'yes';
    cfg.lcmv.fixedori      = 'yes';
    cfg.lcmv.projectnoise  = 'yes';
    cfg.lcmv.weightnorm    = 'nai';
    cfg.lcmv.lambda        = reg{rr};
    sourceall              = ft_sourceanalysis(cfg, avg);
    
    %% Replace .pos field with template_grid.pos
    
    [t, r] = ft_version;
    ddd = load(fullfile(r,'template/sourcemodel/standard_sourcemodel3d5mm.mat'));
    template_grid = ddd.sourcemodel;
    clear ddd
    template_grid = ft_convert_units(template_grid,'mm');
    
    sourceall.pos = template_grid.pos;
    %%
    % Remove cfg field to save memory
    sourceall = rmfield(sourceall,'cfg');
        
    source_pow_post = get_source_pow(data,sourceall,[0.08 0.12]);
    source_pow_pre  = get_source_pow(data,sourceall,[-0.08 -0.04]);
    
    %
    cfg = [];
    cfg.operation = 'subtract';
    cfg.parameter = 'pow';
    sourceR = ft_math(cfg,source_pow_post,source_pow_pre);

    %% Interpolate
    spm_brain = ft_read_mri('D:\scripts\fieldtrip-master\template\anatomy\single_subj_T1.nii');
    cfg              = [];
    cfg.voxelcoord   = 'no';
    cfg.parameter    = 'pow';
    cfg.interpmethod = 'nearest';
    sourceI  = ft_sourceinterpolate(cfg, sourceR, spm_brain);
    
    %%
    % Change the colormap to RdBu
    ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
    cmap = colormap(flipud(brewermap(64,'RdBu'))); % change the colormap
    
    %% Export to nifti formt and use your favourite MRI software to visualise
    cd(save_dir);
    cfg = [];
    cfg.filetype = 'nifti';
    cfg.filename = ['reg_' reg{rr}];
    cfg.parameter = 'pow';
    ft_sourcewrite(cfg,sourceI);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% From run 3: continuous, preprocessed data
load('data_out_si_lp_hp.mat');

cfg                  = [];
cfg.covariance       = 'yes';
cfg.vartrllength     = 2;
cfg.covariancewindow = 'all';
avg                  = ft_timelockanalysis(cfg,data_out_si_lp_hp);

%% Load Source stuff
%  LF for Virtual electrode analysis
%  For simplicity let's just use a single point in auditory cortex
cfg             = [];
cfg.template    = mri;
cfg.nonlinear   = 'yes';
norm            = ft_volumenormalise([],mri);

% Auditory Cortex
pos = [-48 -22 4; 48 -22 4];

% Now we warp the MNI coordinates using the nonlinear warping parameters
posback         = ft_warp_apply(norm.params,pos,'sn2individual');
% xyz positions in individual coordinates
pos_grid        = ft_warp_apply(pinv(norm.initial),posback);

figure; ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:));
ft_plot_mesh(pos_grid,'vertexcolor','r');

%%
% Prepare Leadfield
cfg                 = [];
cfg.method          = 'lcmv';
cfg.channel         = data.label;
cfg.grid.pos        = pos_grid;
cfg.grid.unit       = 'mm';
cfg.headmodel       = headmodel;
cfg.grad            = data.grad;
cfg.reducerank      = 2; %(default = 3 for EEG, 2 for MEG)
cfg.normalize       = 'yes' ; %Normalise Leadfield: 'yes' for beamformer
cfg.normalizeparam  = 1;
lf_2 = ft_prepare_leadfield(cfg);

% % Concat the leadfields
lf_concat = cat(2,lf_2.leadfield{:});
for k = 1:2
    lf.leadfield{k} = lf_concat;
end

%% Regulsarisation Loop
reg = {'0.01%','0.1%','1%','5%','10%','100%','1000%'};

psd_all = [];

for rr = 1:length(reg)
    cfg                    = [];
    cfg.channel            = data.label;
    cfg.grad               = data.grad;
    cfg.method             = 'lcmv';
    cfg.grid               = lf_2;
    cfg.headmodel          = headmodel;
    cfg.lcmv.keepfilter    = 'yes';
    cfg.lcmv.fixedori      = 'yes';
    cfg.lcmv.lambda        = reg{rr};
    sourceall              = ft_sourceanalysis(cfg, avg);
    
% Find filter from Idx point
filter123 = cat(1,sourceall.avg.filter{1,:});

VE          = [];
VE.label    = {'A1'};
VE.fsample  = data_out_si_lp_hp.fsample;
for subs=1:size(data_out_si_lp_hp.trial,2)
    % note that this is the non-filtered "raw" data
    VE.time{subs}       = data_out_si_lp_hp.time{subs};
    VE.trial{subs}(:,:) = filter123(:,:)*data_out_si_lp_hp.trial{subs}(:,:);
end
   
    Fs = VE.fsample;
    L = 1000;
    noverlap = 100;
    [pxx,f,pxxc] = pwelch(VE.trial{1},hamming(L),noverlap,2000,Fs);

    % Add outside the loop
    psd_all(rr,:) = squeeze(pxx);
end

%%
cfg                 = [];
cfg.channel         = 'all'
cfg.trial_length    = 3;
cfg.method          = 'tim';
cfg.foi             = [2 40];
cfg.plot            = 'yes';
cfg.plot_legend     = 'no';
cfg.transparency    = 1;
cfg.plot_ci         = 'no';
cfg.plot_legend     = 'yes';
cfg.plot_mean       = 'no';
cfg.LineWidth = 2;
plot_powspctrm(fliplr(psd_all'),cfg,flipud(reg'),f);
ylim([10e-18 3e-15]); 
print('opti_VE_PSD','-dpng','-r300');

end




