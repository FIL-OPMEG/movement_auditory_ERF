function PSD_after_beamformer(save_dir, headmodel,sourcemodel,mri)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create VE for continuous un-epoched data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load data
cd(save_dir);

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

for run_num = 1:3
    d = load(['data_out_si_hp' num2str(run_num) '.mat']);
%% Low Pass Filter
cfg                 = [];
cfg.lpfilter        = 'yes';
cfg.lpfreq          = 40;
data_out_si_lp_hp   = ft_preprocessing(cfg,d.data_out_si_hp);

%% Compute covariance matrix
cfg                  = [];
cfg.covariance       = 'yes';
cfg.vartrllength     = 2;
cfg.covariancewindow = 'all';
avg                  = ft_timelockanalysis(cfg,data_out_si_lp_hp);

%% Source Analysis
cfg                    = [];
cfg.channel            = data.label;
cfg.grad               = data.grad;
cfg.method             = 'lcmv';
cfg.grid               = lf_2;
cfg.headmodel          = headmodel;
cfg.lcmv.keepfilter    = 'yes';
cfg.lcmv.fixedori      = 'yes';
cfg.lcmv.lambda        = '0.1%';
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

save(['VE_' num2str(run_num) '.mat'],'VE');
clear data_out_si_lp_hp

end

end




