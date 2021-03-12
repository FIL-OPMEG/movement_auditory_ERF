function auditoryERF_preprocess(data_dir,save_dir, run_num, motive_data

%% Hard coded variables (for now)
session_number = '002';
task_name      = 'aef';
subject_name   = 'NA';

%% Start preprocessing.
% Read in the raw data using BIDS
cfg             = [];
cfg.folder      = data_dir;
cfg.precision   = 'single';
cfg.bids.task   = task_name;
cfg.bids.sub    = subject_name;
cfg.bids.ses    = session_number;
if run_num == 1
    cfg.bids.run    = '001';
elseif run_num == 2
    cfg.bids.run    = '002';
elseif run_num == 3
    cfg.bids.run    = '003';
end
rawData         = ft_opm_create(cfg);

%% Resample to 1000Hz
cfg                 = [];
cfg.resamplefs      = 1000;
[rawData]           = ft_resampledata(cfg, rawData);

% % Plot using ft_databrowser
cfg             = [];
cfg.blocksize   = 30;
%cfg.event       = banana;
cfg.channel     = vertcat(ft_channelselection_opm('MEG',rawData));
cfg.viewmode    = 'butterfly';
cfg.colorgroups   = 'allblack';
ft_databrowser(cfg,rawData);

%% Plot PSD
cfg                 = [];
cfg.channel         = vertcat(ft_channelselection_opm('MEG',rawData));
cfg.trial_length    = 10;
cfg.method          = 'tim';
cfg.foi             = [1 100];
cfg.plot            = 'yes';
[pow freq]          = ft_opm_psd(cfg,rawData);
ylim([1 1e4])
title('Raw Data');

%% Load in optitrack data
if run_num == 3
    opti_data = csv2mat_sm(motive_data);
    
    % Plot the rotations
    plot_motive_rotation(opti_data)
    
    % Plot the translations
    plot_motive_translation(opti_data)
    
    % Plot mean marker error
    figure;
    plot(opti_data.time,opti_data.rigidbodies.data(:,8),'LineWidth',2);
    ylabel('Mean Marker Error');xlabel('Time (s)');
    
    %% Sync up opti-track and rawData
    [MovementDataOut, OPMdataOut] = syncOptitrackAndOPMdata(opti_data,...
        rawData,'TriggerChannelName','FluxZ-A');
    save(['MovementDataOut_run' num2str(run_num)], 'MovementDataOut')
    
end

%% Select the OPM Data
cfg             = [];
cfg.channel     = vertcat(ft_channelselection_opm('MEG',rawData));
%rawData_MEG     = ft_selectdata(cfg,OPMdataOut);
rawData_MEG     = ft_selectdata(cfg,rawData);

%% Regress the optitrack data from the
ref                 = MovementDataOut.rigidbodies.data(:,1:6);
[rawData_MEG_reg]   = regress_motive_OPMdata(rawData_MEG,ref,10);
%rawData_MEG_reg     = rawData_MEG;


%% Plot data
cfg             = [];
cfg.blocksize   = 30;
%cfg.channel     = vertcat(ft_channelselection_opm('MEG',rawData));
cfg.viewmode    = 'butterfly';
cfg.colorgroups = 'allblack';
ft_databrowser(cfg,rawData_MEG_reg);

%% MFC
[data_out_mfc, M, chan_inds] = ft_denoise_mfc(rawData_MEG_reg);

%% Plot PSD
cfg                 = [];
cfg.channel         = vertcat(ft_channelselection_opm('MEG',rawData));
cfg.trial_length    = 10;
cfg.method          = 'tim';
cfg.foi             = [0 150];
cfg.plot            = 'yes';
cfg.plot_legend     = 'no';
[pow freq]          = ft_opm_psd(cfg,data_out_mfc);
ylim([1 1e4])

%% Spectral Interpolation
cfg                     = [];
cfg.channel             = 'all';
cfg.dftfilter           = 'yes';
cfg.dftfreq             = [50 83 100 120 150];
cfg.dftreplace          = 'neighbour';
cfg.dftbandwidth        = [2 2 2 3 2];
cfg.dftneighbourwidth   = [1 2 2 2 2];
data_out_si             = ft_preprocessing(cfg,data_out_mfc);

%% Low Pass Filter
cfg                 = [];
cfg.lpfilter        = 'yes';
cfg.lpfreq          = 40;
data_out_si_lp      = ft_preprocessing(cfg,data_out_si);

%% HP-filter
cfg                 = [];
cfg.hpfilter        = 'yes';
cfg.hpfreq          = 2;
cfg.filtord         = 5;
cfg.hpinstabilityfix = 'reduce';
%cfg.hpfilttype = 'fir';
data_out_si_lp_hp     = ft_preprocessing(cfg,data_out_si_lp);

%% Remove DS - channel is bad
cfg                 = [];
cfg.channel         = vertcat(data_out_si_lp.label,'-DS-TAN','-DS-RAD');
data_out_si_lp_hp   = ft_selectdata(cfg,data_out_si_lp_hp);

%% Plot data for Artifact Rejection
cfg             = [];
cfg.blocksize   = 10;
%cfg.channel     = vertcat(ft_channelselection_opm('MEG',rawData));
cfg.viewmode    = 'butterfly';
cfg.colorgroups = 'allblack';
arft            = ft_databrowser(cfg,data_out_si_lp_hp);


%% Turn the highlighted data into nans
arft.artfctdef.reject          = 'nan';
data_out_si_lp_hp_arft = ft_rejectartifact(arft, data_out_si_lp_hp);

cfg             = [];
cfg.blocksize   = 10;
%cfg.channel     = vertcat(ft_channelselection_opm('MEG',rawData));
cfg.viewmode    = 'butterfly';
cfg.colorgroups = 'allblack';
ft_databrowser(cfg,data_out_si_lp_hp_arft);

%% Trial def
% Here I am using a custom trial function which looks for triggers on a
% specific channel (e.g. NI-TRIG) of the unfiltered raw data loaded earlier
cfg                         = [];
cfg.rawData                 = rawData;
cfg.trialdef.trigchan       = 'NI-TRIG';
cfg.trialdef.downsample     = 1000;
cfg.correct_time            = 0.0;
cfg.trialdef.prestim        = 0.2;        % pre-stimulus interval
cfg.trialdef.poststim       = 0.5;        % post-stimulus interval
cfg.trialfun                = 'OPM_trialfun_usemat';
banana                      = ft_definetrial(cfg);

% Correct for optitrack
if run_num == 3
    banana.trl(:,1) = banana.trl(:,1)-rawData_MEG.time{1}(1)*1000;
    banana.trl(:,2) = banana.trl(:,2)-rawData_MEG.time{1}(1)*1000;
end

% Redefines the filtered data
cfg     = [];
data    = ft_redefinetrial(banana,data_out_si_lp_hp_arft);

%% Remove trials with any nan data (i.e. has been marked as artefactual)
trial2keep = [];
trial2reject = [];
count        = 1;
count2       = 1;

for t = 1:length(data.trial)
    result = sum(isnan(data.trial{t}(:)));
    if ~result
        trial2keep(count) = t;
        count=count+1;
    else
        trial2reject(count2) = t;
        count2       = count2+1;
    end
end

if run_num == 3
    save trial2keep trial2keep
    save trial2reject trial2reject
end

% Pick 438 trials (same as run 3)
s           = RandStream('mt19937ar','Seed',99);
trial2keep  = randsample(s,trial2keep,438);

% % remove bad trials
cfg                         = [];
cfg.trials                  = trial2keep;
data = ft_selectdata(cfg,data);

%% Save data
cd(save_dir);
save(['data_run' num2str(run_num) '.mat'],'data');

%% Perform timelockanalysis
% For ease, let's just do this for the TANS
cfg             = [];
cfg.channel    = ft_channelselection_opm('MEG',rawData);
avg_all         = ft_timelockanalysis(cfg,data);

cfg = [];
cfg.baseline = [-0.1 0];
[avg_all] = ft_timelockbaseline(cfg, avg_all);

cfg = [];
%cfg.ylim = [-455 455];
cfg.parameter = 'avg';
cfg.linewidth = 2;
cfg.colorgroups   = 'allblack';
ft_databrowser(cfg,avg_all);

cd(save_dir);
figure; plot(avg_all.time,avg_all.avg,'k','LineWidth',2);
ylabel('fT','FontSize',20);
xlabel('Time (s)','FontSize',20);
xlim([-0.1 0.4]);
if run_num == 3
    ylim([-200 200]);
else
ylim([-200 200]);
end
set(gca,'FontSize',16);
print(['run' num2str(run_num)],'-dpng','-r300');

%% Create and Plot 2D Layout (Fieldtrip)
cfg             = [];
cfg.output      = 'lay_123.mat';
cfg.grad        = data.grad;
cfg.channel     = data.label;
%cfg.headshape   = mesh;
cfg.rotate      = 0;
cfg.center      = 'yes';
cfg.projection  = 'polar';
cfg.channel     =  'all';
cfg.overlap     = 'keep';
lay_123         = ft_prepare_layout(cfg);

%figure;ft_plot_layout(lay_123)

%%
cfg             = [];
cfg.channel    = ft_channelselection_opm('TAN',rawData);
avg_all         = ft_timelockanalysis(cfg,data);

%% Plot Using ft_multiplot_ER
% Change the colormap to RdBu
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
cmap = colormap(flipud(brewermap(64,'RdBu'))); % change the colormap

cfg             = [];
cfg.parameter   = 'avg';
cfg.layout      = lay_123;
cfg.baseline    = [-0.1 0];
cfg.xlim        = [0.08 0.12];
cfg.linewidth   = 2;
%cfg.zlim        = [-200 200];
cfg.showlabels   = 'yes';
cfg.colormap    = cmap;
cfg.comment     = 'no';
figure; set(gcf,'position',[1 1 800 1000]);
ft_topoplotER(cfg,avg_all); hold on;
c = colorbar;
c.Location = 'southoutside';
c.FontSize = 20;
c.Label.String = 'fT';
print(['M100topo_run' num2str(run_num)],'-dpng','-r300');

cfg.xlim        = [0.15 0.25];
cfg.zlim        = [-500 500];
figure; set(gcf,'position',[1 1 800 1000]);
ft_topoplotER(cfg,avg_all); hold on;
c = colorbar;
c.Location = 'southoutside';
c.FontSize = 20;
c.Label.String = 'fT';
print(['M200topo_run' num2str(run_num)],'-dpng','-r300');

end






