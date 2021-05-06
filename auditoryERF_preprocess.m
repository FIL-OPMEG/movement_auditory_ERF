function auditoryERF_preprocess(data_dir,save_dir, run_num,...
    subject_name,motive_data)

%% Start preprocessing.
% Read in the raw data using BIDS
disp('Loading data...');
cfg             = [];
cfg.folder      = data_dir;
cfg.precision   = 'single';
cfg.bids.task   = 'aef';
cfg.bids.sub    = subject_name;
cfg.bids.ses    = '001';
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

% % % Plot using ft_databrowser
% cfg             = [];
% cfg.blocksize   = 30;
% %cfg.event       = banana;
% cfg.channel     = vertcat(ft_channelselection_opm('MEG',rawData));
% cfg.viewmode    = 'butterfly';
% cfg.colorgroups   = 'allblack';
% ft_databrowser(cfg,rawData);

%% Plot PSD
cfg                 = [];
cfg.channel         = vertcat(ft_channelselection_opm('MEG',rawData));
cfg.trial_length    = 10;
cfg.method          = 'tim';
cfg.foi             = [0 120];
cfg.plot            = 'yes';
[pow freq]          = ft_opm_psd(cfg,rawData);
ylim([1 1e4])
title('Raw Data');

%% Load in optitrack data
if run_num == 3
    % Load in optitrack data
    opti_data = csv2mat_sm(motive_data);
    
    % Convert mm to cm
    opti_data = optitrack_to_cm(opti_data);
    
    % Plot the rotations
    plot_motive_rotation(opti_data,'euler')
    
    % Plot the translations
    plot_motive_translation(opti_data,'euler')
    
    % Plot mean marker error (Column 7)
    figure;
    plot(opti_data.time,opti_data.rigidbodies.data(:,7),'LineWidth',2);
    ylabel('Mean Marker Error');xlabel('Time (s)');
    drawnow;
    
    %% Sync up opti-track and rawData
    [MovementDataOut, OPMdataOut] = syncOptitrackAndOPMdata(opti_data,...
        rawData,'TriggerChannelName','FluxZ-A');
    
    cd(save_dir);
    save(['MovementDataOut_run' num2str(run_num)], 'MovementDataOut')
    
end

%% Select the OPM Data
cfg             = [];
cfg.channel     = vertcat(ft_channelselection_opm('MEG',rawData));
if run_num == 3
    rawData_MEG     = ft_selectdata(cfg,OPMdataOut);
else
    rawData_MEG     = ft_selectdata(cfg,rawData);
end

%% Regress the optitrack data from the
if run_num == 3
    ref                 = (MovementDataOut.rigidbodies.data(:,1:6));
    % LP-filter optitrack data
    [ref]          = ft_preproc_lowpassfilter(ref', 1000, 2, 5);
    ref            = ref';
        
    % Regress
    [rawData_MEG_reg]   = regress_motive_OPMdata(rawData_MEG,ref,10);
else
    rawData_MEG_reg     = rawData_MEG;
end

%% MFC
% Please contact t.tierney@ucl.ac.uk for this script
[data_out_mfc, M, chan_inds] = ft_denoise_mfc(rawData_MEG_reg);

%% Plot data
cfg             = [];
cfg.blocksize   = 30;
%cfg.channel     = vertcat(ft_channelselection_opm('MEG',rawData));
cfg.viewmode    = 'butterfly';
cfg.colorgroups = 'allblack';
ft_databrowser(cfg,rawData_MEG);

%% Plot PSD
cfg                 = [];
cfg.channel         = vertcat(ft_channelselection_opm('MEG',rawData));
cfg.trial_length    = 10;
cfg.method          = 'tim';
cfg.foi             = [0 100];
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
cfg                     = [];
cfg.hpfilter            = 'yes';
cfg.hpfreq              = 2;
cfg.filtord             = 5;
cfg.hpinstabilityfix    = 'reduce';
%cfg.hpfilttype = 'fir';
data_out_si_lp_hp       = ft_preprocessing(cfg,data_out_si_lp);
data_out_si_hp          = ft_preprocessing(cfg,data_out_si);

%% Remove DS (and 17 for 002) - channels are bad
cfg                 = [];
if strcmp(subject_name,'001')
    cfg.channel         = vertcat(data_out_si_lp.label,'-DS-TAN','-DS-RAD');
elseif strcmp(subject_name,'002')
    cfg.channel         = vertcat(data_out_si_lp.label,'-17-TAN','-17-RAD','-DS-RAD','-DS-TAN');
end
data_out_si_lp_hp   = ft_selectdata(cfg,data_out_si_lp_hp);
data_out_si_hp      = ft_selectdata(cfg,data_out_si_hp);

if run_num == 3
    save data_out_si_lp_hp data_out_si_lp_hp
end

%% Plot data for Artifact Rejection
cfg             = [];
cfg.blocksize   = 30;
%cfg.channel     = vertcat(ft_channelselection_opm('MEG',rawData));
cfg.viewmode    = 'vertical';
cfg.colorgroups = 'allblack';
arft            = ft_databrowser(cfg,data_out_si_hp);

%% Turn the highlighted data into nans
arft.artfctdef.reject          = 'nan';
data_out_si_lp_hp_arft = ft_rejectartifact(arft, data_out_si_lp_hp);

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
    banana.trl(:,1) = banana.trl(:,1)-round(OPMdataOut.time{1}(1)*1000);
    banana.trl(:,2) = banana.trl(:,2)-round(OPMdataOut.time{1}(1)*1000);
    trl_index       = banana.trl(:,1);
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
    trl_index(trial2reject) = [];
    save trial2keep trial2keep
    save trial2reject trial2reject
    save trl_index trl_index
end

% Pick 550 trials (same as run 3)
s           = RandStream('mt19937ar','Seed',99);
trial2keep  = randsample(s,trial2keep,550);

% % remove bad trials
cfg                         = [];
cfg.trials                  = trial2keep;
data = ft_selectdata(cfg,data);

%% Save data
cd(save_dir);
save(['data_run' num2str(run_num) '.mat'],'data');

%% Perform timelockanalysis
cfg             = [];
cfg.channel    = 'all';
avg_all         = ft_timelockanalysis(cfg,data);

cfg = [];
cfg.baseline = [-0.1 0];
[avg_all] = ft_timelockbaseline(cfg, avg_all);

% Plot in fT
cfg = [];
%cfg.ylim = [-455 455];
cfg.parameter = 'avg';
cfg.linewidth = 2;
cfg.colorgroups   = 'allblack';
ft_databrowser(cfg,avg_all);

% Convert to t-value
epoched_dataset = [];

for i = 1:length(data.trial)
    epoched_dataset(:,:,i) = data.trial{1,i};
end

SE = std(epoched_dataset,[],3)/sqrt(size(epoched_dataset,3));
avg_all.t_value = avg_all.avg./SE;

% Plot t-value
cd(save_dir);
figure; plot(avg_all.time,avg_all.t_value,'k','LineWidth',2);
ylabel('t-value','FontSize',20);
xlabel('Time (s)','FontSize',20);
xlim([-0.1 0.4]);
ylim([-12 12]);
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

%% Select TAN channels
cfg             = [];
cfg.channel    = ft_channelselection_opm('TAN',rawData);
avg_all         = ft_timelockanalysis(cfg,data);

% Select TAN channels from data
cfg             = [];
cfg.channel    = ft_channelselection_opm('TAN',rawData);
data_TAN        = ft_selectdata(cfg,data);

epoched_dataset = [];

for i = 1:length(data_TAN.trial)
    epoched_dataset(:,:,i) = data_TAN.trial{1,i};
end

SE = std(epoched_dataset,[],3)/sqrt(size(epoched_dataset,3));
avg_all.t_value = avg_all.avg./SE;

%% Plot Using ft_multiplot_ER
% Change the colormap to RdBu
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
cmap = colormap(flipud(brewermap(64,'RdBu'))); % change the colormap

cfg             = [];
cfg.parameter   = 't_value';
cfg.layout      = lay_123;
%cfg.baseline    = [-0.1 0];
cfg.xlim        = [0.08 0.12];
cfg.zlim        = [-8 8];
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
c.Label.String = 't-value';
print(['M100topo_run' num2str(run_num)],'-dpng','-r300');

%% Plot Sensor Level PSD for low freqs
cfg                 = [];
cfg.trial_length    = 10;
cfg.method          = 'tim';
cfg.foi             = [2 5];
cfg.plot            = 'yes';
cfg.plot_legend      = 'no';
[pow freq]          = ft_opm_psd(cfg,data_out_si_lp_hp);
print(['PSD_sensor_level' num2str(run_num)],'-dpng','-r300');

end






