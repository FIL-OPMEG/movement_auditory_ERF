%% Paths (RS)
fieldtripDir    = 'D:\scripts\fieldtrip-master';
script_dir      = 'D:\Github\analyse_OPMEG';
data_dir        = 'D:\data\20201208_optitrack';
save_dir        = 'D:\data\20201208_optitrack';
mocap_func      = 'D:\scripts\motioncapture_functions';

% Add Fieldtrip to path
disp('Adding Fieldtrip and analyse_OPMEG to your MATLAB path');
addpath(fieldtripDir)
ft_defaults;

% Add analyse_OPMEG Scripts to path
addpath(genpath(script_dir));
addpath(mocap_func);

% cd to save dir
cd(save_dir)

run_num = 3; %2 or %3

%% (2) Start preprocessing.
% Read in the raw data using BIDS
cfg             = [];
cfg.folder      = data_dir;
cfg.precision   = 'single';
cfg.bids.task   = 'MMF';
cfg.bids.sub    = 'NA';
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

% % Plot using ft_databrowser
cfg             = [];
cfg.blocksize   = 30;
cfg.event       = banana;
%cfg.channel     = vertcat(ft_channelselection_opm('MEG',rawData));
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

%% Load in optitrack data
if run_num == 3
    opti_data = csv2mat_sm('Take 2020-12-08 03.42.02 PM_002_forVisualisation.csv');
    
    figure;
    set(gcf,'Position',[1 1 1000 900]);
    subplot(4,1,1);
    plot(opti_data.time,opti_data.rigidbodies.data(:,1),'k','LineWidth',2);
    set(gca,'FontSize',12);
    title('Rigid Body Rotation: X','FontSize',16);
    subplot(4,1,2);
    plot(opti_data.time,opti_data.rigidbodies.data(:,2),'k','LineWidth',2);
    set(gca,'FontSize',12);
    title('Rigid Body Rotation: Y','FontSize',16);
    subplot(4,1,3);
    plot(opti_data.time,opti_data.rigidbodies.data(:,3),'k','LineWidth',2);
    set(gca,'FontSize',12);
    title('Rigid Body Rotation: Z','FontSize',16);
    subplot(4,1,4);
    plot(opti_data.time,opti_data.rigidbodies.data(:,4),'k','LineWidth',2);
    title('Rigid Body Rotation: W','FontSize',16);
    set(gca,'FontSize',12);
    xlabel('Time (s)','FontSize',16);
    print('opti_rot','-dpng','-r300');
    
    figure;
    set(gcf,'Position',[1 1 1000 900]);
    subplot(3,1,1);
    plot(opti_data.time,opti_data.rigidbodies.data(:,5)/10,'k','LineWidth',2);
    set(gca,'FontSize',12);
    ylabel('cm','FontSize',16);
    title('Rigid Body Position: X','FontSize',16);
    subplot(3,1,2);
    plot(opti_data.time,opti_data.rigidbodies.data(:,6)/10,'k','LineWidth',2);
    set(gca,'FontSize',12);
    ylabel('cm','FontSize',16);
    title('Rigid Body Position: Y','FontSize',16);
    subplot(3,1,3);
    plot(opti_data.time,opti_data.rigidbodies.data(:,7)/10,'k','LineWidth',2);
    set(gca,'FontSize',12);
    ylabel('cm','FontSize',16);
    title('Rigid Body Position: Z','FontSize',16);
    xlabel('Time (s)','FontSize',16);
    print('opti_pos','-dpng','-r300');
    
%     figure;
%     plot(opti_data.time,opti_data.rigidbodies.data(:,8),'LineWidth',2);
%     ylabel('Mean Marker Error');
%     xlabel('Time (s)');
    
    %% Sync up opti-track and rawData
    [MovementDataOut, OPMdataOut] = syncOptitrackAndOPMdata(opti_data,...
        rawData,'TriggerChannelName','FluxZ-A');
    save MovementDataOut MovementDataOut
end

%% Plot PSD
cfg                 = [];
cfg.channel         = vertcat(ft_channelselection_opm('MEG',rawData));
cfg.trial_length    = 10;
cfg.method          = 'tim';
cfg.foi             = [0 10];
cfg.plot            = 'yes';
cfg.plot_legend     = 'no';
[pow freq]          = ft_opm_psd_compare(cfg,rawData_MEG);
ylim([1 1e4])

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
    
%     % Convert quaternions to euler angles
%     eulZYX = quat2eul([MovementDataOut.rigidbodies.data(:,4), MovementDataOut.rigidbodies.data(:,1:3)]);
    
    megind          = rawData_MEG.trial{1};
    megres          = zeros(size(megind));
    %ref             = MovementDataOut.rigidbodies.data(:,5:7);
    ref             = (MovementDataOut.rigidbodies.data(:,1:6));
        
    
    ref_size        = size(ref,1);
    winSize         = size(megind,2);
    
    % add a mean column to the reference regressors
    intercept       = ones(winSize,1);
    
    reference = ref;
    reference       = [reference ones(size(reference,1),1)];
    
    % Start of Window-Loop
    % Calculate window size in terms of number of data points
    wsize = 10*rawData_MEG.fsample;
    fprintf('Using %d data-points per window\n',wsize);
    
    % Create array of zeros for data
    meg_ind_synth_grad  = zeros(size(megind));
    % Create array of zeros for triangular weighting
    a                   = zeros(1,size(megind,2));
    
    % Weighting? I could probably take this out later
    w=ones(size(megind));
    if size(w,1)==1; w=repmat(w,1,size(megind,1)); end
    
    % Start at 0
    offset=0;
    ft_progress('init', 'etf', 'Regressing...')
    
    while true
        ft_progress(offset/size(megind,2))
        
        % Calculate start and stop points for this loop
        start=offset+1;
        stop=min(size(megind,2),offset+wsize);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % These values are grown by a factor of 5
        % Is this needed??
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        counter=0;
        while any (sum(min(w(start:stop),1))) <wsize
            if counter <= 0 ; break; end
            start=max(1,start-wsize/2);
            stop=min(size(megind,2),stop+wsize/2);
            counter=counter-1;
        end
        if rem(stop-start+1,2)==1; stop=stop-1; end
        wsize2=stop-start+1;
        
        % Do gradiometry
        beta    = pinv(reference(start:stop,:))*megind(:,start:stop)';
        yy      = (megind(:,start:stop)' - reference(start:stop,:)*beta)';
        
        % triangular weighting (specified via b variable)
        if start==1
            b=[ones(1,wsize2/2)*wsize2/2, wsize2/2:-1:1];
        elseif stop==size(megind,2)
            b=[1:wsize2/2, ones(1,wsize2/2)*wsize2/2];
        else
            b=[1:wsize2/2, wsize2/2:-1:1];
        end
        
        % Add to meg_ind_synth_grad variable outside the loop, weighted by b
        meg_ind_synth_grad(:,start:stop)=meg_ind_synth_grad(:,start:stop)+bsxfun(@times,yy,b);
        
        % Add triangular weighting to variable outside the loop
        a(1,start:stop)=a(start:stop)+b;
        
        % Adjust offset parameter by window size divided by 5
        offset=offset+wsize/5;
        
        % If we have reached the end of the data BREAK
        if offset>size(megind,2)-wsize/5; break; end
    end
    ft_progress('close')
    
    % Adjust for triangular weighting
    meg_ind_synth_grad=bsxfun(@times,meg_ind_synth_grad,1./a);
    % Find any NaN values and convert to 0
    meg_ind_synth_grad(isnan(meg_ind_synth_grad))=0;
    
    data_out_reg            = rawData_MEG;
    data_out_reg.trial{1}   = meg_ind_synth_grad;
    
else
    data_out_reg = rawData_MEG;
end

%% Plot data
cfg             = [];
cfg.blocksize   = 30;
%cfg.channel     = vertcat(ft_channelselection_opm('MEG',rawData));
cfg.viewmode    = 'butterfly';
cfg.colorgroups = 'allblack';
ft_databrowser(cfg,data_out_reg);

%% MFC
[data_out_mfc, M, chan_inds] = ft_denoise_mfc(data_out_reg);

%% Plot data
cfg             = [];
cfg.blocksize   = 400;
%cfg.channel     = vertcat(ft_channelselection_opm('MEG',rawData));
cfg.viewmode    = 'butterfly';
cfg.colorgroups = 'allblack';
ft_databrowser(cfg,data_out_reg);

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
%data_out_mfc_si_lp  = ft_preprocessing(cfg,data_out_mfc_si);
data_out_si_lp      = ft_preprocessing(cfg,data_out_si);

%% HP-filter
cfg                 = [];
cfg.hpfilter        = 'yes';
cfg.hpfreq          = 2;
cfg.filtord         = 2;
cfg.hpinstabilityfix = 'reduce';
%cfg.hpfilttype = 'fir';
data_out_si_lp_hp     = ft_preprocessing(cfg,data_out_si_lp);

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

%% Trial def
cd('D:\data\20201208_optitrack\sub-NA\ses-001\meg');

cfg = [];
if run_num == 1
    cfg.dataset                 = 'sub-NA_ses-001_task-MMF_run-001_meg.bin';
elseif run_num == 2
    cfg.dataset                 = 'sub-NA_ses-001_task-MMF_run-002_meg.bin';
elseif run_num == 3
    cfg.dataset                 = 'sub-NA_ses-001_task-MMF_run-003_meg.bin';
end
cfg.trialdef.trigchan       = 'NI-TRIG';
cfg.trialdef.downsample     = 1000;
cfg.correct_time            = 0.0;
cfg.trialdef.prestim        = 0.2;        % pre-stimulus interval
cfg.trialdef.poststim       = 0.5;        % post-stimulus interval
cfg.trialfun                = 'OPM_TrialFun_RS';
banana                      = ft_definetrial(cfg);

% Correct for optitrack
if run_num == 3
    banana.trl(:,1) = banana.trl(:,1)-rawData_MEG.time{1}(1)*1000;
    banana.trl(:,2) = banana.trl(:,2)-rawData_MEG.time{1}(1)*1000;
end

% trl_index = int64(banana.trl(:,1));
% MovementDataOut.time(trl_index(1));
% 
% save trl_index trl_index

% Redefines the filtered data
cfg = [];
data = ft_redefinetrial(banana,data_out_si_lp_hp_arft);

%% Remove trials with any nan data
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

%%
reject_cm = [];

for t = 1:length(trial2reject)
    pos_data = MovementDataOut.rigidbodies.data(banana.trl(trial2reject(t),1)...
        :banana.trl(trial2reject(t),2), 5:7)/10;
    
    pos_data_mean = mean(pos_data,1);
    
    distances = pdist2(pos_data, pos_data_mean);
    % Find the max distance
    reject_cm(t) = max(distances(:));
end

keep_cm = [];

for t = 1:length(trial2keep)
    pos_data = MovementDataOut.rigidbodies.data(banana.trl(trial2keep(t),1)...
        :banana.trl(trial2keep(t),2), 5:7)/10;
    
    pos_data_mean = mean(pos_data,1);
    
    distances = pdist2(pos_data, pos_data_mean);
    % Find the max distance
    keep_cm(t) = max(distances(:));
        
%     [index1, index2] = find(distances == maxDistance);
%     pdist2(pos_data(index2(1),:),pos_data(index2(2),:));
end

grp = [zeros(1,length(reject_cm)),ones(1,length(keep_cm))];
figure;
boxplot([reject_cm,keep_cm],grp);

t = table([reject_cm keep_cm]',grp');
writetable(t,'trial_movt.csv','Delimiter',',')

%% Now let's do the same for Nic's speed measure
load('speedRotData_byTrial.mat');
keep_cm         = mean_trlV(trial2keep)'/10;
reject_cm       = mean_trlV(trial2reject)'/10;

grp = [zeros(1,length(reject_cm)),ones(1,length(keep_cm))];
figure;
boxplot([reject_cm,keep_cm],grp);

t = table([reject_cm keep_cm]',grp');
writetable(t,'speed_opti.csv','Delimiter',',')

%% Now let's do the same for Nic's rotation speed measure
keep_cm         = mean_trlQv(trial2keep)';
reject_cm       = mean_trlQv(trial2reject)';

grp = [zeros(1,length(reject_cm)),ones(1,length(keep_cm))];
figure;
boxplot([reject_cm,keep_cm],grp);

t = table([reject_cm keep_cm]',grp');
writetable(t,'rotation_speed_opti.csv','Delimiter',',')

%%
% Do median split based on head movement
indx_low    = find(keep_cm<median(keep_cm));
indx_high   = find(keep_cm>median(keep_cm));

% % remove bad trials
cfg                         = [];
cfg.trials                  = indx_high;
data_high = ft_selectdata(cfg,data);
cfg.trials                  = indx_low;
data_low = ft_selectdata(cfg,data);


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
    ylim([-700 700]);
else
ylim([-300 300]);
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

figure;ft_plot_layout(lay_123)

%%
cfg             = [];
cfg.channel    = ft_channelselection_opm('TAN',rawData);
avg_all         = ft_timelockanalysis(cfg,data);

cfg = [];
cfg.baseline = [-0.1 0];
[avg_all] = ft_timelockbaseline(cfg, avg_all);

%% Plot Using ft_multiplot_ER
% Change the colormap to RdBu
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
cmap = colormap(flipud(brewermap(64,'RdBu'))); % change the colormap

cfg             = [];
cfg.parameter   = 'avg';
cfg.layout      = lay_123;
cfg.baseline    = [-0.2 0];
cfg.linewidth   = 2;
cfg.xlim        = [0.08 0.12];
cfg.zlim        = [-500 500];
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








