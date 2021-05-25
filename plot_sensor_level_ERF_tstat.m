function plot_sensor_level_ERF_tstat(save_dir,sens_name)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to load the raw and processed sensor-level data from all 
% three runs, calculate the auditory ERF from the OPM channel with 
% max t-value and plot.
%
% This corresponds to Supplementary Figure 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(save_dir);

avg_all_d = [];

for i = 1:3
    disp('Loading data...');
    d = load(['data_unprocessed' num2str(i) '.mat'])
    data = d.data
    clear d;
    %% Perform timelockanalysis
    cfg             = [];
    cfg.channel    = 'all';
    avg_all         = ft_timelockanalysis(cfg,data);
    
    cfg = [];
    cfg.baseline = [-0.1 0];
    [avg_all] = ft_timelockbaseline(cfg, avg_all);
    
    % Convert to t-value
    epoched_dataset = [];
    
    for j = 1:length(data.trial)
        epoched_dataset(:,:,j) = data.trial{1,j};
    end
    
    SE = std(epoched_dataset,[],3)/sqrt(size(epoched_dataset,3));
    avg_all.t_value = avg_all.avg./SE;
    
    avg_all_d{i} = avg_all;
    clear avg_all
end


figure;
cfg = [];
cfg.channel = {sens_name};
cfg.parameter = 't_value';
cfg.baseline = [-0.1 0];
cfg.showlegend    = 'no';
cfg.xlim    = [-0.1 0.4];
cfg.linecolor = [0.4275    0.9804    0.3922;0.9804    0.5686    0.3843;
    0.2824    0.3137    0.9804];
cfg.linewidth = 2;
cfg.ylim = [-10 10];
ft_singleplotER(cfg,avg_all_d{1},avg_all_d{2},avg_all_d{3});
set(gca,'FontSize',18);
xlabel('Time (s)','FontSize',20);
ylabel('t-value','FontSize',20);
title('');
print('run123_unprocessed_ERF','-dpng','-r300');
legend({'Sitting';'Standing';'Standing + Moving'},'Location','SouthOutside');

avg_all_d = [];

for i = 1:3
    disp('Loading data...');
    d = load(['data_run' num2str(i) '.mat'])
    data = d.data
    clear d;
    %% Perform timelockanalysis
    cfg             = [];
    cfg.channel    = 'all';
    avg_all         = ft_timelockanalysis(cfg,data);
    
    cfg = [];
    cfg.baseline = [-0.1 0];
    [avg_all] = ft_timelockbaseline(cfg, avg_all);
    
    % Convert to t-value
    epoched_dataset = [];
    
    for j = 1:length(data.trial)
        epoched_dataset(:,:,j) = data.trial{1,j};
    end
    
    SE = std(epoched_dataset,[],3)/sqrt(size(epoched_dataset,3));
    avg_all.t_value = avg_all.avg./SE;
    
    avg_all_d{i} = avg_all;
    clear avg_all
end


figure;
cfg = [];
cfg.channel = {sens_name};
cfg.parameter = 't_value';
cfg.baseline = [-0.1 0];
cfg.showlegend    = 'no';
cfg.xlim    = [-0.1 0.4];
cfg.linecolor = [0.4275    0.9804    0.3922;0.9804    0.5686    0.3843;
    0.2824    0.3137    0.9804];
cfg.linewidth = 2;
cfg.ylim = [-10 10];
ft_singleplotER(cfg,avg_all_d{1},avg_all_d{2},avg_all_d{3});
set(gca,'FontSize',18);
xlabel('Time (s)','FontSize',20);
ylabel('t-value','FontSize',20);
title('');
print('run123_sensor_level_ERF','-dpng','-r300');
legend({'Sitting';'Standing';'Standing + Moving'},'Location','SouthOutside');
end



