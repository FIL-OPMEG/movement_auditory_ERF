function plot_psd_preprocessing(results_dir1,results_dir2)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to plot the PSD from 2-10Hz for all three runs, using:
% 1. Raw Data
% 2. Pre-processed sensor-level data
% 3. Beamformed data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Raw Data
% Subject 001
cd(results_dir1);
disp('Loading data...');
d{1} = load('rawData_MEG1.mat');
d{2} = load('rawData_MEG2.mat');
d{3} = load('rawData_MEG3.mat');

cfg                 = [];
cfg.channel         = {'N3-TAN'};
cfg.trial_length    = 10;
cfg.method          = 'tim';
cfg.foi             = [0 40];
cfg.plot            = 'no';
cfg.plot_legend      = 'no';
[p{1} freq]      = ft_opm_psd(cfg,d{1}.rawData_MEG);
[p{2} freq]      = ft_opm_psd(cfg,d{2}.rawData_MEG);
[p{3} freq]      = ft_opm_psd(cfg,d{3}.rawData_MEG);

% Subject 002
cd(results_dir2);
disp('Loading data...');
d{1} = load('rawData_MEG1.mat');
d{2} = load('rawData_MEG2.mat');
d{3} = load('rawData_MEG3.mat');

cfg                 = [];
cfg.channel         = {'1A-TAN'};
cfg.trial_length    = 10;
cfg.method          = 'tim';
cfg.foi             = [0 40];
cfg.plot            = 'no';
cfg.plot_legend     = 'no';
[p{4} freq]         = ft_opm_psd(cfg,d{1}.rawData_MEG);
[p{5} freq]         = ft_opm_psd(cfg,d{2}.rawData_MEG);
[p{6} freq]         = ft_opm_psd(cfg,d{3}.rawData_MEG);

% Plot
figure;
cols = [0.4275    0.9804    0.3922;0.9804    0.5686    0.3843;
                0.2824    0.3137    0.9804];

for i=1:3
po      = (mean(p{i}(:,:,:),3)+ mean(p{i+3}(:,:,:),3))/2;
po_mean = squeeze(mean(po,2));

set(gcf,'Position',[100 100 800 800]);
fig= gcf;
fig.Color=[1,1,1]; hold on;

% Plot all channels
h = plot(freq,po_mean','Color',cols(i,:),'LineWidth',3);
h.Color(4) = 0.9;
end
grid on
ax = gca; % current axes
ax.FontSize = 30;
ax.TickLength = [0.02 0.02];
set(gca, 'YScale', 'log'); hold on;

xlabel('Frequency (Hz)','FontSize',37)
labY = ['$$PSD (' 'fT' ' \sqrt[-1]{Hz}$$)'];
ylabel(labY,'interpreter','latex','FontSize',37)
xlim([0 100]);
ylim([10 1e+4]);
print('PSD_rawData_0_100' ,'-dpng','-r300');
xlim([2 10]);
print('PSD_rawData_0_10' ,'-dpng','-r300');

%% HP-Filter
% Subject 001
cd(results_dir1);
disp('Loading data...');
d{1} = load('data_out_si_hp1.mat');
d{2} = load('data_out_si_hp2.mat');
d{3} = load('data_out_si_hp3.mat');

cfg                 = [];
cfg.channel         = {'N3-TAN'};
cfg.trial_length    = 10;
cfg.method          = 'tim';
cfg.foi             = [0 40];
cfg.plot            = 'no';
cfg.plot_legend      = 'no';
[p{1} freq]      = ft_opm_psd(cfg,d{1}.data_out_si_hp);
[p{2} freq]      = ft_opm_psd(cfg,d{2}.data_out_si_hp);
[p{3} freq]      = ft_opm_psd(cfg,d{3}.data_out_si_hp);

% Subject 002
cd(results_dir2);
disp('Loading data...');
d{1} = load('data_out_si_hp1.mat');
d{2} = load('data_out_si_hp2.mat');
d{3} = load('data_out_si_hp3.mat');

cfg                 = [];
cfg.channel         = {'1A-TAN'};
cfg.trial_length    = 10;
cfg.method          = 'tim';
cfg.foi             = [0 40];
cfg.plot            = 'no';
cfg.plot_legend     = 'no';
[p{4} freq]         = ft_opm_psd(cfg,d{1}.data_out_si_hp);
[p{5} freq]         = ft_opm_psd(cfg,d{2}.data_out_si_hp);
[p{6} freq]         = ft_opm_psd(cfg,d{3}.data_out_si_hp);

% Plot
figure;
cols = [0.4275    0.9804    0.3922;0.9804    0.5686    0.3843;
                0.2824    0.3137    0.9804];

for i=1:3
po      = (mean(p{i}(:,:,:),3)+ mean(p{i+3}(:,:,:),3))/2;
po_mean = squeeze(mean(po,2));

set(gcf,'Position',[100 100 800 800]);
fig= gcf;
fig.Color=[1,1,1]; hold on;

% Plot all channels
h = plot(freq,po_mean','Color',cols(i,:),'LineWidth',3);
h.Color(4) = 0.9;
end
grid on
ax = gca; % current axes
ax.FontSize = 30;
ax.TickLength = [0.02 0.02];
set(gca, 'YScale', 'log'); hold on;

xlabel('Frequency (Hz)','FontSize',37)
labY = ['$$PSD (' 'fT' ' \sqrt[-1]{Hz}$$)'];
ylabel(labY,'interpreter','latex','FontSize',37)
xlim([0 100]);
ylim([10 1e+4]);
print('PSD_sensor_level_0_100' ,'-dpng','-r300');
xlim([2 10]);
print('PSD_sensor_level_0_10' ,'-dpng','-r300');

%% After Beamforming
disp('Loading data...');
cd(results_dir1);
d{1} = load('VE_1.mat');
d{2} = load('VE_2.mat');
d{3} = load('VE_3.mat');

cfg                 = [];
cfg.channel         = 'all';
cfg.trial_length    = 10;
cfg.method          = 'tim';
cfg.foi             = [0 40];
cfg.plot            = 'no';
cfg.plot_legend      = 'no';
[p{1} freq]      = ft_opm_psd(cfg,d{1}.VE);
[p{2} freq]      = ft_opm_psd(cfg,d{2}.VE);
[p{3} freq]      = ft_opm_psd(cfg,d{3}.VE);

% Subject 002
disp('Loading data...');
cd(results_dir2);
d{1} = load('VE_1.mat');
d{2} = load('VE_2.mat');
d{3} = load('VE_3.mat');

cfg                 = [];
cfg.channel         = 'all';
cfg.trial_length    = 10;
cfg.method          = 'tim';
cfg.foi             = [0 40];
cfg.plot            = 'no';
cfg.plot_legend      = 'no';
[p{4} freq]      = ft_opm_psd(cfg,d{1}.VE);
[p{5} freq]      = ft_opm_psd(cfg,d{2}.VE);
[p{6} freq]      = ft_opm_psd(cfg,d{3}.VE);


% Plot
figure;
cols = [0.4275    0.9804    0.3922;0.9804    0.5686    0.3843;
                0.2824    0.3137    0.9804];
            
for i=1:3
po      = (mean(p{i}(:,:,:),3) + mean(p{i+3}(:,:,:),3))/2;
po_mean = squeeze(mean(po,2));

set(gcf,'Position',[100 100 800 800]);
fig= gcf;
fig.Color=[1,1,1]; hold on;

% Plot all channels
h = plot(freq,po_mean','Color',cols(i,:),'LineWidth',3);
h.Color(4) = 0.9;
end
grid on
ax = gca; % current axes
ax.FontSize = 30;
ax.TickLength = [0.02 0.02];
set(gca, 'YScale', 'log'); hold on;

xlabel('Frequency (Hz)','FontSize',37)
labY = ['$$PSD (' 'dB' '/{Hz}$$)'];
ylabel(labY,'interpreter','latex','FontSize',37)
ylim([3e-9 2.6e-8]);
%print('PSD_beamformer_0_40' ,'-dpng','-r300');
xlim([2 10]);
print('PSD_beamformer_0_10' ,'-dpng','-r300');
end
