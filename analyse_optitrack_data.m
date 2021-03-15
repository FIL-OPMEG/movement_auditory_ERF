function analyse_optitrack_data(MovementDataOut,trl_index,log_array,trial2keep)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Script for analysing the optitrack data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
Fs  = 1000; % Sampling Rate

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Position Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot the position data
pos_data = MovementDataOut.rigidbodies.data(:,4:6)/10; % Convert to cm
pos_data = pos_data(trl_index(1):trl_index(end)+500,:); % Trim to start of the tones

%% Here we are converting to the MRI coordinate system where:
% - X = Left-Right
% - Y = Forward-Back
% - Z = Up-Down
pos_data = [pos_data(:,3) pos_data(:,1)*-1 pos_data(:,2)];

%%
% Make the first point 0 0 0
first_point = repmat(pos_data(1,:),size(pos_data,1),1);
pos_data = pos_data-first_point; clear first_point

t = [0:1:length(pos_data)-1]/Fs;
%log_array = log_array(trl_index(1):trl_index(end)+500);

%% Generate colors
cols = [100, 100, 250; 142 185 57; 198 61 61]./255;

%% Plot
coord = {'X','Y','Z'};
coord_name = {'X: Right-Left','Y:Forward-Back','Z: Up-Down'};

figure;
set(gcf,'Position',[1 1 1200 500]);

for c = 1:size(pos_data,2)
    m = pos_data(:,c);
    %m(log_array,:) = NaN;
    p1 = plot(t,m,'LineWidth',2,'Color',cols(c,:)); hold on;
    p1.Color(4) = 0.7; 
end
%plot(t,m,'r','LineWidth',2);
ylabel(['Distance (cm)'],'FontSize',24);
set(gca,'FontSize',18);
xlabel('Time (s)','FontSize',24);
%legend(coord_name,'Location','EastOutside');
print('opti_pos_reject','-dpng','-r300');

%% Get euclidian distance from start point
dist_from_start = zeros(length(t),1);

for time = 1:length(t)
    dist_from_start(time) = pdist2(pos_data(1,:),pos_data(time,:));
end

dist_from_start = pdist2(pos_data(1,:),pos_data(:,:));

figure;
set(gcf,'Position',[1 1 1200 600]);
m = dist_from_start(:,:);
p1 = plot(t,m,'LineWidth',2,'Color','k'); hold on;
p1.Color(4) = 0.8;
set(gca,'FontSize',18);
ylabel(['Euclidian Distane (cm)'],'FontSize',22);
xlabel('Time (s)','FontSize',22);
print('euclidian_distance_from_start','-dpng','-r300');

%% Plot some histograms

for c = 1:3
    figure; h = histfit(pos_data(:,c),50); xlim([-50 50]);
    h(1).FaceColor = cols(c,:);
    h(2).Color = [.2 .2 .2];
    set(gca,'FontSize',18);
    xlabel('Distance (cm)','FontSize',22);
    ylabel('Probability','FontSize',22);
    yt = get(gca, 'YTick');
    set(gca,'YTick',linspace(yt(1),yt(end),4));
    yt = get(gca, 'YTick');
    norm_vals = round((yt/length(pos_data(:,c))),2);
    xt = get(gca, 'XTick');
    set(gca,'XTick',linspace(xt(1),xt(end),5));
    set(gca, 'YTick', yt, 'YTickLabel',norm_vals);
    print(['histogram' num2str(c)],'-dpng','-r300');
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotation Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the rotation data and trim to the start of the experiment
degXYZ = (MovementDataOut.rigidbodies.data(:,1:3));
degXYZ = degXYZ(trl_index(1):trl_index(end)+500,:); % Trim to start of the tones

%% Here we are converting to the MRI coordinate system where:
% - X = Pitch
% - Y = Yaw
% - Z = Roll

degXYZ = [degXYZ(:,2)*-1 degXYZ(:,1)*-1 degXYZ(:,3)*-1];

%%
% Make the first point 0 0 0
first_point = repmat(degXYZ(1,:),size(degXYZ,1),1);
degXYZ = degXYZ-first_point; clear first_point

%% Plot
cols = [235 210 0; 230 0 99; 23 196 230]/255;

coord = {'X','Y','Z'};
coord_name = {'X: Pitch','Y: Yaw','Z: Roll'};

figure;
set(gcf,'Position',[1 1 1200 500]);

for c = 1:size(degXYZ,2)
    m = degXYZ(:,c);
    %m(log_array,:) = NaN;
    p1 = plot(t,m,'LineWidth',2,'Color',cols(c,:)); hold on;
    p1.Color(4) = 0.7; 
end
%plot(t,m,'r','LineWidth',2);
ylabel(['Degrees (°)'],'FontSize',24);
set(gca,'FontSize',18);
xlabel('Time (s)','FontSize',24);
%legend(coord_name,'Location','EastOutside');
print('opti_rot_reject','-dpng','-r300');

%% Polar Histogram
for c = 1:3
    figure; ax = polaraxes;
    h = polarhistogram(deg2rad(degXYZ(:,c)),50,'Normalization','probability');
    thetalim([-30 30]);
    ax.FontSizeMode = 'manual'
    ax.FontSize = 18
    ax.RAxis.FontSize = 14;
    ax.ThetaTickLabel = {'-30°','-15°','0°','15°','30°'}
    
    h(1).FaceColor = cols(c,:);
    %h(2).Color = [.2 .2 .2];
    %xlabel('Degrees (°)','FontSize',22);
    %     ylabel('Sample Count','FontSize',22);
    print(['polarhistogram' num2str(c)],'-dpng','-r300');
end

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate XYZ per trial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trl_index = trl_index-trl_index(1);
trl_index = trl_index(trial2keep);

keep_cm = [];

% Position
for k = 1:length(trial2keep)
    pos_data_t = pos_data(trl_index(k):trl_index(k)+500,:);
    for d = 1:3
        distances = pdist2(pos_data_t(:,d), mean(pos_data_t(:,d),2));
        % Find the max distance
        keep_cm(k,d) = max(distances(:));
    end
end

% Degrees
for k = 1:length(trial2keep)
    deg_data_t = degXYZ(trl_index(k):trl_index(k)+500,:);
    for d = 1:3
        distances = pdist2(deg_data_t(:,d), mean(deg_data_t(:,d),2));
        % Find the max distance
        keep_cm(k,d+3) = max(distances(:));
    end
end

% Euclidian Distance
for k = 1:length(trial2keep)
    euclidian_data_t = dist_from_start(trl_index(k):trl_index(k)+500);
    distances = pdist2(euclidian_data_t', euclidian_data_t');
    % Find the max distance
    keep_cm(k,7) = max(distances(:));
end

%
figure;
boxplot(keep_cm)

t = array2table(keep_cm,'VariableNames',{'PosX','PosY','PosZ','RotX',...
'RotY','RotZ','EuclidianDist'});
writetable(t,'trial_movt.csv','Delimiter',',')







%% Load in NA mesh
% Load in the MRI
scannercast_loc = 'D:\Github\scannercast\examples\NA';
cd(scannercast_loc);

mri = ft_read_mri('NA.nii');
mri = ft_determine_coordsys(mri);
mri.coordsys = 'neuromag';

%% Extract Scalp Surface from the MRI and create a mesh
cfg                     = [];
cfg.output              = 'scalp';
cfg.scalpsmooth         = 5;
cfg.scalpthreshold      = 0.08; % Change this value if the mesh looks weird
scalp                   = ft_volumesegment(cfg, mri);

% Create mesh
cfg                     = [];
cfg.method              = 'isosurface';
cfg.numvertices         = 5000;
mesh                    = ft_prepare_mesh(cfg,scalp);
mesh                    = ft_convert_units(mesh,'cm');

figure;
ft_plot_mesh(mesh,'facecolor',[238,206,179]./255,'EdgeColor','none',...
    'facealpha',0.8,'facecolor','skin'); camlight;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Let's make an animated plot!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Bit of code to determine XYZ
figure;
set(gcf,'Position',[100 500 1200 600]);
h = ft_plot_mesh(mesh,'EdgeColor','none',...
    'facealpha',1,'facecolor','skin'); camlight;
view([130 37]);
box on;
ax = gca;
ax.XLim = [-80 80];
ax.YLim = [-80 80];
ax.ZLim = [-80 80]; 
ax.Color = 'w';
t = hgtransform('Parent',ax);
set(h,'Parent',t);

pos_data_ds_k = [0 0 0];
indx = 1;

for k = 1:1000
    pos_data_ds_k(indx) = pos_data_ds_k(indx)+0.1;
    
    Txy     = makehgtform('translate', pos_data_ds_k);
    set(t,'Matrix',Txy);
    drawnow;
    pause(0.001);
    
end

%%
figure;
set(gcf,'Position',[100 500 1200 600]);
h = ft_plot_mesh(mesh,'facecolor',[238,206,179]./255,'EdgeColor','none',...
    'facealpha',0.8,'facecolor','skin'); camlight;
view([130 37]);
box on;
ax = gca;
ax.XLim = [-80 80];
ax.YLim = [-80 80];
ax.ZLim = [-80 80];
ax.Color = 'w';
t = hgtransform('Parent',ax);
set(h,'Parent',t);

deg_data_ds_k = [0 0 0];
indx = 3;

for k = 1:1000
    deg_data_ds_k(indx) = deg_data_ds_k(indx)+0.01;
    Rrotate = makehgtform('zrotate', deg_data_ds_k(indx));
    set(t,'Matrix',Rrotate);
    drawnow;
    pause(0.001);
end



%% Downsample the data to 10Hz
t_ds = downsample(([0:1:length(pos_data)-1]/Fs),100);
pos_data_ds = downsample(pos_data,100);
deg_data_ds = downsample(deg2rad(degXYZ),100);
deg_data_ds2 = downsample((degXYZ),100);

%% Plot Mesh
figure;
set(gcf,'Position',[100 500 1200 600]);
gif('animate_NA_movement.gif','DelayTime',0.1) 
subplot(2,2,[1 3]);
h = ft_plot_mesh(mesh,'EdgeColor','none',...
    'facealpha',1,'facecolor','skin'); camlight;
view([-130 20]);
box on;
ax = gca;
ax.XLim = [-50 50];
ax.YLim = [-50 50];
ax.ZLim = [-50 50];
ax.Color = 'w';
ylabel('Forward-Back (cm)');xlabel('Left-Right (cm)'); zlabel('Up-Down (cm)');
grid on
hold on;
%set(gca,'BoxStyle','full')

%plotcube(L,O,.1,[0.1 0.1 0.1]);   % use function plotcube
t = hgtransform('Parent',ax);
set(h,'Parent',t)

view_angle = 0;
subplot(2,2,2);
cols = [235 123 1; 76 0 230; 23 230 55]./255;
rrr = animatedline('Color',cols(1,:),'LineWidth',2);
bbb = animatedline('Color',cols(2,:),'LineWidth',2);
ggg = animatedline('Color',cols(3,:),'LineWidth',2);
xlim([0 t_ds(1)+30]);
xlabel('Time (s)');
ylim([-max(abs(pos_data_ds(:))) max(abs(pos_data_ds(:)))]);
ylabel('Position (cm)');
% legend({'Right-Left';'Forward-Back';'Up-Down'}, 'Location','NorthEast');
% legend('boxoff')


subplot(2,2,4);
cols = [235 210 0; 230 0 99; 23 196 230]/255;
rrr_r = animatedline('Color',cols(1,:),'LineWidth',2);
bbb_r = animatedline('Color',cols(2,:),'LineWidth',2);
ggg_r = animatedline('Color',cols(3,:),'LineWidth',2);
xlim([0 t_ds(1)+30]);
xlabel('Time (s)');
ylim([-max(abs(deg_data_ds2(:))) max(abs(deg_data_ds2(:)))]);
ylabel('Degree (°)');
% legend({'Pitch','Yaw','Roll'}, 'Location','NorthEast');
% legend('boxoff')

for k = 1:200
    subplot(2,2,[1 3]);
    pos_data_ds_k = [pos_data_ds(k,1) pos_data_ds(k,2) pos_data_ds(k,3)];
    deg_data_ds_k = [deg_data_ds(k,1) deg_data_ds(k,2) deg_data_ds(k,3)];
    Xrotate = makehgtform('xrotate', deg_data_ds_k(1));
    Yrotate = makehgtform('yrotate', deg_data_ds_k(2));
    Zrotate = makehgtform('zrotate', deg_data_ds_k(3));
    Txy     = makehgtform('translate', pos_data_ds_k);
    set(t,'Matrix',Xrotate*Yrotate*Zrotate*Txy);
    %set(t,'Matrix',Xrotate);
    title([num2str(t_ds(k)) 's']);

    subplot(2,2,2)
    addpoints(rrr,t_ds(k),pos_data_ds(k,1));
    addpoints(bbb,t_ds(k),pos_data_ds(k,2)); hold on
    addpoints(ggg,t_ds(k),pos_data_ds(k,3)); hold on;
    if t_ds(k) > 15
        xlim([t_ds(k)-15 t_ds(k)+15]);
    end

    subplot(2,2,4)
    addpoints(rrr_r,t_ds(k),deg_data_ds2(k,1));
    addpoints(bbb_r,t_ds(k),deg_data_ds2(k,2)); hold on
    addpoints(ggg_r,t_ds(k),deg_data_ds2(k,3)); hold on;
    if t_ds(k) > 15
        xlim([t_ds(k)-15 t_ds(k)+15]);
    end
    drawnow;
    gif;
end






%%
%%
figure;
view(3)

[x,y,z] = cylinder([.2 0]);
h(1) = surface(x,y,z,'FaceColor','red');

t = hgtransform('Parent',ax);
set(h,'Parent',t)

pos1 = [pos_data_ds(1,1) pos_data_ds(2,1) pos_data_ds(3,1)]; 

for k = 1:length(t_ds)
    pos_data_ds_k = [pos_data_ds(k,1) pos_data_ds(k,1) pos_data_ds(k,1)];
    Txy = makehgtform('translate',[pos1-pos_data_ds_k]);
    set(t,'Matrix',Txy)
    title([num2str(t_ds(k)) 's']);
    drawnow;
    pause(0.001);
end

pause(1)

%%
figure;
gif('animate_NA_movement.gif','DelayTime',0.1) 

hold on;
view([-63 20]);
    ylim([-200 200]);
    xlim([-200 200]);
    zlim([-200 200]);
plotcube(L,O,.1,[0.1 0.1 0.1]);   % use function plotcube
camlight;

v = 0;

h = animatedline('MaximumNumPoints', 30,'Color','r','LineWidth',2);
for k = 1:length(pos_data_ds(:,3))
    addpoints(h,pos_data_ds(k,3),pos_data_ds(k,2),pos_data_ds(k,1));
    title([num2str(t_ds(k)) 's']);
    drawnow;
    pause(0.01);
    view([v 20]);
    if v == 359.5
        v = 0;
    else
        v = v+0.5;
    end
    gif
end






