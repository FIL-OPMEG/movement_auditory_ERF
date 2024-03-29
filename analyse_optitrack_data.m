function analyse_optitrack_data(save_dir)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for analysing and plotting the motion capture data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
Fs  = 1000; % Sampling Rate

%% Load Data
% Bit fudged
cd(save_dir);
disp('Loading data...');
m               = load('MovementDataOut_run3.mat');
MovementDataOut = m.MovementDataOut;
clear m

m               = load('trl_index.mat');
trl_index       = m.trl_index;
clear m

m               = load('trial2keep.mat');
trial2keep      = m.trial2keep;
clear m

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Position Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot the position data
pos_data = MovementDataOut.rigidbodies.data(:,4:6); % Convert to cm
pos_data = pos_data(trl_index(1):trl_index(end)+500,:); % Trim to start of the tones

% Convert from Y-up to Z-up right-handed coordinate system
pos_data = ft_warp_apply([1 0 0 0; 0 0 -1 0; 0 1 0 0 ; 0 0 0 1], pos_data);

pos_data(:,2) = pos_data(:,2).*-1;
pos_data(:,1) = pos_data(:,1).*-1;
% Here we are assuming:

% - X = Left-Right
% - Y = Forward-Back
% - Z = Up-Down

%%
% Make the first point 0 0 0
pos_data = zero_optitrack_data(pos_data);

% Make time array starting at 0
t = [0:1:length(pos_data)-1]/Fs;

%% Plot
coord = {'X','Y','Z'};
cols = [100, 100, 250; 142 185 57; 198 61 61]./255;
coord_name = {'X: Left-Right','Y:Forward-Back','Z: Up-Down'};

figure;
set(gcf,'Position',[1 1 1200 500]);

for c = 1:size(pos_data,2)
    m = pos_data(:,c);
    p1 = plot(t,m,'LineWidth',2,'Color',cols(c,:)); hold on;
    p1.Color(4) = 0.7; 
end
ylabel(['Distance (cm)'],'FontSize',24);
set(gca,'FontSize',18);
xlabel('Time (s)','FontSize',24);
legend(coord_name,'Location','EastOutside');
print('opti_pos_reject','-dsvg','-r300');

%% Get euclidian distance from start point and plot
dist_from_start = pdist2(pos_data,pos_data(1,:));

% Plot the euclidian distance over time
figure;
set(gcf,'Position',[1 1 1200 600]);
m = dist_from_start(:,:);
p1 = plot(t,m,'LineWidth',2,'Color','k'); hold on;
p1.Color(4) = 0.8;
set(gca,'FontSize',18);
ylabel(['Euclidian Distane (cm)'],'FontSize',22);
xlabel('Time (s)','FontSize',22);
print('euclidian_distance_from_start','-dsvg','-r300');

%% Plot some histograms

% For X, Y, Z...
for c = 1:3
    % Use make_pos_hist
    make_pos_hist(pos_data(:,c),cols(c,:));
    
    % Change the Y ticks
    yt = get(gca, 'YTick');
    set(gca,'YTick',linspace(yt(1),yt(end),4));
    yt = get(gca, 'YTick');
    
    % Change X lim
    if c == 2
        xlim([-100 100]);
    elseif c == 1
        xlim([-100 100]);
    else
        xlim([-75 25]);
    end
    
    % Change X ticks
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
% - Y = Roll
% - Z = Yaw
degXYZ = ft_warp_apply([1 0 0 0; 0 0 -1 0; 0 1 0 0 ; 0 0 0 1], degXYZ);

degXYZ(:,2) = degXYZ(:,2).*-1;
degXYZ(:,1) = degXYZ(:,1).*-1;

%%
% Make the first point 0 0 0
degXYZ = zero_optitrack_data(degXYZ);

%% Interpolate angles that change abruptly

for ang = 1:3
    
    rot2= degXYZ(:,ang);
    [TF,S1] = ischange(rot2,'Threshold',500);
    indx = find(TF==1);
    mm      = movmedian(rot2,1000);
    
    
    for i = 1:length(indx)
        try
            rot2(indx(i)-200:indx(i)+200) = mm(indx(i)-200:indx(i)+200);
        catch
            rot2(indx(i)-200:end) = mm(indx(i)-200:end);
        end
    end
    degXYZ(:,ang) = rot2;
end

%% Plot the Rotation Data
cols = [235 210 0; 230 0 99; 23 196 230]/255;

coord = {'X','Y','Z'};
coord_name = {'X: Pitch','Y: Roll','Z: Yaw'};

figure;
set(gcf,'Position',[1 1 1200 500]);

for c = 1:size(degXYZ,2)
    m = degXYZ(:,c);
    %m(log_array,:) = NaN;
    p1 = plot(t,m,'LineWidth',2,'Color',cols(c,:)); hold on;
    p1.Color(4) = 0.7; 
end
%plot(t,m,'r','LineWidth',2);
ylabel(['Degrees (�)'],'FontSize',24);
set(gca,'FontSize',18);
xlabel('Time (s)','FontSize',24);
legend(coord_name,'Location','EastOutside');
print('opti_rot_reject','-dsvg','-r300');

%% Polar Histogram
for c = 1:3
    figure; ax = polaraxes;
    h = polarhistogram(deg2rad(degXYZ(:,c)),50,...
        'Normalization','probability');
    thetalim([-60 60]);
    ax.FontSizeMode = 'manual'
    ax.FontSize = 18
    ax.RAxis.FontSize = 14;
    ax.ThetaTickLabel = {'-60�','-30�','0�','30�','60�'}
    
    h(1).FaceColor = cols(c,:);

    print(['polarhistogram' num2str(c)],'-dpng','-r300');
end

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Data per trial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trl_index_for_euc = trl_index-trl_index(1);
trl_index_for_euc(1) = 1;
%trl_index = trl_index(trial2keep);

keep_cm = [];

% Position
for k = 1:length(trial2keep)
    pos_data_t = pos_data(trl_index_for_euc(k):trl_index_for_euc(k)+500,:);
    for d = 1:3
        distances = pdist2(pos_data_t(:,d), mean(pos_data_t(:,d),2));
        % Find the max distance
        keep_cm(k,d) = max(distances(:));
    end
end

% Degrees
for k = 1:length(trial2keep)
    deg_data_t = degXYZ(trl_index_for_euc(k):trl_index_for_euc(k)+500,:);
    for d = 1:3
        distances = pdist2(deg_data_t(:,d), mean(deg_data_t(:,d),2));
        % Find the max distance
        keep_cm(k,d+3) = max(distances(:));
    end
end

% Euclidian Distance
for k = 1:length(trial2keep)
    euclidian_data_t = dist_from_start(trl_index_for_euc(k):trl_index_for_euc(k)+500);
    distances = pdist2(euclidian_data_t, euclidian_data_t);
    % Find the max distance
    keep_cm(k,7) = max(distances(:));
end

%% Plot Boxplot for Quick View
figure;
boxplot(keep_cm)

%% Export to csv file for plotting with Python
t = array2table(keep_cm,'VariableNames',{'PosX','PosY','PosZ','RotX',...
'RotY','RotZ','EuclidianDist'});
writetable(t,'trial_movt.csv','Delimiter',',')
















%% Make the Mesh
% This will only work on my computer
load('D:\tess_head.mat');
mesh = [];
mesh.pos = Vertices;
mesh.tri = Faces;
mesh.unit = 'm';
mesh.coordsys = 'mni';
mesh = ft_convert_units(mesh,'cm');

% Rotate by 90deg to get into MATLAB coordinate system
ttt = cos(90*(pi/180));
rrr = sin(90*(pi/180));

trans = [ttt -rrr 0 0;
    rrr ttt 0 0;
    0 0 1 0;
    0 0 0 1];

mesh.pos = ft_warp_apply(trans,mesh.pos);

% Plot
figure; ft_plot_mesh(mesh); camlight;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Let's make an animated plot!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Downsample the data to 10Hz
t_ds = downsample(([0:1:length(pos_data)-1]/Fs),100);
pos_data_ds = downsample(pos_data,100);
deg_data_ds = downsample(deg2rad(degXYZ),100);
deg_data_ds2 = downsample((degXYZ),100);

figure; plot(deg_data_ds2);

%% Plot Mesh
figure;
set(gcf,'Position',[100 500 1200 600]);
gif('example_movt_video.gif','frame',gcf,'DelayTime',0.1) 
subplot(2,2,[1 3]);
h = ft_plot_mesh(mesh,'EdgeColor','none',...
    'facealpha',1,'facecolor','skin'); camlight;
view([-130 20]);
box on;
ax = gca;
ax.XLim = [-100 100];
ax.YLim = [-100 100];
ax.ZLim = [-100 100];
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
cols = [100, 100, 250; 142 185 57; 198 61 61]./255;
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
ylabel('Degree (�)');
% legend({'Pitch','Yaw','Roll'}, 'Location','NorthEast');
% legend('boxoff')

for k = 1401:1801
    subplot(2,2,[1 3]);
    pos_data_ds_k = [pos_data_ds(k,1) pos_data_ds(k,2) pos_data_ds(k,3)];
    deg_data_ds_k = [deg_data_ds(k,1) deg_data_ds(k,2) deg_data_ds(k,3)];
    Xrotate = makehgtform('xrotate', deg_data_ds_k(1));
    Yrotate = makehgtform('yrotate', deg_data_ds_k(2));
    Zrotate = makehgtform('zrotate', deg_data_ds_k(3));
    Txy     = makehgtform('translate', pos_data_ds_k);
    set(t,'Matrix',Xrotate*Yrotate*Zrotate*Txy);
    %set(t,'Matrix',Zrotate);    
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


%% SPARE BITS OF CODE:

% %% Bit of code to determine XYZ
% figure;
% set(gcf,'Position',[100 500 1200 600]);
% h = ft_plot_mesh(mesh,'EdgeColor','none',...
%     'facealpha',1,'facecolor','skin'); camlight;
% view([130 37]);
% box on;
% ax = gca;
% ax.XLim = [-80 80];
% ax.YLim = [-80 80];
% ax.ZLim = [-80 80]; 
% ax.Color = 'w';
% t = hgtransform('Parent',ax);
% set(h,'Parent',t);
% 
% pos_data_ds_k = [0 0 0];
% indx = 1;
% 
% for k = 1:1000
%     pos_data_ds_k(indx) = pos_data_ds_k(indx)+0.1;
%     
%     Txy     = makehgtform('translate', pos_data_ds_k);
%     set(t,'Matrix',Txy);
%     drawnow;
%     pause(0.001);
%     
% end
% 
% %%
% figure;
% set(gcf,'Position',[100 500 1200 600]);
% h = ft_plot_mesh(mesh,'facecolor',[238,206,179]./255,'EdgeColor','none',...
%     'facealpha',0.8,'facecolor','skin'); camlight;
% view([130 37]);
% box on;
% ax = gca;
% ax.XLim = [-80 80];
% ax.YLim = [-80 80];
% ax.ZLim = [-80 80];
% ax.Color = 'w';
% t = hgtransform('Parent',ax);
% set(h,'Parent',t);
% 
% deg_data_ds_k = [0 0 0];
% indx = 3;
% 
% for k = 1:1000
%     deg_data_ds_k(indx) = deg_data_ds_k(indx)+0.01;
%     Rrotate = makehgtform('zrotate', deg_data_ds_k(indx));
%     set(t,'Matrix',Rrotate);
%     drawnow;
%     pause(0.001);
% end



