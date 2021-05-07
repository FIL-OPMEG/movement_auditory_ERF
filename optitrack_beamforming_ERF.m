function optitrack_beamforming_ERF(save_dir, atlas_dir,...
    run_num,mri, headmodel,sourcemodel)

%% Hardcoded for now
compute_VE      = 'single_points';

%% Load data
cd(save_dir);
load(['data_run' num2str(run_num) '.mat']);

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

% %%
% % the data consists of fewer channels than the precomputed
% % leadfields, the following chunk of code takes care of this
% [a,b] = match_str(data.label, lf.label);
% for k = 1:numel(lf.leadfield)
%     if ~isempty(lf.leadfield{k})
%         tmp = lf.leadfield{k};
%         tmp = tmp(b,:);
%         tmp = tmp-repmat(mean(tmp,1),[size(tmp,1) 1]); % average re-ref
%         lf.leadfield{k} = tmp;
%     end
% end
% lf.label = lf.label(b);

%% Compute covariance matrix
cfg                  = [];
cfg.covariance       = 'yes';
cfg.vartrllength     = 2;
cfg.covariancewindow = [0 0.5];
avg                  = ft_timelockanalysis(cfg,data);

%% Source Analysis
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
cfg.lcmv.lambda        = '0.1%';
sourceall              = ft_sourceanalysis(cfg, avg);
% % Remove extra .mom
% for k = 1:size(sourceall.pos,1)
%     if ~isempty(sourceall.avg.mom{k})
%         sourceall.avg.mom{k}(4:6,:) = [];
%     end
% end

%%
% Replace .pos field with template_grid.pos

[t, r] = ft_version;
ddd = load(fullfile(r,'template/sourcemodel/standard_sourcemodel3d5mm.mat'));
template_grid = ddd.sourcemodel;
clear ddd
template_grid = ft_convert_units(template_grid,'mm');

sourceall.pos = template_grid.pos;

%%
% Remove cfg field to save memory
sourceall = rmfield(sourceall,'cfg');

%%
ERF_name = {'M100'};
ERF_toi  = [0.08 0.12]

for ERF = 1:length(ERF_name)
    
    source_pow_post = get_source_pow(data,sourceall,[ERF_toi(ERF,1) ERF_toi(ERF,2)]);
    source_pow_pre  = get_source_pow(data,sourceall,[-0.08 -0.04]);
    % %
    % % source_pow_post.avg.pow = source_pow_post.avg.pow./source_pow_post.avg.noise;
    % % source_pow_pre.avg.pow = source_pow_pre.avg.pow./source_pow_pre.avg.noise;
    % %
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
    
    % Mask bits outside the brain
    %sourceI.anat_mask = spm_brain_seg.brain .* double(sourceI.anatomy);
    
    % Plot
    cfg                 = [];
    cfg.funparameter    = 'pow';
    cfg.funcolormap        = cmap;
    cfg.funcolorlim     = 'maxabs';
    %cfg.maskparameter   = 'anat_mask';
    ft_sourceplot(cfg,sourceI);
    title([ERF_name{ERF}]);
    
    %% Export to nifti formt and use your favourite MRI software to visualise
    cd(save_dir);
    cfg = [];
    cfg.filetype = 'nifti';
    cfg.filename = [ERF_name{ERF} '_run' num2str(run_num)];
    cfg.parameter = 'pow';
    ft_sourcewrite(cfg,sourceI);
%     
%     
%         %% Export to connectome workbench (specfic to my computer)
%         try
%             system(['C:\wtcnapps\workbench\bin_windows64\wb_command -volume-to-surface-mapping D:\data\auditory_moving_ERF\RS_results\' ...
%                 ['M100_run' num2str(run_num) '.nii'] ' D:\scripts\Conte69_atlas-v2.LR.32k_fs_LR.wb\Conte69.L.midthickness.32k_fs_LR.surf.gii D:\data\auditory_moving_ERF\RS_results\' ['M100_run' num2str(run_num)] '_LEFT.shape.gii -trilinear'])
%             system(['C:\wtcnapps\workbench\bin_windows64\wb_command -volume-to-surface-mapping D:\data\auditory_moving_ERF\RS_results\' ...
%                 ['M100_run' num2str(run_num) '.nii'] ' D:\scripts\Conte69_atlas-v2.LR.32k_fs_LR.wb\Conte69.R.midthickness.32k_fs_LR.surf.gii D:\data\auditory_moving_ERF\RS_results\' ['M100_run' num2str(run_num)] '_RIGHT.shape.gii -trilinear'])
%         catch
%             disp('Could not convert to gifti format');
%         end
end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VE analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch compute_VE
    case 'parcel'
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Now the VE analysis: Parcel
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
        %atlas_HCPMMP.coordsys          = 'mni';
        
        %% Do beamforming without NAI weight normalisation
        disp('Beamforming...');
        cfg                 = [];
        cfg.channel         = data.label;
        cfg.method          = 'lcmv';
        cfg.grid            = lf;
        cfg.headmodel       = headmodel;
        cfg.keeptrials      = 'no';
        cfg.lcmv.keepfilter = 'yes';
        cfg.lcmv.fixedori   = 'yes';
        cfg.lcmv.lambda     = '0.1%';
        cfg.lcmv.projectnoise  = 'yes';
        sourceavg           = ft_sourceanalysis(cfg, avg);
        sourceavg           = rmfield(sourceavg,'cfg');
        sourceavg.pos       = template_grid.pos;
        
        %% Create VE
        [VE] = atlas2VE(atlas_HCPMMP,template_grid,...
            {'LH_Early_Auditory_Cortex'}, data,...
            sourceavg, avg);
        
        % Select only the Early Auditory Cortex Parcel
        cfg             = [];
        cfg.channel     = {'LH_Early_Auditory_Cortex'};
        VE_A1  = ft_selectdata(cfg,VE);
        
        % Save the data
        save(['VE_A1_run_' num2str(run_num)],'VE_A1');
        
        %% Perform timelockanalysis
        cfg             = [];
        avg_VE      	= ft_timelockanalysis([],VE_A1);
        
        %%
        epoched_dataset = [];
        
        for i = 1:length(data.trial)
            epoched_dataset(:,:,i) = VE_A1.trial{1,i};
        end
        
        SE = std(epoched_dataset,[],3)/sqrt(size(epoched_dataset,3));
        avg_VE.t_value = avg_VE.avg./SE;
        
        %% Plot
        cfg = [];
        cfg.channel = avg_VE.label;
        cfg.parameter = 't_value';
        cfg.baseline = [-0.1 0];
        cfg.showlegend    = 'yes';
        cfg.xlim    = [-0.1 0.4];
        cfg.linecolor = 'k';
        cfg.linewidth = 2;
        %cfg.ylim = [-2.8e-7 -2.8e-7];
        figure; ft_singleplotER(cfg,avg_VE)
        set(gca,'FontSize',18);
        xlabel('Time (s)','FontSize',20);
        ylabel('t-value','FontSize',20);
        title('');
        print(['run_' num2str(run_num) '_VE'],'-dpng','-r300');
        
        %% Plot all on the same graph
        if run_num == 3
            figure;
            avg_VE = [];
            for r = 1:3
                load(['VE_A1_run_' num2str(r) '.mat']);
                
                %% Perform timelockanalysis
                cfg             = [];
                avg_VE{r} = ft_timelockanalysis([],VE_A1);
                
                % Convert to t-value
                epoched_dataset = [];
                
                for i = 1:length(VE_A1.trial)
                    epoched_dataset(:,:,i) = VE_A1.trial{1,i};
                end
                
                SE = std(epoched_dataset,[],3)/sqrt(size(epoched_dataset,3));
                avg_VE{r}.t_value = avg_VE{r}.avg./SE;
                clear epoched_dataset
                subplot(1,3,r); imagesc(VE_A1.time{1},1:length(VE_A1.trial),vertcat(VE_A1.trial{:}));
                zlim([-8e-27 8e-27]);
                title(['Run ' num2str(r)]);
            end
            
            %% Plot
            figure;
            cfg = [];
            cfg.channel = avg_VE{1}.label;
            cfg.parameter = 't_value';
            cfg.baseline = [-0.1 0];
            cfg.showlegend    = 'no';
            cfg.xlim    = [-0.1 0.4];
            cfg.linecolor = [0.4275    0.9804    0.3922;0.9804    0.5686    0.3843;
                0.2824    0.3137    0.9804];
            cfg.linewidth = 2;
            cfg.ylim = [-10 10];
            ft_singleplotER(cfg,avg_VE{1},avg_VE{2},avg_VE{3});
            set(gca,'FontSize',18);
            xlabel('Time (s)','FontSize',20);
            ylabel('t-value','FontSize',20);
            title('');
            print('run123_VE_ERF','-dpng','-r300');
            legend({'Sitting';'Standing';'Standing + Moving'},'Location','SouthOutside');
            
        end
        
        
        
    case 'single_points'
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Alternatively we can compute the VE using two single points
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% Load Source stuff
        %  LF for Virtual electrode analysis
        cfg             = [];
        cfg.template    = mri;
        cfg.nonlinear   = 'yes';
        norm            = ft_volumenormalise([],mri);
        
        % Auditory Cortex
        pos = [-48 -22 4; 52 -22 8];
        
        % Now we warp the MNI coordinates using the nonlinear warping parameters
        posback         = ft_warp_apply(norm.params,pos,'sn2individual');
        % xyz positions in individual coordinates
        pos_grid        = ft_warp_apply(pinv(norm.initial),posback);
        
        figure; ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:));
        ft_plot_mesh(pos_grid,'vertexcolor','r');
        
        % % Find location of closest vertex on cortical mesh
        % Idx = knnsearch(lf.pos, pos_grid)
        
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
        
        %%
        cfg                    = [];
        cfg.channel            = data.label;
        cfg.grad               = data.grad;
        cfg.method             = 'lcmv';
        cfg.grid               = lf_2;
        cfg.headmodel          = headmodel;
        cfg.lcmv.keepfilter    = 'yes';
        cfg.lcmv.fixedori      = 'yes';
        cfg.lcmv.projectnoise  = 'yes';
        %cfg.lcmv.weightnorm    = 'nai';
        cfg.lcmv.lambda        = '1%';
        sourceall              = ft_sourceanalysis(cfg, avg);
        
        % Find filter from Idx point
        filter123 = cat(1,sourceall.avg.filter{1,:});
        
        VE          = [];
        VE.label    = {'A1'};
        VE.fsample  = data.fsample;
        for subs=1:size(data.trial,2)
            % note that this is the non-filtered "raw" data
            VE.time{subs}       = data.time{subs};
            VE.trial{subs}(:,:) = filter123(:,:)*data.trial{subs}(:,:);
        end
        
        %% Save the data
        save(['VE_A1_run_' num2str(run_num)],'VE');
        
        %% Perform timelockanalysis
        cfg             = [];
        avg_VE      	= ft_timelockanalysis([],VE);
        
        %%
        epoched_dataset = [];
        
        for i = 1:length(data.trial)
            epoched_dataset(:,:,i) = VE.trial{1,i};
        end
        
        SE = std(epoched_dataset,[],3)/sqrt(size(epoched_dataset,3));
        avg_VE.t_value = avg_VE.avg./SE;
        
        %% Plot
        cfg = [];
        cfg.channel = avg_VE.label;
        cfg.parameter = 't_value';
        cfg.baseline = [-0.1 0];
        cfg.showlegend    = 'yes';
        cfg.xlim    = [-0.1 0.4];
        cfg.linecolor = 'k';
        cfg.linewidth = 2;
        cfg.ylim = [-10 10];
        figure; ft_singleplotER(cfg,avg_VE)
        set(gca,'FontSize',18);
        xlabel('Time (s)','FontSize',20);
        ylabel('t-value','FontSize',20);
        title('');
        print(['run_' num2str(run_num) '_VE'],'-dpng','-r300');
        
        
        %% Plot all the data on the same axis
        if run_num == 3
            figure;
            avg_VE = [];
            for r = 1:3
                load(['VE_A1_run_' num2str(r) '.mat']);
                
                %% Perform timelockanalysis
                cfg             = [];
                avg_VE{r} = ft_timelockanalysis([],VE);
                
                % Convert to t-value
                epoched_dataset = [];
                
                for i = 1:length(VE.trial)
                    epoched_dataset(:,:,i) = VE.trial{1,i};
                end
                
                SE = std(epoched_dataset,[],3)/sqrt(size(epoched_dataset,3));
                avg_VE{r}.t_value = (avg_VE{r}.avg./SE);
                clear epoched_dataset
                subplot(1,3,r); imagesc(VE.time{1},1:length(VE.trial),vertcat(VE.trial{:}));
                zlim([-1e-27 1e-27]);
                title(['Run ' num2str(r)]);
            end
            
            %% Plot
            figure;
            cfg = [];
            cfg.channel = VE.label;
            cfg.parameter = 't_value';
            cfg.baseline = [-0.1 0];
            cfg.showlegend    = 'no';
            cfg.xlim    = [-0.1 0.4];
            cfg.linecolor = [0.4275    0.9804    0.3922;0.9804    0.5686    0.3843;
                0.2824    0.3137    0.9804];
            cfg.linewidth = 2;
            cfg.ylim = [-10 10];
            ft_singleplotER(cfg,avg_VE{1},avg_VE{2},avg_VE{3});
            set(gca,'FontSize',18);
            xlabel('Time (s)','FontSize',20);
            ylabel('t-value','FontSize',20);
            title('');
            print('run123_VE_ERF','-dpng','-r300');
            legend({'Sitting';'Standing';'Standing + Moving'},'Location','SouthOutside');
            
        end
end

end







