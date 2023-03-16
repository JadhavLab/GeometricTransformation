%---------------------------------------------------------------%
%  This is the main script for calculating the distance of      %
%  neural states between N' and F neural manifolds.             %
%  The UMAP MATLAB toolbox (V4.1) is needed (https://www.       %
%  mathworks.com/matlabcentral/fileexchange/                    %
%  71902-uniform-manifold-approximation-and-projection-umap).   %                         %
%  -- Wenbo Tang (Jan 07, 2023)                                 %
%---------------------------------------------------------------%
clc
clear all
close all
%%
animalprefix_list = {'AM2','JS17','JS21','ZT2','JS34'}; % all animals
famnov_dir = '/Volumes/NovelFamiliar'; % set directory
day = 1;
ep = 6;%epoch, epoch 6 = N' session
CTXHP = 1;% CA1 = 0 , PFC = 1;
pos_interp = (0:0.005:1)';% position bins, normalized
cellthresh = 5;% minimal number of active cells per bin needed
%%
% calculate the trial-averaged firing rates
% the code to get UMAP manifold is the same as in "cal_UMAP.m"
% get trial number for each animal
trialnum_all = [];
for theAnimal = 1:length(animalprefix_list)
    animalprefix = animalprefix_list{theAnimal};
    dir = sprintf('%s/%s_direct/',famnov_dir,animalprefix); % set directory
    %%
    trajinfo = loaddatastruct(dir, animalprefix, 'trajinfo', day); % get trajinfo
    %---------trajtime--------%
    trajbound = trajinfo{day}{ep}.trajbound;
    trajtime_raw = trajinfo{day}{ep}.trajtime;
    wellstand = trajinfo{day}{ep}.wellstend;
    traj_correct = trajinfo{day}{ep}.rewarded;


    trajtime{1} = trajtime_raw(find(wellstand(:,1)== 1 & wellstand(:,2)== 2),:);
    trajtime{2} = trajtime_raw(find(wellstand(:,1)== 2 & wellstand(:,2)== 1),:);
    trajtime{3} = trajtime_raw(find(wellstand(:,1)== 1 & wellstand(:,2)== 3),:);
    trajtime{4} = trajtime_raw(find(wellstand(:,1)== 3 & wellstand(:,2)== 1),:);

    outtraj_left = length(find(wellstand(:,1)== 1 & wellstand(:,2)== 2));% outbound: center to left
    outtraj_right = length(find(wellstand(:,1)== 1 & wellstand(:,2)== 3));% outbound: center to right
    intraj_left = length(find(wellstand(:,1)== 2 & wellstand(:,2)== 1));
    intraj_right = length(find(wellstand(:,1)== 3 & wellstand(:,2)== 1));

    trialnum_all = [trialnum_all;outtraj_left,intraj_left,outtraj_right,intraj_right];
end
trialnum = min(trialnum_all);
%%
% get firing rates for each trial
spike_matrix = [];
spike_matrix_label = [];
for theAnimal = 1:length(animalprefix_list)
    animalprefix = animalprefix_list{theAnimal};
    dir = sprintf('%s/%s_direct/',famnov_dir,animalprefix); % set directory
    % interneuron list
    if strcmp(animalprefix,'AM2')
        exclude_list = [6,9];%AM2
    elseif strcmp(animalprefix,'JS21')
        exclude_list = [6,3;6,6;6,7;6,8;7,1;21,1;22,1;25,1];%JS21
    elseif strcmp(animalprefix,'JS17')
        exclude_list = [6,1;7,1];%JS17
    elseif strcmp(animalprefix,'ZT2')
        exclude_list = [31,1;31,5;32,2;27,1;14,10];%ZT2
    elseif strcmp(animalprefix,'JS34')
        exclude_list = [];%JS34
    end
    % files for single-trial firing rates can also be found in ./SingleTrial_ratemaps/
    load(sprintf('%s%slinfields_singletrial0%d.mat',dir,animalprefix,day)); % get linearized place fields
    %%
    spike_matrix_animal = [];
    spike_matrix_label = [];
    %-----get cell indices-----%
    if CTXHP
        [cellidx, ~] = matchidx_acrossep(dir, animalprefix, day,exclude_list); %(tet, cell)
    else
        [~, cellidx] = matchidx_acrossep(dir, animalprefix, day,exclude_list); %(tet, cell)
    end
    % for analyzing CA1 cells, simply replace ctxidx with hpidx below
    cellnum = length(cellidx(:,1));
    for tracknum = 1:4 % 4 different trajctories
            temp_label = pos_interp + tracknum -1;
            temp = zeros(length(pos_interp),cellnum,trialnum(tracknum)); % reset
            for tj = 1:trialnum(tracknum)
                temp_traj = zeros(length(pos_interp),cellnum); % reset

                for i = 1:cellnum
                    cind = cellidx(i,:);
                    if (length(linfields{day}{ep})>= cind(1))
                        if  (length(linfields{day}{ep}{cind(1)})>= cind(2))
                            linfield1 = linfields{day}{ep}{cind(1)}{cind(2)};
                            if length(linfield1)>= tracknum
                                temp1 = linfield1{tracknum}{tj};
                                if ~isempty(temp1)
                                    pos1 = temp1(:,1);
                                    linfield_hp = temp1(:,5);% occupancy normalized firing rate 
                                    pos1 = pos1./max(pos1);
                                    % make sure the two linearized fields are the
                                    % same length
                                    linfield_hp(isnan(linfield_hp)) = 0;
                                    linfield_hp = interp1(pos1,linfield_hp,pos_interp,'nearest');

                                else
                                    linfield_hp = zeros(size(pos_interp));
                                end
                            else
                                linfield_hp = zeros(size(pos_interp));
                            end
                        else
                            linfield_hp = zeros(size(pos_interp));
                        end
                    else
                        linfield_hp = zeros(size(pos_interp));
                    end
                    temp_traj(:,i) = linfield_hp';
                end
                temp(:,:,tj) = temp_traj;
            end
            temp = nanmean(temp,3);
            spike_matrix_animal = [spike_matrix_animal;temp];
            spike_matrix_label = [spike_matrix_label;temp_label];
    end
    spike_matrix = [spike_matrix,spike_matrix_animal];
end
%%
% exclude bins with less than 5 cells active
validid = [];
for i = 1:length(spike_matrix(:,1))
    temp = spike_matrix(i,:);
    id = find(temp > 0);
    if length(id) >= cellthresh
        validid = [validid;i];
    end
end

spike_matrix = spike_matrix(validid,:);
spike_matrix_label = spike_matrix_label(validid);
%%
% remove the file folder that has functions with the same name as those in the UMAP toolbox
rmpath(genpath('/Users/wenbotang/Src_Matlab'))
% use template
if CTXHP
    spike_umap = run_umap(spike_matrix,'min_dist',0.6,'n_neighbors',50,'metric','cosine','n_components',3,'template_file','all_PFCumap_trials_ep6.mat');
else
    spike_umap = run_umap(spike_matrix,'min_dist',0.6,'n_neighbors',50,'metric','cosine','n_components',3,'template_file','all_CA1umap_trials_ep6.mat');
end
% UMAP transformed trajectories
spike_umap_avg = spike_umap;
spike_umap_avg_label = spike_matrix_label;
%%
% add back the file folder needed
addpath(genpath('/Users/wenbotang/Src_Matlab'))
% gather data from the F session
ep = 4; %epoch 4, F
trialnum_all = [];
for theAnimal = 1:length(animalprefix_list)
    animalprefix = animalprefix_list{theAnimal};
    dir = sprintf('%s/%s_direct/',famnov_dir,animalprefix); % set directory
    %%
    trajinfo = loaddatastruct(dir, animalprefix, 'trajinfo', day); % get trajinfo
    %---------trajtime--------%
    trajbound = trajinfo{day}{ep}.trajbound;
    trajtime_raw = trajinfo{day}{ep}.trajtime;
    wellstand = trajinfo{day}{ep}.wellstend;
    traj_correct = trajinfo{day}{ep}.rewarded;


    trajtime{1} = trajtime_raw(find(wellstand(:,1)== 1 & wellstand(:,2)== 2),:);
    trajtime{2} = trajtime_raw(find(wellstand(:,1)== 2 & wellstand(:,2)== 1),:);
    trajtime{3} = trajtime_raw(find(wellstand(:,1)== 1 & wellstand(:,2)== 3),:);
    trajtime{4} = trajtime_raw(find(wellstand(:,1)== 3 & wellstand(:,2)== 1),:);

    outtraj_left = length(find(wellstand(:,1)== 1 & wellstand(:,2)== 2));% outbound: center to left
    outtraj_right = length(find(wellstand(:,1)== 1 & wellstand(:,2)== 3));% outbound: center to right
    intraj_left = length(find(wellstand(:,1)== 2 & wellstand(:,2)== 1));
    intraj_right = length(find(wellstand(:,1)== 3 & wellstand(:,2)== 1));

    trialnum_all = [trialnum_all;outtraj_left,intraj_left,outtraj_right,intraj_right];
end
trialnum = min(trialnum_all);


spike_matrix = [];
spike_matrix_label = [];
for theAnimal = 1:length(animalprefix_list)
    animalprefix = animalprefix_list{theAnimal};
    dir = sprintf('%s/%s_direct/',famnov_dir,animalprefix); % set directory
    % interneuron list
    if strcmp(animalprefix,'AM2')
        exclude_list = [6,9];%AM2
    elseif strcmp(animalprefix,'JS21')
        exclude_list = [6,3;6,6;6,7;6,8;7,1;21,1;22,1;25,1];%JS21
    elseif strcmp(animalprefix,'JS17')
        exclude_list = [6,1;7,1];%JS17
    elseif strcmp(animalprefix,'ZT2')
        exclude_list = [31,1;31,5;32,2;27,1;14,10];%ZT2
    elseif strcmp(animalprefix,'JS34')
        exclude_list = [];%JS34
    end
    load(sprintf('%s%slinfields_singletrial0%d.mat',dir,animalprefix,day)); % get linearized place fields
    %%
    spike_matrix_animal = [];
    spike_matrix_label = [];
    %-----get cell indices-----%
    if CTXHP
        [cellidx, ~] = matchidx_acrossep(dir, animalprefix, day,exclude_list); %(tet, cell)
    else
        [~, cellidx] = matchidx_acrossep(dir, animalprefix, day,exclude_list); %(tet, cell)
    end
    cellnum = length(cellidx(:,1));
    for tracknum = 1:4
            for tj = 1:trialnum(tracknum)
                temp_label = pos_interp + tracknum -1;
                temp = zeros(length(pos_interp),cellnum); % reset
                for i = 1:cellnum
                    cind = cellidx(i,:);
                    if (length(linfields{day}{ep})>= cind(1))
                        if  (length(linfields{day}{ep}{cind(1)})>= cind(2))
                            linfield1 = linfields{day}{ep}{cind(1)}{cind(2)};
                            if length(linfield1)>= tracknum
                                temp1 = linfield1{tracknum}{tj};
                                if ~isempty(temp1)
                                    pos1 = temp1(:,1);
                                    linfield_hp = temp1(:,5);% occupancy normalized firing rate 
                                    pos1 = pos1./max(pos1);
                                    % make sure the two linearized fields are the
                                    % same length
                                    linfield_hp(isnan(linfield_hp)) = 0;
                                    linfield_hp = interp1(pos1,linfield_hp,pos_interp,'nearest');

                                else
                                    linfield_hp = zeros(size(pos_interp));
                                end
                            else
                                linfield_hp = zeros(size(pos_interp));
                            end
                        else
                            linfield_hp = zeros(size(pos_interp));
                        end
                    else
                        linfield_hp = zeros(size(pos_interp));
                    end
                    temp(:,i) = linfield_hp';
                end
                spike_matrix_animal = [spike_matrix_animal;temp];
                spike_matrix_label = [spike_matrix_label;temp_label];
            end
            
    end
    spike_matrix = [spike_matrix,spike_matrix_animal];
end
%%
% exclude bins with less than 5 cells active
validid = [];
for i = 1:length(spike_matrix(:,1))
    temp = spike_matrix(i,:);
    id = find(temp > 0);
    if length(id) >= cellthresh
        validid = [validid;i];
    end
end

spike_matrix = spike_matrix(validid,:);
spike_matrix_label = spike_matrix_label(validid);
%%
% remove the file folder that has functions with the same name as those in the UMAP toolbox
rmpath(genpath('/Users/wenbotang/Src_Matlab'))
% project the F's trajectories using the template from N'
if CTXHP
    spike_umap = run_umap(spike_matrix,'min_dist',0.6,'n_neighbors',50,'metric','cosine','n_components',3,'template_file','all_PFCumap_trials_ep6.mat');
else
    spike_umap = run_umap(spike_matrix,'min_dist',0.6,'n_neighbors',50,'metric','cosine','n_components',3,'template_file','all_CA1umap_trials_ep6.mat');
end

spike_umap_ep4 = spike_umap;
spike_umap_ep4_label = spike_matrix_label;
%%
% calculate the within-cluster distance
umap_dist_ep4 = zeros(length(spike_umap_avg_label),1);
for i = 1:length(spike_umap_avg_label)
    current_label = spike_umap_avg_label(i);
    label_ids = find(spike_umap_ep4_label == current_label);
    spike_umap_ep4_tmp = spike_umap_ep4(label_ids,:);
    dist_tmp = sum((spike_umap_ep4_tmp'-spike_umap_avg(i,:)').^2).^0.5; % Euclidean distance
    umap_dist_ep4(i) = nanmean(dist_tmp);
end
%%
%-------plot results-------%
% load shuffled result
if CTXHP
    load('PFCumap_dist_shuf-N.mat')
else
    load('CA1umap_dist_shuf-N.mat')
end
figure('Position',[300,300,200,300]),
groupnames = [ones(size(umap_dist_ep4));2 * ones(size(umap_dist_shuf))];
boxplot([umap_dist_ep4;umap_dist_shuf],groupnames,'Labels',{'F-N2','Shuffle-N2'})
ylim([0,20])
ylabel('Distance (a.u.)')
if CTXHP
   title('PFC')
else
   title('CA1')
end
h=findobj(gca,'tag','Outliers');delete(h) % not plotting outliers
h = findobj(gca,'Tag','Box');

% adjust whiskers to 10%-90%
up_adj = findobj(gca,'type', 'line', 'tag', 'Upper Whisker'); 
up_adjYData = {up_adj.YData};
up_adjYData{1}(1,2)= quantile(umap_dist_shuf,0.9);% Setting the upper whiskers
up_adjYData{2}(1,2)= quantile(umap_dist_ep4,0.9);% Setting the upper whiskers
[up_adj.YData] = deal(up_adjYData{:});  

up_adj = findobj(gca,'type', 'line', 'tag', 'Upper Adjacent Value'); 
up_adjYData = {up_adj.YData};
up_adjYData{1}(:)= quantile(umap_dist_shuf,0.9);% Setting the upper whiskers
up_adjYData{2}(:)= quantile(umap_dist_ep4,0.9);% Setting the upper whiskers
[up_adj.YData] = deal(up_adjYData{:});  

low_adj = findobj(gca,'type', 'line', 'tag', 'Lower Whisker'); 
low_adjYData = {low_adj.YData};
low_adjYData{1}(1,1)= quantile(umap_dist_shuf,0.1);% Setting the lower whiskers
low_adjYData{2}(1,1)= quantile(umap_dist_ep4,0.1);% Setting the lower whiskers
[low_adj.YData] = deal(low_adjYData{:});  

low_adj = findobj(gca,'type', 'line', 'tag', 'Lower Adjacent Value'); 
low_adjYData = {low_adj.YData};
low_adjYData{1}(:)= quantile(umap_dist_shuf,0.1);% Setting the lower whiskers
low_adjYData{2}(:)= quantile(umap_dist_ep4,0.1);% Setting the lower whiskers
[low_adj.YData] = deal(low_adjYData{:});  
