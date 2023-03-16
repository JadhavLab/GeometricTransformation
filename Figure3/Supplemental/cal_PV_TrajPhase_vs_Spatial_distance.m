%---------------------------------------------------------------%
%  This is the main script for calculating the distance of      %
%  neural states of the same spatial location vs. the same      %
%  trajectory phase in the original state space.                %
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
CTXHP = 0;% CA1 = 0 , PFC = 1;
pos_interp = (0:0.005:1)';% position bins, normalized
cellthresh = 5;% minimal number of active cells per bin needed
%%
trialnum_all = [];
for theAnimal = 1:length(animalprefix_list)
    animalprefix = animalprefix_list{theAnimal};
    dir = sprintf('%s/%s_direct/',famnov_dir,animalprefix); % set directory
    %%
    trajinfo = loaddatastruct(dir, animalprefix, 'trajinfo', day); % get trajinfo

    %%
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
    for tracknum = 1:4 % 4 trajectory types
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
% create an unique vector of trajectory labels
traj_label = [];
for tracknum = 1:4
    temp_label = pos_interp + tracknum -1;
    traj_label = [traj_label,temp_label];
end
%%
% remove the file folder that has functions with the same name as the
% default pdist function
rmpath(genpath('/Users/wenbotang/Src_Matlab'))
% calculate distance of neural state aligned by spatial location
PV_dist_spatial = [];
for i = 1:length(traj_label(:,1))
    for track = [1,3]
        current_label = traj_label(i,track);
        label_ids = find(spike_matrix_label == current_label);
        spike_currenttraj= spike_matrix(label_ids,:);
        if track == 1
           alt_label = traj_label(i,2);%spatial location aligned; INL-OUTL
        else
           alt_label = traj_label(i,4);%spatial location aligned; INR-OUTR
        end
        label_ids_alt = find(spike_matrix_label == alt_label);
        spike_alttraj= spike_matrix(label_ids_alt,:);
        dist_tmp = pdist2(spike_alttraj,spike_currenttraj,'euclidean'); % Euclidean distance
        PV_dist_spatial = [PV_dist_spatial;dist_tmp(:)];
    end
end
%%
% calculate distance of neural state aligned by trajectory phase
PV_dist_phase = [];
for i = 1:length(traj_label(:,1))
    for track = [1,3]
        current_label = traj_label(i,track);
        label_ids = find(spike_matrix_label == current_label);
        spike_currenttraj= spike_matrix(label_ids,:);
        if track == 1
           alt_label = traj_label(end-i+1,2);%phase aligned; INL-OUTL
        else
           alt_label = traj_label(end-i+1,4);%phase aligned; INR-OUTR
        end
        label_ids_alt = find(spike_matrix_label == alt_label);
        spike_alttraj= spike_matrix(label_ids_alt,:);
        dist_tmp = pdist2(spike_alttraj,spike_currenttraj,'euclidean');% Euclidean distance
        PV_dist_phase = [PV_dist_phase;dist_tmp(:)];
    end
end
%%
%-------plot results-------%
figure('Position',[300,300,200,300]),
groupnames = [ones(size(PV_dist_spatial));2 * ones(size(PV_dist_phase))];
boxplot([PV_dist_spatial;PV_dist_phase],groupnames,'Labels',{'Spatial','Traj-phase'})
ylim([20,180])
ylabel('PV distance (a.u.)')
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
up_adjYData{1}(1,2)= quantile(PV_dist_phase,0.9);% Setting the upper whiskers
up_adjYData{2}(1,2)= quantile(PV_dist_spatial,0.9);% Setting the upper whiskers
[up_adj.YData] = deal(up_adjYData{:});  

up_adj = findobj(gca,'type', 'line', 'tag', 'Upper Adjacent Value'); 
up_adjYData = {up_adj.YData};
up_adjYData{1}(:)= quantile(PV_dist_phase,0.9);% Setting the upper whiskers
up_adjYData{2}(:)= quantile(PV_dist_spatial,0.9);% Setting the upper whiskers
[up_adj.YData] = deal(up_adjYData{:});  

low_adj = findobj(gca,'type', 'line', 'tag', 'Lower Whisker'); 
low_adjYData = {low_adj.YData};
low_adjYData{1}(1,1)= quantile(PV_dist_phase,0.1);% Setting the lower whiskers
low_adjYData{2}(1,1)= quantile(PV_dist_spatial,0.1);% Setting the lower whiskers
[low_adj.YData] = deal(low_adjYData{:});  

low_adj = findobj(gca,'type', 'line', 'tag', 'Lower Adjacent Value'); 
low_adjYData = {low_adj.YData};
low_adjYData{1}(:)= quantile(PV_dist_phase,0.1);% Setting the lower whiskers
low_adjYData{2}(:)= quantile(PV_dist_spatial,0.1);% Setting the lower whiskers
[low_adj.YData] = deal(low_adjYData{:});  
