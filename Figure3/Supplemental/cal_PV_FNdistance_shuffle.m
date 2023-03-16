%---------------------------------------------------------------%
%  This is the main script for calculating the distance of      %
%  neural states between N' and shuffled neural activity in the %
%  original state space.   -- Wenbo Tang (Jan 07, 2023)         % 
%---------------------------------------------------------------%
clc
clear all
close all
%%
animalprefix_list = {'AM2','JS17','JS21','ZT2','JS34'}; % all animals
famnov_dir = '/Volumes/NovelFamiliar'; % set directory
day = 1;
ep = 4;
CTXHP = 1;% CA1 = 0 , PFC = 1;
pos_interp = (0:0.005:1)';% position bins, normalized
cellthresh = 5;% minimal number of active cells per bin needed
%%
% calculate the averaged traj
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
    load(sprintf('%s%slinfields_singletrial0%d.mat',dir,animalprefix,day)); % get single-trial firing rates
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
spike_matrix_avg = spike_matrix;
spike_matrix_avg_label = spike_matrix_label;
%%
% calculate shuffle matrix
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
% gather shuffled data
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
    load(sprintf('%s%slinfields_singletrial0%d.mat',dir,animalprefix,day)); % get single-trial firing rates
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
                %----shuffle cell IDs----%
                cellid_shuf = randperm(cellnum);
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
                    temp(:,cellid_shuf(i)) = linfield_hp';
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

spike_matrix_shuf = spike_matrix;
spike_matrix_shuf_label = spike_matrix_label;
%%
% calculate the within-cluster distance
PV_dist_shuf = zeros(length(spike_matrix_avg_label),1);
for i = 1:length(spike_matrix_avg_label)
    current_label = spike_matrix_avg_label(i);
    label_ids = find(spike_matrix_shuf_label == current_label);
    spike_matrix_shuf_tmp = spike_matrix_shuf(label_ids,:);
    dist_tmp = sum( (spike_matrix_shuf_tmp'-spike_matrix_avg(i,:)').^2).^0.5;% Euclidean distance
    PV_dist_shuf(i) = nanmedian(dist_tmp);
end
%%
% save the shuffled the result for plotting with the real distance
if CTXHP
    save('PFCPV_dist_shuf-N.mat','PV_dist_shuf')
else
    save('CA1PV_dist_shuf-N.mat','PV_dist_shuf')
end
