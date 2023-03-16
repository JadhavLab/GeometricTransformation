%---------------------------------------------------------------%
%  This is the main script for calculating INSeq vs. OUTSeq     %
%  trajectory distance in the original state space.             %
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
    for tracknum = 1:4 % 4 trajectory types
            for tj = 1:trialnum(tracknum)
                temp_label = [(pos_interp + tracknum -1),tracknum * ones(length(pos_interp),1), tj*ones(length(pos_interp),1)]; %[position, trajectory type, trial number];
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
spike_matrix_label = spike_matrix_label(validid,:);
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
% calculate INSeq distance
PV_dist_inseq = [];% reset matrix
for track = [1,3]
    if track == 1
        alt_track = 4; %INSeq pair 1, INR-OUTL
    else
        alt_track = 2; %INSeq pair 2, INL-OUTR
    end
    trial_num1 = trialnum(track);
    trial_num2 = trialnum(alt_track);
    %-----------get pairs of trials----------%
    C_1 = repmat(1:trial_num1,trial_num2,1);
    C_1 = reshape(C_1,1,[]);
    C_2 = repmat(1:trial_num2,1,trial_num1);
    pairind = [C_1' C_2'];
    
    %-----------calculate distance for each trial pair----------%
    for pair = 1:length(pairind(:,1))
        tr1 = pairind(pair,1);
        tr2 = pairind(pair,2);
        dist_tmp = [];
        for i = 1:length(traj_label(:,1))
            current_label1 = traj_label(i,track);
            label_ids = find(spike_matrix_label(:,1) == current_label1 & spike_matrix_label(:,2) == track & spike_matrix_label(:,3) == tr1);
            
            if ~CTXHP
                current_label2 = traj_label(i,alt_track);% if CA1, aligned by spatial location
            else
                current_label2 = traj_label(end-i+1,alt_track); % if PFC, aligned by trajetory phase (distance from the start)
            end

            label_ids_alt = find(spike_matrix_label(:,1) == current_label2 & spike_matrix_label(:,2) == alt_track & spike_matrix_label(:,3) == tr2);
            
            if ~isempty(label_ids) && ~isempty(label_ids_alt)
                spike_currenttraj= spike_matrix(label_ids,:);
                spike_alttraj= spike_matrix(label_ids_alt,:);
                dist_tmp = [dist_tmp;sum( (spike_currenttraj'-spike_alttraj').^2).^0.5];% Euclidean distance
            end
        end
        PV_dist_inseq = [PV_dist_inseq;nanmean(dist_tmp(:))];
    end
end

%%
% calculate OUTSeq distance
PV_dist_outseq = [];% reset matrix
for track = [1,3]
    if track == 1
        alt_track = 2;%OUTSeq pair 1, INL-OUTL
    else
        alt_track = 4;%OUTSeq pair R, INR-OUTR
    end
    trial_num1 = trialnum(track);
    trial_num2 = trialnum(alt_track);
    %-----------get pairs of trials----------%
    C_1 = repmat(1:trial_num1,trial_num2,1);
    C_1 = reshape(C_1,1,[]);
    C_2 = repmat(1:trial_num2,1,trial_num1);
    pairind = [C_1' C_2'];
    
    %-----------calculate distance for each trial pair----------%
    for pair = 1:length(pairind(:,1))
        tr1 = pairind(pair,1);
        tr2 = pairind(pair,2);
        dist_tmp = [];
        for i = 1:length(traj_label(:,1))
            current_label1 = traj_label(i,track);
            label_ids = find(spike_matrix_label(:,1) == current_label1 & spike_matrix_label(:,2) == track & spike_matrix_label(:,3) == tr1);
            
            if ~CTXHP
                current_label2 = traj_label(i,alt_track);% if CA1, aligned by spatial location
            else
                current_label2 = traj_label(end-i+1,alt_track);% if PFC, aligned by trajetory phase (distance from the start)
            end

            label_ids_alt = find(spike_matrix_label(:,1) == current_label2 & spike_matrix_label(:,2) == alt_track & spike_matrix_label(:,3) == tr2);
            
            if ~isempty(label_ids) && ~isempty(label_ids_alt)
                spike_currenttraj= spike_matrix(label_ids,:);
                spike_alttraj= spike_matrix(label_ids_alt,:);
                dist_tmp = [dist_tmp;sum( (spike_currenttraj'-spike_alttraj').^2).^0.5]; % Euclidean distance
            end
        end
        PV_dist_outseq = [PV_dist_outseq;nanmean(dist_tmp(:))];
    end
end
%%
%-------plot results-------%
figure('Position',[300,300,200,300]),
groupnames = [ones(size(PV_dist_inseq));2 * ones(size(PV_dist_outseq))];
boxplot([PV_dist_inseq;PV_dist_outseq],groupnames,'Labels',{'InSeq','OutSeq'})
ylim([40,120])
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
up_adjYData{1}(1,2)= quantile(PV_dist_outseq,0.9);% Setting the upper whiskers
up_adjYData{2}(1,2)= quantile(PV_dist_inseq,0.9);% Setting the upper whiskers
[up_adj.YData] = deal(up_adjYData{:});  

up_adj = findobj(gca,'type', 'line', 'tag', 'Upper Adjacent Value'); 
up_adjYData = {up_adj.YData};
up_adjYData{1}(:)= quantile(PV_dist_outseq,0.9);% Setting the upper whiskers
up_adjYData{2}(:)= quantile(PV_dist_inseq,0.9);% Setting the upper whiskers
[up_adj.YData] = deal(up_adjYData{:});  

low_adj = findobj(gca,'type', 'line', 'tag', 'Lower Whisker'); 
low_adjYData = {low_adj.YData};
low_adjYData{1}(1,1)= quantile(PV_dist_outseq,0.1);% Setting the lower whiskers
low_adjYData{2}(1,1)= quantile(PV_dist_inseq,0.1);% Setting the lower whiskers
[low_adj.YData] = deal(low_adjYData{:});  

low_adj = findobj(gca,'type', 'line', 'tag', 'Lower Adjacent Value'); 
low_adjYData = {low_adj.YData};
low_adjYData{1}(:)= quantile(PV_dist_outseq,0.1);% Setting the lower whiskers
low_adjYData{2}(:)= quantile(PV_dist_inseq,0.1);% Setting the lower whiskers
[low_adj.YData] = deal(low_adjYData{:});      
