%---------------------------------------------------------------%
%  This is the main script for UMAP transformation of neuronal  %
%  population activity. The UMAP MATLAB toolbox (V4.1) is       %
%  needed (https://www.mathworks.com/matlabcentral/fileexchange/%
%  71902-uniform-manifold-approximation-and-projection-umap).   %
%  -- Wenbo Tang (Jan 07, 2023)                                 %
%---------------------------------------------------------------%
clc
clear all
close all
%%
animalprefix_list = {'AM2','JS17','JS21','ZT2','JS34'};% all animals
famnov_dir = '/Volumes/NovelFamiliar'; % set directory
day = 1;
ep = 6; %epoch, F = 4, N' = 6
CTXHP = 1;% CA1 = 0 , PFC = 1;
pos_interp = (0:0.005:1)';% position bins, normalized
cellthresh = 5; % minimal number of active cells per bin needed
usetemplate = 1; % a template of the UMAP transformation using N' has been save. Use it to get the exact image in the paper (=1), or run a new UMAP transformation on the data yourself (=0), which should generate a quantitatively similar result.
%%
% get trial number for each animal (the minimal number of trials across animals 
% starting from the first trial in a given session was used as the final trial
% number).
trialnum_all = []; % reset matrix
for theAnimal = 1:length(animalprefix_list)
    animalprefix = animalprefix_list{theAnimal};
    dir = sprintf('%s/%s_direct/',famnov_dir,animalprefix); % set directory
    %%
    % load trial information
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
trialnum = min(trialnum_all); %the minimal number of trials across animals was used

%%
% get firing rates for each trial
spike_matrix = []; % reset matrix for neural responses
spike_matrix_label = []; % reset matrix for trial labels

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
    load(sprintf('%s%slinfields_singletrial0%d.mat',dir,animalprefix,day)); % get single-trial firing rates
    %%
    spike_matrix_animal = []; % matrix for each animal
    spike_matrix_label = [];
    %-----get cell indices-----%
    if CTXHP
        [cellidx, ~] = matchidx_acrossep(dir, animalprefix, day,exclude_list); %cellidx = (tet, cell)
    else
        [~, cellidx] = matchidx_acrossep(dir, animalprefix, day,exclude_list); %cellidx = (tet, cell)
    end
    cellnum = length(cellidx(:,1));
    for tracknum = 1:4 % 4 different trajctories
        for tj = 1:trialnum(tracknum)
            temp_label = pos_interp + tracknum -1; % trial label
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
rmpath(genpath('/Users/wenbotang/Src_Matlab')) % remove the file folder that has functions with the same name as those in the UMAP toolbox
if CTXHP
    % PFC
    if usetemplate % use template?
        spike_umap = run_umap(spike_matrix,'min_dist',0.6,'n_neighbors',50,'metric','cosine','n_components',3,'template_file','all_PFCumap_trials_ep6.mat');
    else
        spike_umap = run_umap(spike_matrix,'min_dist',0.6,'n_neighbors',50,'metric','cosine','n_components',3,'save_template_file','all_PFCumap_trials_newtemplate.mat');
    end

else
    % CA1
    if usetemplate
        spike_umap = run_umap(spike_matrix,'min_dist',0.6,'n_neighbors',50,'metric','cosine','n_components',3,'template_file','all_CA1umap_trials_ep6.mat');
    else
        spike_umap = run_umap(spike_matrix,'min_dist',0.6,'n_neighbors',50,'metric','cosine','n_components',3,'save_template_file','all_PFCumap_trials_newtemplate.mat');
    end
end
%%
%----plot UMAP results----%
% use the matplotlib toolbox for colormaps
traj1 = inferno(length(pos_interp));
traj2 = parula(length(pos_interp));
traj3 = viridis(length(pos_interp));
traj4 = plasma(length(pos_interp));
traj_colormap = [traj1;traj2(end:-1:1,:);traj3;traj4(end:-1:1,:)];% align colors from the trajectory start to end

figure,
scatter3(spike_umap(:,1),spike_umap(:,2),spike_umap(:,3),30*ones(size(spike_matrix_label)),spike_matrix_label,'filled')
colormap(traj_colormap)
if CTXHP
    ylim([-10,10])
    xlim([-10,10])
    zlim([-15,0])
    view([156,36])
    title('PFC')
else
    zlim([-15,5])
    ylim([-25,10])
    xlim([-15,15])
    view([78,60])
    title('CA1')
end
xlabel('UMAP Dim 1')
ylabel('UMAP Dim 2')
zlabel('UMAP Dim 3')
