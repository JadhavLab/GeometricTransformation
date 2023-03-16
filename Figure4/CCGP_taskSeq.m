%---------------------------------------------------------------%
%  This is the main script for CCGP of task sequences           %
%  using linear SVMs; Libsvm MATLAB toolbox (V3.22) is needed   %
%  (https://www.csie.ntu.edu.tw/~cjlin/libsvm/).                %
%  -- Wenbo Tang (Jan 07, 2023)                                 %
%---------------------------------------------------------------%
clc
clear all
close all
addpath(genpath('/Users/wenbotang/MATLAB/libsvm-3.22')); % add the libsvm toolbox in path
%%
pos_interp = (0:0.005:1)';% position bins, normalized
CTXHP = 0;% CA1 = 0 , PFC = 1;
figopt = 1; % plot  CCGP result?
seg_id = [0.3,0.8]; % phase id separate early vs. late phases; 30% and 80%
% separate clusters into different sets
testing_sets = [1,3,0.3;4,2,0.3;1,3,0.8;4,2,0.8]; %[traj_id1,traj_id2, seg_id]
% load clusters (dichotomies)
if CTXHP
    load('PFC_clusters_data');
else
    load('CA1_clusters_data');
end
%%
% CCGP using linear SVMs
training_ids = combnk(1:4,3); % leave one data set out
% combination loop
for current_set = 1:length(testing_sets(:,1))
    % gather data
    Features = [];
    trajlabel = [];
    Features_test = [];
    trajlabel_test = [];
    % testing sets were not used in training the classifer
    if current_set == 1 || current_set == 3
       training_sets = testing_sets([2,4],:);
    else
       training_sets = testing_sets([1,3],:);
    end
    
    % label training sets, task sequence 1 vs. 2
    for i = 1:2
        current_tr = training_sets(i,1);
        if training_sets(i,3) < seg_id(2)
            tempid = find(spike_matrix_label_seg < (seg_id(1)+current_tr-1) & spike_matrix_label_seg >= (current_tr-1));
        else
            tempid = find(spike_matrix_label_seg < (seg_id(2)+current_tr-1) & spike_matrix_label_seg >= (seg_id(1)+current_tr-1));
        end
        current_Features = spike_umap_seg(tempid,:);
        current_trajlabel = 1*ones(length(tempid),1);
        
        Features = [Features;current_Features];
        trajlabel = [trajlabel;current_trajlabel];
        
        current_tr = training_sets(i,2);
        if training_sets(i,3) < seg_id(2)
            tempid = find(spike_matrix_label_seg < (seg_id(2)+current_tr-1) & spike_matrix_label_seg >= (seg_id(1)+current_tr-1));
        else
            tempid = find(spike_matrix_label_seg < (seg_id(1)+current_tr-1) & spike_matrix_label_seg >= (current_tr-1));
        end
        current_Features = spike_umap_seg(tempid,:);
        current_trajlabel = 2*ones(length(tempid),1);
        
        Features = [Features;current_Features];
        trajlabel = [trajlabel;current_trajlabel];
        
    end
    % training sets
    dat_all.Features = Features';
    dat_all.trajlabel = trajlabel';
    
    % label testing sets, task sequence 1 vs. 2
    for i = 1
        current_tr = testing_sets(current_set,1);
        if testing_sets(current_set,3) < seg_id(2)
            tempid = find(spike_matrix_label_seg < (seg_id(1)+current_tr-1) & spike_matrix_label_seg >= (current_tr-1));
        else
            tempid = find(spike_matrix_label_seg < (seg_id(2)+current_tr-1) & spike_matrix_label_seg >= (seg_id(1)+current_tr-1));
        end
        current_Features = spike_umap_seg(tempid,:);
        current_trajlabel = 1*ones(length(tempid),1);
      
        Features_test = [Features_test;current_Features];
        trajlabel_test = [trajlabel_test;current_trajlabel];
        
        current_tr = testing_sets(current_set,2);
        if testing_sets(current_set,3) < seg_id(2)
            tempid = find(spike_matrix_label_seg < (seg_id(2)+current_tr-1) & spike_matrix_label_seg >= (seg_id(1)+current_tr-1));
        else
            tempid = find(spike_matrix_label_seg < (seg_id(1)+current_tr-1) & spike_matrix_label_seg >= (current_tr-1));
        end
        current_Features = spike_umap_seg(tempid,:);
        current_trajlabel = 2*ones(length(tempid),1);
      
        Features_test = [Features_test;current_Features];
        trajlabel_test = [trajlabel_test;current_trajlabel];
        
    end
    dat_all.Features_test = Features_test';
    dat_all.trajlabel_test = trajlabel_test';
    
    %%
    %--- plot the labeling results, binary clusters (dichotomy),task sequence 1 vs. 2---%
    if figopt
        %plot the training sets
        figure('position',[400,400,350,300]),
        id1 = find(trajlabel == 1);
        scatter3(Features(id1,1),Features(id1,2),Features(id1,3),25,'x','MarkerEdgeColor',[153,153,153]/255,'MarkerFaceColor',[100,100,100]/255,...
            'MarkerEdgeAlpha',.5,'MarkerFaceAlpha',.2,'linewidth',2)
        hold on
        id1 = find(trajlabel == 2);
        scatter3(Features(id1,1),Features(id1,2),Features(id1,3),25,'x','MarkerEdgeColor',[255,193,194]/255,'MarkerFaceColor',[255,102,102]/255,...
            'MarkerEdgeAlpha',.5,'MarkerFaceAlpha',.2,'linewidth',2)
        %plot the testing sets
        id1 = find(trajlabel_test == 1);
        scatter3(Features_test(id1,1),Features_test(id1,2),Features_test(id1,3),70,'MarkerFaceColor',[100,100,100]/255,'MarkerEdgeColor','k',...
            'MarkerEdgeAlpha',.8,'MarkerFaceAlpha',.5)

        hold on
        id1 = find(trajlabel_test == 2);
        scatter3(Features_test(id1,1),Features_test(id1,2),Features_test(id1,3),70,'MarkerFaceColor',[255,193,194]/255,'MarkerEdgeColor',[255,102,102]/255,...
            'MarkerEdgeAlpha',.8,'MarkerFaceAlpha',.5)

    end
    %%
    %--- CCGP using linear SVMs ---%
    [~,accuracy_test,accuracy_shuffle,~,~] = decoding_CCGP_linearSVM_simple(dat_all,figopt);
    accuracy_all(current_set) = accuracy_test;
    accuracy_shuf{current_set} = accuracy_shuffle;
end
%%
% gather results of trial-shuffled data
accuracy_shuf_all = [];
for i = 1:4
    accuracy_shuf_all = [accuracy_shuf_all;accuracy_shuf{i}'];
end
%%
%-------plot real accuracy against shuffles-------%
figure('Position',[300,300,200,300]),
boxplot(accuracy_shuf_all,'Labels',{'Task seq.'})
hold on
scatter([0.95,0.975,1.025,1.05],accuracy_all,40,'k^')
ylim([0,1])
ylabel('CCGP (%)')
if CTXHP
    title('PFC')
else
    title('CA1')
end
