%---------------------------------------------------------------%
%  This is the main script for decoding all dichotomies using   %
%  linear SVMs; Libsvm MATLAB toolbox (V3.22) is needed         %
%  (https://www.csie.ntu.edu.tw/~cjlin/libsvm/).                %
%  -- Wenbo Tang (Jan 08, 2023)                                 %
%---------------------------------------------------------------%
clc
clear all
close all
addpath(genpath('/Users/wenbotang/MATLAB/libsvm-3.22')); % add the libsvm toolbox in path
%%
pos_interp = (0:0.005:1)';% position bins, normalized
CTXHP = 0;% CA1 = 0 , PFC = 1;
seg_id = [0.3,0.8]; % phase id separate early vs. late phases; 30% and 80%

% feature groups, [traj_type, seg_id, seq_id]
f_g = [1,0.3,0;2,0.8,1;3,0.3,1;4,0.8,0;1,0.8,0;2,0.3,1;3,0.8,1;4,0.3,0];

% define all 35 dichotomies
C_g = combnk(1:8,4); %randomly take 4 clusters out of 8
C_g1 = C_g(1:length(C_g)/2,:);
C_g2 = C_g(end:-1:1+length(C_g)/2,:);

% load clusters (dichotomies)
if CTXHP
    load('PFC_clusters_data');
else
    load('CA1_clusters_data');
end
%%
% dichotomy combination loop
for g = 1:length(C_g1(:,1))
    Features = []; % reset matrix
    trajlabel = [];
  
    % gather data
    training_group1 = f_g(C_g1(g,:),:);
    training_group2 = f_g(C_g2(g,:),:);
    % group 1, 4 clusters
    for i = 1:4
        current_tr = training_group1(i,1);
        current_seg = training_group1(i,2);
        if current_seg < seg_id(2)
            tempid = find(spike_matrix_label_seg < (seg_id(1)+current_tr-1) & spike_matrix_label_seg >= (current_tr-1));
        else
            tempid = find(spike_matrix_label_seg < (seg_id(2)+current_tr-1) & spike_matrix_label_seg >= (seg_id(1)+current_tr-1));
        end
        current_Features = spike_umap_seg(tempid,:);
        current_trajlabel = 1*ones(length(tempid),1);
        
        Features = [Features;current_Features];
        trajlabel = [trajlabel;current_trajlabel];
    end
    % group 2, 4 clusters
    for i = 1:4
        current_tr = training_group2(i,1);
        current_seg = training_group2(i,2);
        if current_seg < seg_id(2)
            tempid = find(spike_matrix_label_seg < (seg_id(1)+current_tr-1) & spike_matrix_label_seg >= (current_tr-1));
        else
            tempid = find(spike_matrix_label_seg < (seg_id(2)+current_tr-1) & spike_matrix_label_seg >= (seg_id(1)+current_tr-1));
        end
        current_Features = spike_umap_seg(tempid,:);
        current_trajlabel = 2*ones(length(tempid),1);
        
        Features = [Features;current_Features];
        trajlabel = [trajlabel;current_trajlabel];
    end
    
    dat_all.Features = Features';
    dat_all.trajlabel = trajlabel';
    %%
    %--- predictions using linear SVMs ---%
    [accuracy,~,accuracy_shuffle] = decoding_dichotomy_linearSVM(dat_all);
    accuracy_all(g) = accuracy; % accuracy for real data
    accuracy_shuf{g} = accuracy_shuffle;% accuracy for trial-label shuffled data
end
%%
% gather results of trial-shuffled data
accuracy_shuf_all = [];
for i = 1:35
    accuracy_shuf_all = [accuracy_shuf_all;accuracy_shuf{i}'];
end
%%
%-------plot real accuracy against shuffles-------%
xi = 0.95:(1.05-0.95)/34:1.05;% x values for stacking scatter plot
figure('Position',[300,300,200,300]),
if CTXHP
    boxplot(accuracy_shuf_all,'Labels',{'PFC'})
else
    boxplot(accuracy_shuf_all,'Labels',{'CA1'})
end
hold on
scatter(xi,accuracy_all,40,'ko')
ylim([0,100])
ylabel('Decoding accuracy (%)')
