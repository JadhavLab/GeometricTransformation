%---------------------------------------------------------------%
%  This is the main script for calculating population vector    %
%  (PV) similarity  -- Wenbo Tang (Jan 05, 2023)                %
%---------------------------------------------------------------%
clc
close all
clear all;
%%
% set directory
famnov_dir = '/Volumes/NovelFamiliar/';
%%
% the animal list
animal_list = {'AM2','JS17','JS21','ZT2','JS34'};
day = 1;
pos_interp = (0:0.01:1)';%position bins; 2-cm spatial bin
%%
PVcorr_all_PFC_FN = [];
PVcorr_all_CA1_FN = [];
PVcorr_all_PFC_NN = [];
PVcorr_all_CA1_NN = [];
for FN = [0,1] % FN = 1, novel-familiar enviroments; 0, novel-novel environments
    if FN == 1
        eps = [4,6];
    else
        eps = [2,6];
    end
    for CTXHP = [0,1]% CTXHP = 1, PFC; CTXHP = 0, CA1
        for animal = 1:length(animal_list)
            animalprefix = animal_list{animal};
            % set animal directory, and exclude interneurons
            dir = sprintf('%s/%s_direct/',famnov_dir,animalprefix);
            if strcmp(animalprefix,'AM2')
                exclude_list = [6,9];%AM2
            elseif strcmp(animalprefix,'JS21')
                exclude_list = [6,3;6,6;6,7;6,8;7,1;21,1;22,1;25,1];%JS21
            elseif strcmp(animalprefix,'JS17')
                exclude_list = [6,1;7,1];%JS17
            elseif strcmp(animalprefix,'ZT2')
                exclude_list = [31,1;31,5;32,2;27,1;14,10];%JS17
            elseif strcmp(animalprefix,'JS34')
                exclude_list = [];%JS34
            end

            % apply cell inclusion criterion
            if CTXHP
                [cellind, ~] = poolcells_acrossep(dir, animalprefix, day, eps, exclude_list); %cellind = (tet, cell)
            else
                [~, cellind] = poolcells_acrossep(dir, animalprefix, day, eps, exclude_list); %cellind = (tet, cell)
            end
            cellnum = length(cellind(:,1));
            %%
            %-----load the rate maps-----%
            load(sprintf('%s%slinfields0%d.mat',dir,animalprefix,day)); % get linearized place fields
            nn = 0;
            for track = 1:4 % trajectory 1-4
                for e = 1:length(eps)
                    ep = eps(e);
                    rm1 = []; % ratemap matrix
                    for i = 1:cellnum
                        cind = cellind(i,:);
                        if (length(linfields{day}{ep})>= cind(1))
                              if  (length(linfields{day}{ep}{cind(1)})>= cind(2))
                                  linfield1 = linfields{day}{ep}{cind(1)}{cind(2)};
                              else 
                                  linfield1 =[];
                              end
                        else
                              linfield1=[];
                        end

                        if ~isempty(linfield1)
                          temp1 = linfield1{track};
                          pos1 = temp1(:,1);
                          linfield_hp = temp1(:,5);

                          a = find(isnan(linfield_hp));
                          %pad nan
                          if ~isempty(a)
                             [lo,hi]= findcontiguous(a);  %find contiguous NaNs
                             for ii = 1:length(lo) 
                                 if lo(ii) > 1 & hi(ii) < length(linfield_hp)
                                     fill = linspace(linfield_hp(lo(ii)-1), ...
                                     linfield_hp(hi(ii)+1), hi(ii)-lo(ii)+1);
                                     linfield_hp(lo(ii):hi(ii)) = fill;
                                 end
                             end
                          end

                          linfield_hp = interp1(pos1./max(pos1),linfield_hp,pos_interp,'nearest');
                        else
                          linfield_hp = zeros(size(pos_interp));
                        end
                        rm1 = [rm1;linfield_hp'];
                    end
                    rm{ep}{track} = rm1;
                    pos{ep}{track} = pos1;

                end 
            end
            %%
            % gather all 4 trajectories
            rm_all1 = [rm{eps(1)}{1},rm{eps(1)}{2},rm{eps(1)}{3},rm{eps(1)}{4}];
            rm_all2 = [rm{eps(2)}{1},rm{eps(2)}{2},rm{eps(2)}{3},rm{eps(2)}{4}];

            validid = [];
            for i = 1:length(rm_all2(:,1))
                max2 = max(rm_all2(i,:));
                max1 = max(rm_all1(i,:));
                if max1 > 3 || max2 > 3 % if the cell has a peak rate larger than 3 Hz
                    validid = [validid,i];
                end
            end

            rm_all1 = rm_all1(validid,:);
            rm_all2 = rm_all2(validid,:);

            for i = 1:length(rm_all2(1,:))
                PVcorr(i) = corr(rm_all1(:,i),rm_all2(:,i));
            end
            if FN
                if CTXHP
                    PVcorr_all_PFC_FN = [PVcorr_all_PFC_FN,PVcorr'];
                else
                    PVcorr_all_CA1_FN = [PVcorr_all_CA1_FN,PVcorr'];
                end
            else
                if CTXHP
                    PVcorr_all_PFC_NN = [PVcorr_all_PFC_NN,PVcorr'];
                else
                    PVcorr_all_CA1_NN = [PVcorr_all_CA1_NN,PVcorr'];
                end
            end
            clear PVcorr
        end
    end
end
%%
% histograms
xbins = -0.4:0.05:1; % bins for histogram
hist_peranimal_PFC_FN = [];
hist_peranimal_PFC_NN = [];
hist_peranimal_CA1_FN = [];
hist_peranimal_CA1_NN = [];
for animal = 1:length(animal_list)
    current_hist = hist(PVcorr_all_PFC_FN(:,animal),xbins);
    current_hist = current_hist./length(PVcorr_all_PFC_FN(:,animal)); %normalize
    hist_peranimal_PFC_FN = [hist_peranimal_PFC_FN;current_hist];
    
    
    current_hist = hist(PVcorr_all_PFC_NN(:,animal),xbins);
    current_hist = current_hist./length(PVcorr_all_PFC_NN(:,animal)); %normalize
    hist_peranimal_PFC_NN = [hist_peranimal_PFC_NN;current_hist];
    
    current_hist = hist(PVcorr_all_CA1_FN(:,animal),xbins);
    current_hist = current_hist./length(PVcorr_all_CA1_FN(:,animal)); %normalize
    hist_peranimal_CA1_FN = [hist_peranimal_CA1_FN;current_hist];
    
    current_hist = hist(PVcorr_all_CA1_NN(:,animal),xbins);
    current_hist = current_hist./length(PVcorr_all_CA1_NN(:,animal)); %normalize
    hist_peranimal_CA1_NN = [hist_peranimal_CA1_NN;current_hist];
end

%%
% histograms pooling across animals
current_hist = hist(PVcorr_all_PFC_FN(:),xbins);
current_hist = current_hist./length(PVcorr_all_PFC_FN(:)); %normalize
hist_peranimal_PFC_FN_all = current_hist;
    
    
current_hist = hist(PVcorr_all_PFC_NN(:),xbins);
current_hist = current_hist./length(PVcorr_all_PFC_NN(:)); %normalize
hist_peranimal_PFC_NN_all = current_hist;
    
current_hist = hist(PVcorr_all_CA1_FN(:),xbins);
current_hist = current_hist./length(PVcorr_all_CA1_FN(:)); %normalize
hist_peranimal_CA1_FN_all = current_hist;
    
current_hist = hist(PVcorr_all_CA1_NN(:),xbins);
current_hist = current_hist./length(PVcorr_all_CA1_NN(:)); %normalize
hist_peranimal_CA1_NN_all = current_hist;
%%
% plot results
% plot PV similarity for each animal
figure('position',[1150 500 800 200])
newcolors_FN = [0.8 0.8 0.9; 0.6 0.6 0.9; 0.4 0.4 0.8; 0.2 0.2 0.8; 0 0 1]; % set colors for different animals
newcolors_NN = [0.9 0.8 0.8; 0.9 0.6 0.6; 0.8 0.4 0.4; 0.8 0.2 0.2; 1 0 0]; % set colors for different animals
subplot(131)
h = plot(xbins,hist_peranimal_PFC_FN,'linewidth',2);
set(h, {'color'}, num2cell(newcolors_FN,2));       
hold on
h = plot(xbins,hist_peranimal_PFC_NN,'linewidth',2);
set(h, {'color'}, num2cell(newcolors_NN,2));  
xlabel('PV similarity')
ylabel('Fraction')
title('PFC per animal')

subplot(132)
h = plot(xbins,hist_peranimal_CA1_FN,'linewidth',2);
set(h, {'color'}, num2cell(newcolors_FN,2));       
hold on
h = plot(xbins,hist_peranimal_CA1_NN,'linewidth',2);
set(h, {'color'}, num2cell(newcolors_NN,2));  
xlabel('PV similarity')
ylabel('Fraction')  
title('CA1 per animal')

%mean value for each animal
subplot(133)
newcolors = [0.83 0.14 0.14
             1.00 0.54 0.00
             0.47 0.25 0.80
             0.25 0.80 0.54
             0 1 0.8]; % set colors for different animals
PVcorr_PFC_FN_mean = nanmean(PVcorr_all_PFC_FN);
PVcorr_PFC_NN_mean = nanmean(PVcorr_all_PFC_NN);

PVcorr_CA1_FN_mean = nanmean(PVcorr_all_CA1_FN);
PVcorr_CA1_NN_mean = nanmean(PVcorr_all_CA1_NN);

PVcorr_mean = [PVcorr_CA1_FN_mean',PVcorr_CA1_NN_mean',...
    PVcorr_PFC_FN_mean',PVcorr_PFC_NN_mean'];
PVcorr_mean = PVcorr_mean';
h = plot(PVcorr_mean);
set(h, {'color'}, num2cell(newcolors,2));
hold on
h = plot(PVcorr_mean,'o');
set(h, {'color'}, num2cell(newcolors,2));

ylabel('PV similarity')
xticks(1:4)
xticklabels({'FN (CA1)','NN (CA1)','FN (PFC)','NN (PFC)'})
ylim([-0.5,1])
%%
% plot PV similarity across all animals
figure('position',[1150 500 800 300])
subplot(121)
plot(xbins,hist_peranimal_PFC_FN_all,'color',newcolors_FN(end,:),'linewidth',2);
hold on
plot(xbins,hist_peranimal_PFC_NN_all,'color',newcolors_NN(end,:),'linewidth',2);
xlabel('PV similarity')
ylabel('Fraction')
title('PFC')

subplot(122)
plot(xbins,hist_peranimal_CA1_FN_all,'color',newcolors_FN(end,:),'linewidth',2);
hold on
plot(xbins,hist_peranimal_CA1_NN_all,'color',newcolors_NN(end,:),'linewidth',2);
xlabel('PV similarity')
ylabel('Fraction')
title('CA1')
legend('F-N','N-N')
