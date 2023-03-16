%---------------------------------------------------------------%
%  This is the script for plotting the linearized rate maps for %
%  spatially-tuned cells  -- Wenbo Tang (Jan 05, 2023)          %
%---------------------------------------------------------------%
clc
close all
clear all;
%%
% set directory
famnov_dir = '/Volumes/NovelFamiliar/';
%%
% the animal list
animalprefix_list = {'AM2','JS17','JS21','ZT2','JS34'};
day = 1;
eps = [2,4,6];
CTXHP = 0;% CTXHP = 1, PFC; CTXHP = 0, CA1
pos_interp = (0:0.01:1)';%position bins; 2-cm spatial bin
%%
% gather rate maps 
for e = 1:length(eps)
    ep = eps(e);
    
    rm1 = []; % ratemap matrix
    peakrate1 = []; % peak rate matrix
    peakloc1 = []; % peak location matrix
    
    for theAnimal = 1:length(animalprefix_list)
        animalprefix = animalprefix_list{theAnimal};   
        dir = sprintf('%s/%s_direct/',famnov_dir,animalprefix);
        % get the cell list for interneurons
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
        if CTXHP
            [cellind, ~] = poolcells_acrossep(dir, animalprefix, day, eps, exclude_list); %cellind = (tet, cell)
        else
            [~, cellind] = poolcells_acrossep(dir, animalprefix, day, eps, exclude_list); %cellind = (tet, cell)
        end
        cellnum = length(cellind(:,1));
        %%
        %-----create the ratemaps-----%
        load(sprintf('%s%slinfields0%d.mat',dir,animalprefix,day)); % get linearized place fields
        nn = 0; % reset counter
        for i = 1:cellnum
              cind = cellind(i,:);
              rm_cell = [];
              peakrate_cell = [];
              peakloc_cell = [];
              if (length(linfields{day}{ep})>= cind(1))
                    if  (length(linfields{day}{ep}{cind(1)})>= cind(2))
                        linfield1 = linfields{day}{ep}{cind(1)}{cind(2)};
                    else 
                        linfield1 =[];
                    end
              else
                    linfield1=[];
              end
              for track = 1:4 % 4 different trajectories
                  if ~isempty(linfield1)
                      temp1 = linfield1{track};
                      pos1 = temp1(:,1);
                      linfield_hp = temp1(:,5);
                      if mod(track,2)
                          linfield_hp = interp1(pos1./max(pos1),linfield_hp,pos_interp,'nearest');
                      else
                          linfield_hp = interp1(pos1./max(pos1),linfield_hp(end:-1:1),pos_interp,'nearest');
                      end
                      
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
                      
                  else
                       linfield_hp = zeros(size(pos_interp));
                  end
                  [rm1_peak,rm1_peakloc] = max(linfield_hp);
                  rm1_peakloc = pos_interp(rm1_peakloc);
                  linfield_hp(isnan(linfield_hp)) = 0;
                  idx = find(linfield_hp > 0.25*rm1_peak);
                  rm1 = [rm1;linfield_hp'];
                  peakrate1 = [peakrate1;rm1_peak];
                  peakloc1 = [peakloc1;rm1_peakloc];
              end
         end
     end
     rm{ep}.ratemap = rm1;
     rm{ep}.peakrate = peakrate1;
     rm{ep}.peakloc = peakloc1;
end
%%
% plot cells with a stable peak in the novel environment
validid = find(abs(rm{6}.peakloc - rm{2}.peakloc) < 0.2);

% sort cells by the peak firing rate in N'
[peakloc_sorted,sortedid] = sort(rm{6}.peakloc(validid));

% normalize by the peak firing rate 
ratemap3 = rm{6}.ratemap(validid(sortedid),:)./rm{6}.peakrate(validid(sortedid));
ratemap2 = rm{4}.ratemap(validid(sortedid),:)./rm{4}.peakrate(validid(sortedid));
ratemap1 = rm{2}.ratemap(validid(sortedid),:)./rm{2}.peakrate(validid(sortedid));
% remove NaNs
ratemap1(isnan(ratemap1)) = 0;
ratemap2(isnan(ratemap2)) = 0;
ratemap3(isnan(ratemap3)) = 0;

% sort the cell without a valid rate map in N'
zerocellid = find(peakloc_sorted == 0);
allrate = sum(ratemap3(zerocellid,:),2);
[~,indx2]  = sort(allrate);
ratemap3(zerocellid,:) = ratemap3(zerocellid(indx2),:);
ratemap2(zerocellid,:) = ratemap2(zerocellid(indx2),:);
ratemap1(zerocellid,:) = ratemap1(zerocellid(indx2),:);
%%
% plot the result
figure('position',[1150 500 800 300])

m=100;
cm_viridis= inferno(m); % use the matplotlib toolbox

subplot(1,4,1)
imagesc(ratemap1)
colormap(cm_viridis)
caxis([0.2 1])
subplot(1,4,2)
imagesc(ratemap2)
colormap(cm_viridis)
caxis([0.2 1])
subplot(1,4,3)
imagesc(ratemap3)
caxis([0.2 1])
colormap(cm_viridis)

subplot(1,4,4)
% sort F rate maps
[peakloc_sorted,sortedid] = sort(rm{4}.peakloc(validid));
ratemap2 = rm{4}.ratemap(validid(sortedid),:)./rm{4}.peakrate(validid(sortedid));
ratemap2(isnan(ratemap2)) = 0;
zerocellid = find(peakloc_sorted == 0);
allrate = sum(ratemap2(zerocellid,:),2);
[~,indx2]  = sort(allrate);
ratemap2(zerocellid,:) = ratemap2(zerocellid(indx2),:);
imagesc(ratemap2)
colormap(cm_viridis)
caxis([0.2 1])
