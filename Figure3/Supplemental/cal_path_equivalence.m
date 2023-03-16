%---------------------------------------------------------------%
%  This is the script for calculating path-equivalent           %
%  coefficients for all cells    -- Wenbo Tang (Jan 07, 2023)   %
%---------------------------------------------------------------%
clc
clear all;
close all;
%%
animalprefix_list = {'AM2','JS17','JS21','ZT2','JS34'}; % all animals
famnov_dir = '/Volumes/NovelFamiliar'; % set directory
day = 1;
eps = [2,4,6]; % epochs = N, F, N'
CTXHP = 0; %PFC = 1, CA1 = 0
%%
% gather path-equivalent coefficients
path_corr = []; % reset matrix
for theAnimal = 1:length(animalprefix_list)
    animalprefix = animalprefix_list{theAnimal};
    dir = sprintf('%s/%s_direct/',famnov_dir,animalprefix);
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
    %%
    % get cell indices, exclude interneurons
    if CTXHP
        [cellind, ~] = matchidx_acrossep(dir, animalprefix, day, exclude_list); %(tet, cell)
    else
        [~, cellind] = matchidx_acrossep(dir, animalprefix, day, exclude_list); %(tet, cell)
    end
    cellnum = length(cellind(:,1));
    
    %%
    load(sprintf('%s%slinfields0%d.mat',dir,animalprefix,day)); % get linearized place fields
    %%
    temp = zeros(cellnum,length(eps)); %reset matrix
    %-----gather rate maps-----%
    for i = 1:cellnum
          cind = cellind(i,:);
          for ep = eps
              if (length(linfields{day}{ep})>= cind(1))
                    if  (length(linfields{day}{ep}{cind(1)})>= cind(2))
                        linfield1 = linfields{day}{ep}{cind(1)}{cind(2)};
                        linfield_all = [];
                        for track = 1:4
                            temp1 = linfield1{track};
                            pos1 = temp1(:,1);
                            pos1 = pos1./max(pos1);
                            if track == 1
                                pos_match = pos1;
                            end
                            occnormrate1 = temp1(:,5);
                            a = find(~isnan(occnormrate1));

                            if ~mod(track,2) % align trajectory from the start to end
                                pos1 = 1 - pos1;
                            end
                            occnormrate_tmp = nan(size(occnormrate1));
                            % make trajectories the same length
                            occnormrate_tmp = interp1(pos1(a),occnormrate1(a),pos_match);
                            linfield_all = [linfield_all, occnormrate_tmp]; 
                        end
                        % calculate r between a pair of trajectories
                        pairind = combnk(1:4,2);
                        nonnanbins = find(~isnan(sum(linfield_all,2)));
                        for p = 1:length(pairind(:,1))
                            rval(p) = corr(linfield_all(nonnanbins,pairind(p,1)),linfield_all(nonnanbins,pairind(p,2)));
                        end
                        % take the median across all pairs
                        temp(i,ep/2) = nanmedian(rval);
                    else
                        temp(i,ep/2) = nan; % no valid rate map found for the cell
                    end
              else
                  temp(i,ep/2) = nan; % no valid rate map found for the cell
              end
          end
    end
    path_corr = [path_corr;temp]; % gather data
end
%%
% remove the file folder that has functions with the same name as the
% plotting functions
rmpath(genpath('/Users/wenbotang/Src_Matlab'))
%-------plot results-------%
figure('Position',[300,300,400,300]),
groupnames = [ones(size(path_corr(:,1)));2 * ones(size(path_corr(:,1)));3 * ones(size(path_corr(:,1)))];
boxplot([path_corr(:,1);path_corr(:,2);path_corr(:,3)],groupnames,'Labels',{'N','F','N2'})
ylim([-0.5,1])
ylabel('Path-equivalent coefficient (R)')
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
for i=1:3
    up_adjYData{i}(1,2)= quantile(path_corr(:,end-i+1),0.9);% Setting the upper whiskers
end
[up_adj.YData] = deal(up_adjYData{:});  

up_adj = findobj(gca,'type', 'line', 'tag', 'Upper Adjacent Value'); 
up_adjYData = {up_adj.YData};
for i=1:3
    up_adjYData{i}(:)= quantile(path_corr(:,end-i+1),0.9);% Setting the upper whiskers
end
[up_adj.YData] = deal(up_adjYData{:});  

low_adj = findobj(gca,'type', 'line', 'tag', 'Lower Whisker'); 
low_adjYData = {low_adj.YData};
for i=1:3
    low_adjYData{i}(1,1)= quantile(path_corr(:,end-i+1),0.1);% Setting the lower whiskers
end
[low_adj.YData] = deal(low_adjYData{:});  

low_adj = findobj(gca,'type', 'line', 'tag', 'Lower Adjacent Value'); 
low_adjYData = {low_adj.YData};
for i=1:3
    low_adjYData{i}(:)= quantile(path_corr(:,end-i+1),0.1);% Setting the lower whiskers
end
[low_adj.YData] = deal(low_adjYData{:});  
