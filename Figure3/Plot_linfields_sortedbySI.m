%---------------------------------------------------------------%
%  This is the script for plotting all linearized rate maps     %
%  sorted by trajectory selectivity index                       %
%  -- Wenbo Tang (Jan 06, 2023)                                 %
%---------------------------------------------------------------%
clc
close all
clear all;
%%
animalprefix_list = {'AM2','JS17','JS21','ZT2','JS34'}; % all animals
famnov_dir = '/Volumes/NovelFamiliar/'; % set directory
day = 1;
eps = [2,4,6];
CTXHP = 0; % CTXHP = 1, PFC; CTXHP = 0, CA1
pos_interp = (0:0.01:1)';% position bins, normalized, 2-cm spatial bin
%%
% gather linearized rate maps
for e = 1:length(eps)
    ep = eps(e);
    rm1 = []; % ratemap matrix
    peakrate1 = []; % peak rate matrix
    peakloc1 = []; % peak location matrix
    SImean_in1 = []; % SI matrix

    for theAnimal = 1:length(animalprefix_list)
        animalprefix = animalprefix_list{theAnimal};   
        disp(['Animal: ',animalprefix,' Epoch:',num2str(ep)]) % show progress
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

        % get cell indices, exclude interneurons
        if CTXHP
            [cellind, ~] = poolcells_acrossep(dir, animalprefix, day, eps, exclude_list); %cellind = (tet, cell)
        else
            [~, cellind] = poolcells_acrossep(dir, animalprefix, day, eps, exclude_list); %cellind = (tet, cell)
        end
        cellnum = length(cellind(:,1));
        %%
        %-----create the ratemaps-----%
        load(sprintf('%s%slinfields0%d.mat',dir,animalprefix,day)); % get linearized place fields
        nn = 0;% reset counts
       
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
                      % arrange the positions from the start to the end of
                      % the trajectory
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
                  rm1_peakloc = pos_interp(rm1_peakloc)+track-1;
                  linfield_hp(isnan(linfield_hp)) = 0;

                  rm_cell = [rm_cell,linfield_hp'];
                  peakrate_cell = [peakrate_cell,rm1_peak];
                  peakloc_cell = [peakloc_cell,rm1_peakloc];
              end
              meanrate1 = nanmean(rm_cell(102:202)); % INL
              meanrate2 = nanmean(rm_cell(304:404)); % INR
              SImean_in_cell = (meanrate1-meanrate2)./(meanrate1+meanrate2); % SI for inbounds
             
              SImean_in1 = [SImean_in1;SImean_in_cell];
              rm1 = [rm1;rm_cell];
              peakrate1 = [peakrate1;peakrate_cell];
              peakloc1 = [peakloc1;peakloc_cell];
        end
     end
     rm{ep}.ratemap = rm1;
     rm{ep}.peakrate = peakrate1;
     rm{ep}.peakloc = peakloc1;
     rm{ep}.SImean_in = SImean_in1;

end
%%
% gather data
SI_all = [];
ratemap_all = [];
for ep = eps
    validid = find(max(rm{ep}.peakrate(:,[2,4])') > 3 & max(rm{ep}.peakrate(:,[1,3])') > 3 & abs(rm{ep}.SImean_in') > 0.1);
    SImean_in_tmp = rm{ep}.SImean_in(validid);
    ratemap_tmp = rm{ep}.ratemap(validid,:)./max(rm{ep}.peakrate(validid,:)')'; % normalized by the peak firing rate
    ratemap_tmp(isnan(ratemap_tmp)) = 0;
    SI_all = [SI_all;SImean_in_tmp];
    ratemap_all = [ratemap_all;ratemap_tmp];
end
%%
% sort cell by SI(IN)
[SI_sorted,sortedid] = sort(SI_all);
ratemap_all = ratemap_all(sortedid,:);
ratemap_all = ratemap_all(:,[304:404,1:303]);
%%
% plot results
figure('position',[300,300,150,400]),
m =100;
cm_viridis= inferno(m); % use the matplotlib toolbox
imagesc(ratemap_all)
hold on
% plot the line for SI = 0
if CTXHP
    plot([1,404],[120.5,120.5],'w--')%PFC
    title('PFC')
else
    plot([1,404],[121.5,121.5],'w--')%CA1
    title('CA1')
end
for track = 1:3
    plot([1.5+100*track,1.5+100*track],[1,length(ratemap_all(:,1))],'w--')%trajectory boundary
end
colormap(cm_viridis)
caxis([0.2,1])
ylabel('Cell sample # (sorted by SI(IN))')
xlabel('Position bins')
