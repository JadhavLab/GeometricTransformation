function [decodinginfo,errorpos_all] = decoding_position_CA1_novelfamiliar_trial(animalprefix,day,ep,eprun,timerange,cellthresh)
%---------------------------------------------------------------%
%  This is the function for CA1 decoding of spatial locations   %
%  using the rate maps from one epoch to decode position in a   %
%  different epoch on a trial-by-trial basis                    %
%  -- Wenbo Tang (Jan 05, 2023)                                 %
%---------------------------------------------------------------%

% INPUTS:
%
%    animalprefix = animal prefix.
%    day = experimental day.
%    ep = epoch to be decoded.
%    eprun = epoch used for rate-map templates.
%    timerange = start and end time for the testing trial.
%    cellthresh = number of active cells needed for a time bin
%
% OUTPUTS:
%
%    decodinginfo, decoding information for the time bins to be decoded during the testing trial
%    errorpos_all, decoding error fin each time bin of the testing trial
%%
famnov_dir = '/Volumes/NovelFamiliar/'; % set directory
dir = sprintf('%s/%s_direct/',famnov_dir,animalprefix);


if (day<10)
    daystring = ['0',num2str(day)];
else
    daystring = num2str(day);
end

if (ep<10)
    epochstring = ['0',num2str(ep)];
else
    epochstring = num2str(ep);
end

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
%-----use neurons active across epochs-----%
[~, hpidx] = matchidx_acrossep(dir, animalprefix, day,exclude_list); %(tet, cell)
hpnum = length(hpidx(:,1));
tBinSz = 200; %default temporal bin in ms;
%%
%-----create the ratemaps [nPosBin x nHPCells] -----%
rm = []; % ratemap matrix
pm = []; % position matrix
tm = []; % track matrix
cellidxm = []; % cell indices

load(sprintf('%s%slinfields0%d.mat',dir,animalprefix,day)); % get linearized place fields
for i = 1:hpnum
      cind = hpidx(i,:);
      if (length(linfields{day}{eprun})>= cind(1))
            if  (length(linfields{day}{eprun}{cind(1)})>= cind(2))
                linfield1 = linfields{day}{eprun}{cind(1)}{cind(2)};
            else 
                linfield1 =[];
            end
      else
            linfield1=[];
      end
            
      if ~isempty(linfield1)
           linfield_hp = [];
           lintrack_hp = [];
           pos_hp = [];
                      
           for track = 1:4
                temp1 = linfield1{track};
                pos1 = temp1(:,1);
                lintrack1 = ones(size(pos1))*track;
                occnormrate1 = temp1(:,5);
                
                linfield_hp = [linfield_hp;occnormrate1];

                pos_hp = [pos_hp;pos1];
                lintrack_hp = [lintrack_hp;lintrack1];
           end
           if (max(linfield_hp) >= 3) % peak firing rate max larger than 3 Hz
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
               
               rm = [rm;linfield_hp'];
               pm = [pm;pos_hp'];
               tm = [tm;lintrack_hp'];
               cellidxm = [cellidxm; cind];
           end
      end
end
rm = rm'; %[nPosBin x nHPCells]

pm = pm';
tm = tm';
trajinfo = mean(tm,2);% trajactory number
rm = rm+ (eps.^8); %Add a small number so there are no zeros

expecSpk =rm.*tBinSz./1000; %[nPos x nCells] Expected number of spikes per bin

hpnum = length(rm(1,:)); % update
%%
%-----load spikes, linpos and trajinfo-----%
linpos = loaddatastruct(dir, animalprefix, 'linpos', day); % get linpos

lintimes = linpos{day}{ep}.statematrix.time;
traj = linpos{day}{ep}.statematrix.traj;
lindist = linpos{day}{ep}.statematrix.lindist;
timevec = timerange(1):tBinSz/1000:timerange(2);
%%
spikes = loaddatastruct(dir, animalprefix, 'spikes', day); % get spikes
decodinginfo = struct;
errorpos_all = [];
for i = 1:length(timevec)
    bin = timevec(i);
    avglindist_bin = find(lintimes >= (bin - tBinSz/2000) & lintimes < (bin + tBinSz/2000));
    traj_bin = traj(avglindist_bin);
    if ~isempty(traj_bin)
        if all(traj_bin-traj_bin(1) == 0) && (traj_bin(1) > 0)%valid bin
           current_tr = traj_bin(1);
           avglindist = nanmean(lindist(avglindist_bin));
           celldata = [];
           for n = 1:hpnum
              index = [day,ep,cellidxm(n,:)] ;
              if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
                  spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
              else
                  spiketimes = [];
              end
              if ~isempty(spiketimes)
                  spikebins = spiketimes(find(spiketimes >= (bin - tBinSz/2000) & spiketimes <= (bin + tBinSz/2000)));
                  spikecount(n) = length(spikebins);
                  tmpcelldata = [length(spikebins),n];
              else
                  tmpcelldata = [0,n];
                  spikecount(n) = 0;
              end
              celldata = [celldata;tmpcelldata];
           end
           cellcounts = sum((spikecount > 0));
           trajidx = find(trajinfo == current_tr);
           pm_tr = pm(trajidx,1);
           if (cellcounts > cellthresh) && (avglindist > 15) && (avglindist < max(pm_tr)-15) % exclude 15-cm round reward wells
              cellsi = find(spikecount > 0);
              rm_tr = rm(trajidx,cellsi);

              spkPerBin = celldata(cellsi,1)';
              expecSpk_tr = expecSpk(trajidx,cellsi);

              nPBin = size(rm_tr,1); %N positional bin
              wrking = bsxfun(@power, rm_tr, spkPerBin); %[nPos x nTbin x nCell]
              wrking = prod(wrking,2);
              expon = exp(-sum(expecSpk_tr,2)); %Exponent of equation.
              post = bsxfun(@times,wrking, expon); %[nPos x nTbin x nCell]
              post(isnan(post)) = 0;  
              post = post./sum(post);


              [~,id] = max(post);
              errorpos = abs(pm_tr(id) - avglindist);
              errorpos_all = [errorpos_all;errorpos];

              decodinginfo(i).bin = bin;
              decodinginfo(i).traj = current_tr;              
              decodinginfo(i).actualpos = avglindist;
              decodinginfo(i).errorpos = errorpos;
           else
              decodinginfo(i).bin = bin;
              decodinginfo(i).traj = NaN;              
              decodinginfo(i).actualpos = NaN;
              decodinginfo(i).errorpos = NaN;
           end

        else
           decodinginfo(i).bin = bin;
           decodinginfo(i).traj = NaN;              
           decodinginfo(i).actualpos = NaN;
           decodinginfo(i).errorpos = NaN;
        end
    else
       decodinginfo(i).bin = bin;
       decodinginfo(i).traj = NaN;              
       decodinginfo(i).actualpos = NaN;
       decodinginfo(i).errorpos = NaN;
    end
end

       


       
          

              
            
            


        
        
        
        
