function decodinginfo = decoding_position_CA1_novelfamiliar_cross(animalprefix,day,ep,linfields,decoded_trajtime)
%---------------------------------------------------------------%
%  This is the function for CA1 decoding of spatial locations   %
%  during the test trial (cross-validation procedure)           %
%  -- Wenbo Tang (Jan 05, 2023)                                 %
%---------------------------------------------------------------%

% INPUTS:
%
%    animalprefix = animal prefix.
%    day = experimental day.
%    ep = epoch.
%    linfields = linearized place fields with a standardized structure (leave out the testing trial).
%    decoded_bins = time bins to be decoded during the testing trial
%
% OUTPUTS:
%
%    decodinginfo for the time bins to be decoded during the testing trial
%%
% set animal directory
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

tBinSz = 50; % moving window, in ms
tBinSz_sm = 200; %default temporal bin in ms;

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
%-----get cell information for the current epoch-----%
tetinfo = loaddatastruct(dir, animalprefix, 'tetinfo');
cellinfo = loaddatastruct(dir,animalprefix,'cellinfo'); 
spks = loaddatastruct(dir, animalprefix, 'spikes', day); % get spikes
% current epoch
ep_curr = ep;
detect_tet = [];
for ttt = 1:length(tetinfo{day}{ep_curr})
     if ~isempty(tetinfo{day}{ep_curr}{ttt})
        if isfield(tetinfo{day}{ep_curr}{ttt},'area')
            if isfield(tetinfo{day}{ep_curr}{ttt},'descrip')
                if strcmp(tetinfo{day}{ep_curr}{ttt}.descrip,'riptet')
                    detect_tet = [detect_tet ttt];
                end
            end
        end
     end
end
cellind = []; % (tet, unit)
for tet = detect_tet  % tetrode iterations
    spikes = [];
    cind = [];
    if max(length(spks{day}{ep_curr}) >= tet)
        spikes = spks{day}{ep_curr}{tet};
    end
    if ~isempty(spikes)
        if  strcmp(tetinfo{day}{ep_curr}{tet}.area,'CA1')
            hpneurons = find(~cellfun('isempty',spikes));
            nn = 0;
            for neuron = hpneurons
                if ~isempty(spikes{neuron}.data)
                    if isfield(cellinfo{day}{ep_curr}{tet}{neuron},'meanrate')
                       if (cellinfo{day}{ep_curr}{tet}{neuron}.meanrate < 7)
                           if isfield(cellinfo{day}{ep_curr}{tet}{neuron},'numspikes')
                               if (cellinfo{day}{ep_curr}{tet}{neuron}.numspikes > 100)
                                   nn = nn+1;
                                   cind(nn) = neuron;
                               end
                           end
                        end
                    end
                end
            end
        else
                cind = [];
        end
     else
        cind = [];
    end
    temp = [tet.*ones(length(cind),1),cind'];
    cellind = [cellind;temp];
end
%%
% exclude interneurons
if ~isempty(exclude_list)
    hpidx = setdiff(cellind,exclude_list,'rows');
else
    hpidx = cellind;
end
hpnum = length(hpidx(:,1));
%%
%-----create the ratemaps [nPosBin x nHPCells]-----%
rm = []; % ratemap matrix
pm = []; % position matrix
pm_all = []; % position matrix
tm = []; % track matrix
cellidxm = []; % cell indices

for i = 1:hpnum
      cind = hpidx(i,:);
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
           linfield_hp = [];
           lintrack_hp = [];
           pos_hp = [];
           pos_hp_all = [];
           posrange = [];
           tmp = [0,0];
           for track = 1:4
                temp1 = linfield1{track};
                pos1 = temp1(:,1);
                lintrack1 = ones(size(pos1))*track;
                occnormrate1 = temp1(:,5);
                linfield_hp = [linfield_hp;occnormrate1];
                posrange = [posrange;tmp+[pos1(1),pos1(end)]];
                tmp = [posrange(end,2),posrange(end,2)];
                pos_hp_all = [pos_hp_all;pos1+posrange(track,1)];
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
               pm_all = [pm_all;pos_hp_all'];
               tm = [tm;lintrack_hp'];
               cellidxm = [cellidxm; cind];
           end
      end
end
rm = rm'; %[nPosBin x nHPCells]
pm = pm';
pm_all = pm_all';
tm = tm';
trajinfo = mean(tm,2);% trajactory number
rm = rm+ (eps.^8); %Add a small number so there are no zeros
expecSpk =rm.*tBinSz_sm./1000; %[nPos x nCells] Expected number of spikes per bin
hpnum = length(rm(1,:)); % update
%%
%-----load spikes, linpos and trajinfo-----%
linpos = loaddatastruct(dir, animalprefix, 'linpos', day); % get linpos
pos = loaddatastruct(dir, animalprefix, 'pos', day); % get pos

lintimes = linpos{day}{ep}.statematrix.time;
traj = linpos{day}{ep}.statematrix.traj;
lindist = linpos{day}{ep}.statematrix.lindist;
timevec = lintimes(1):tBinSz/1000:lintimes(end);

load(sprintf('%s%sthetatime0%d.mat',dir,animalprefix,day)); % get theta time
thetaep = thetatime{day}{ep};
thetalist(:,1) = thetaep.starttime;
thetalist(:,2) = thetaep.endtime;
thetavec = list2vec_binedges(thetalist,timevec,tBinSz/2000);
timevec = timevec(find(thetavec == 1));
%%
%-----decode animal's current position in each time bin-----%
decoded_bins = find(timevec >= decoded_trajtime(1) & timevec <= decoded_trajtime(2)); 

spikes = loaddatastruct(dir, animalprefix, 'spikes', day); % get spikes
decodinginfo = struct;
for i = 1:length(decoded_bins)
    bin = timevec(decoded_bins(i));
    avglindist_bin = find(lintimes >= (bin - tBinSz_sm/2000) & lintimes < (bin + tBinSz_sm/2000));
    traj_bin = traj(avglindist_bin);
    if all(traj_bin-traj_bin(1) == 0) && (traj_bin(1) > 0)%valid bin
       current_tr = traj_bin(1);
       avglindist = nanmean(lindist(avglindist_bin));
       avglindist_all = avglindist + posrange(current_tr,1);
       celldata = [];
       for n = 1:hpnum
          index = [day,ep,cellidxm(n,:)] ;
          if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
              spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
          else
              spiketimes = [];
          end
          if ~isempty(spiketimes)
              spikebins = spiketimes(find(spiketimes >= (bin - tBinSz_sm/2000) & spiketimes <= (bin + tBinSz_sm/2000)));
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
       pm_tr_all = pm_all(:,1);
       if (cellcounts > 0) && (avglindist > 15) && (avglindist < max(pm_tr)-15) % exclude 15-cm round reward wells
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
          
          rm_tr = rm(:,cellsi);
          nPBin = size(rm_tr,1); %N positional bin
          wrking = bsxfun(@power, rm_tr, spkPerBin); %[nPos x nTbin x nCell]
          wrking = prod(wrking,2);
          expon = exp(-sum(expecSpk,2)); %Exponent of equation.
          post_all = bsxfun(@times,wrking, expon); %[nPos x nTbin x nCell]
          post_all(isnan(post_all)) = 0;  
          post_all = post_all./sum(post_all);
          [~,id] = max(post_all);
          
          errorpos = pm_tr_all(id) - avglindist_all;
          decodinginfo(i).bin = bin;
          decodinginfo(i).traj = current_tr;              
          decodinginfo(i).actualpos = avglindist_all;
          decodinginfo(i).decodedpos = pm_tr_all(id);
          decodinginfo(i).error = errorpos;
          decodinginfo(i).pMat = post;
          decodinginfo(i).pMat_all = post_all;
          decodinginfo(i).posrange = posrange;
          decodinginfo(i).posvec = pos_hp_all;

       else
          decodinginfo(i).bin = bin;
          decodinginfo(i).traj = NaN;              
          decodinginfo(i).actualpos = NaN;
          decodinginfo(i).decodedpos = NaN;
          decodinginfo(i).error = NaN;
          decodinginfo(i).pMat = [];
          decodinginfo(i).pMat_all = [];
          decodinginfo(i).posrange = [];
          decodinginfo(i).posvec = [];
       end

    else
       decodinginfo(i).bin = bin;
       decodinginfo(i).traj = NaN;              
       decodinginfo(i).actualpos = NaN;
       decodinginfo(i).decodedpos = NaN;
       decodinginfo(i).error = NaN;
       decodinginfo(i).pMat = [];
       decodinginfo(i).pMat_all = [];
       decodinginfo(i).posrange = [];
       decodinginfo(i).posvec = [];
    end
end

       


       
          

              
            
            


        
        
        
        
