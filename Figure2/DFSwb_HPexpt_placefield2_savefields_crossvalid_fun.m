function linfields = DFSwb_HPexpt_placefield2_savefields_crossvalid_fun(animalprefix,ep,exclude_trajtime, savedir)

%---------------------------------------------------------------%
%  This is the function for calculating linearized place fields %
%  for cross-validation                                         %
%  Previously written by Shantanu, modified by Wenbo            %
%  -- Wenbo Tang (Sep 11, 2019)                                 %
%---------------------------------------------------------------%

% INPUTS:
%
%    animalprefix = animal prefix
%    ep = epoch. If multiday training is used, modifications will be needed (using day as an additional inputs).
%    exclude_trajtime = periods not used for estimation, usually the testing trial
%    savedir = directory for saving files
%
% OUTPUTS:
%
%    linfields are the estimated linearized place fields. 
%%
% codes below used the time filter framework

% set parameters
runscript = 1;
savelinfields = 0; % To save trajdata in [prefix-linfields-day] file for each day

minabsvel = 3;  % cm/sec - Most conservative for runs and place fields
minlinvel = 5;
savedata = 0;

if runscript == 1
    %Animal selection
    %-----------------------------------------------------
    animals = {animalprefix};
    %--------------------------------------------------------
    
    % Epoch filter
    % -------------
    runepochfilter = 'isequal($environment, ''wtr1'') ';
    
    % Cell filter
    % -----------
    % All cells with Nspk consition
    % ------
    cellfilter = '( strcmp($area, ''CA1'') && ($numspikes > 100) ) ||  (strcmp($area, ''PFC'') && ($numspikes > 100) ) ';
    % Time filter
    % -----------
    
    % Either use tetlist for riple detection for ripfilter,
    % or use ripple tet filter based on tag in tetinfo marking elec as ripple tet
    
    riptetfilter = '(isequal($descrip, ''riptet''))';
    timefilter = { {'DFTFsj_getlinstate', '(($state ~= -1) & (abs($linearvel) >= 5))', 6}, {'DFTFsj_getriptimes','($nripples == 0)','tetfilter',riptetfilter,'minthresh',3}};
    
    % Iterator
    % --------
    iterator = 'singlecellanal';
    
    % Filter creation
    % ----------------
    psf = createfilter('animal',animals,'epochs',runepochfilter,'cells',cellfilter,'excludetime', timefilter,'iterator', iterator);
    
    psf.epochs{1} = [1,ep];
    tmp = psf.excludetime{1, 1}{ep/2};
    tmp = [tmp;exclude_trajtime];
    psf.excludetime{1}{1} = tmp;
    psf.data{1}{1}  =  psf.data{1}{ep/2};
    psf.excludetime{1}(2:end) = [];
    psf.data{1}(2:end)  = [];
    
    % Set analysis function
    % ----------------------
    psf = setfilterfunction(psf, 'DFAsj_filtercalclinfields_tf', {'spikes', 'linpos'}, 'binsize', 2);
    
    
    disp('Finished filter creation');
    
    % Run analysis
    % ------------
    psf = runfilter(psf);  % Place Field Stability
    
    disp('Finished running filter script');
    %--------------------- Finished Filter Function Run -------------------
    
    if savedata == 1
        clear figopt1 runscript plotdays plotanimidx savedata savelinfields
        save(savefile);
    end
end  
%%
% Get trajdata and days and epochs
trajdata = []; index=[];
allanimindex=[]; allmaps=[]; alltrajs=[];

for an = 1:length(psf)
    for i=1:length(psf(an).output{1}),
        index{an}(i,:)=psf(an).output{1}(i).index;
        alltrajdata{an}{i}=psf(an).output{1}(i).trajdata;
        % Only indexes
        animindex=[an psf(an).output{1}(i).index]; % Put animal index in front
        allanimindex = [allanimindex; animindex]; % Collect all Anim Day Epoch Tet Cell Index
        alltrajs{i} = psf(an).output{1}(i).trajdata;
    end
end

for an = 1:length(psf)
    prefix = psf(an).animal{1};
    animdirect = psf(an).animal{2};
    if (animdirect(end) == '/')
        animdirect = animdirect(1:end-1);
    end
    
    curranimidxs = index{an}; % All indexes for cells from current animal. day-ep-tet-cell
    uniquedays = unique(curranimidxs(:,1));
    
    for d = 1:length(uniquedays)
        day = uniquedays(d);
        dayidxs = find(curranimidxs(:,1)==day);
        curranimdayidxs = curranimidxs(dayidxs,:);
        linfields = []; 
        for c = 1:length(curranimdayidxs)
            curridx = curranimdayidxs(c,:);
            currmapidx = find( index{an}(:,1)==curridx(1) & index{an}(:,2)==curridx(2) & index{an}(:,3)==curridx(3) & index{an}(:,4)==curridx(4));
            linfields{curridx(1)}{curridx(2)}{curridx(3)}{curridx(4)}=alltrajdata{an}{currmapidx};
        end
        
        % Save for current day
         if savelinfields==1
            savefile = sprintf('%s/%slinfields_crossvalid%02d.mat', savedir, prefix, day);
            save(savefile,'linfields');
         end  
    end % end day
end % end an

    
