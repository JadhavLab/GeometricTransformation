function [ctxidx, hpidx] = matchidx_acrossep(animaldir, animalprefix, day, exclude_list)
eps = [2 4 6];

%%
tetinfo = loaddatastruct(animaldir, animalprefix, 'tetinfo');
cell = loaddatastruct(animaldir,animalprefix,'cellinfo'); 
spks = loaddatastruct(animaldir, animalprefix, 'spikes', day); % get spikes

%%
for ep = eps
    % get the tetrode number
    detect_tet = [];
    for ttt = 1:length(tetinfo{day}{ep})
                if ~isempty(tetinfo{day}{ep}{ttt})
                    if isfield(tetinfo{day}{ep}{ttt},'area')
                        if strcmp(tetinfo{day}{ep}{ttt}.area,'PFC') % Cortical neurons first
                            detect_tet = [detect_tet ttt];
                        elseif isfield(tetinfo{day}{ep}{ttt},'descrip')
                            if strcmp(tetinfo{day}{ep}{ttt}.descrip,'riptet')
                                detect_tet = [detect_tet ttt];
                            end
                        end
                    end
                end
    end
    % get the cell number
    % hp
    cellind = []; % (tet, unit)
    for tet = detect_tet  % tetrode iterations
        spikes = [];
        cind = [];
        if max(length(spks{day}{ep}) >= tet)
            spikes = spks{day}{ep}{tet};
        end
        if ~isempty(spikes)
            if  strcmp(tetinfo{day}{ep}{tet}.area,'CA1') || strcmp(tetinfo{day}{ep}{tet}.area,'iCA1')  % CA1 neurons only
                hpneurons = find(~cellfun('isempty',spikes));
                nn = 0;
                for neuron = hpneurons
                    if ~isempty(spikes{neuron}.data)
%                         if isfield(cell{day}{ep}{tet}{neuron},'tag2') %no tag2 for now
%                             if strcmp(cell{day}{ep}{tet}{neuron}.tag2,'CA1Pyr') || strcmp(cell{day}{ep}{tet}{neuron}.tag2,'iCA1Pyr') % putative pyramidal neurons only
    %                         if strcmp(cell{day}{ep}{tet}{neuron}.tag2,'CA1Pyr')
                                if isfield(cell{day}{ep}{tet}{neuron},'meanrate')
                                    if (cell{day}{ep}{tet}{neuron}.meanrate < 7)
                                        if isfield(cell{day}{ep}{tet}{neuron},'numspikes')
                                            if (cell{day}{ep}{tet}{neuron}.numspikes > 100)
                                                nn = nn+1;
                                                cind(nn) = neuron;
                                            end
                                        end
                                    end
                                end
%                             end
%                         end
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
    % use exclude list
    if ~isempty(exclude_list)
        cellind_raw = cellind;
        cellind = setdiff(cellind_raw,exclude_list,'rows');
    end
    %ctx
    cellind_ctx = [];
    for tet = detect_tet  % tetrode iterations
        spikes = [];
        cind = [];
        if max(length(spks{day}{ep}) >= tet)
            spikes = spks{day}{ep}{tet};
        end
        if ~isempty(spikes)
            if  strcmp(tetinfo{day}{ep}{tet}.area,'PFC')  % PFC neurons only
                ctxneurons = find(~cellfun('isempty',spikes));
                nn = 0;
                for neuron = ctxneurons
                    if ~isempty(spikes{neuron}.data)
                        if isfield(cell{day}{ep}{tet}{neuron},'tag2')
    %                         if strcmp(cell{day}{ep}{tet}{neuron}.tag2,'CA1Pyr') || strcmp(cell{day}{ep}{tet}{neuron}.tag2,'iCA1Pyr') % putative pyramidal neurons only
                            if strcmp(cell{day}{ep}{tet}{neuron}.tag2,'PFC')
                                    if isfield(cell{day}{ep}{tet}{neuron},'numspikes')
                                        if (cell{day}{ep}{tet}{neuron}.numspikes > 100)
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
        cellind_ctx = [cellind_ctx;temp];
    end
    % use exclude list
    if ~isempty(exclude_list)
        cellind_raw = cellind_ctx;
        cellind_ctx = setdiff(cellind_raw,exclude_list,'rows');
    end
    ctxidx_ep{ep} = cellind_ctx;
    hpidx_ep{ep} = cellind;
    
    % match across epochs
    if ep == eps(1)
        ctxidx_temp = ctxidx_ep{eps(1)}; % use the first ep as template
        hpidx_temp = hpidx_ep{eps(1)};
    else
        if ~isempty(ctxidx_ep{ep}) && ~isempty(hpidx_ep{ep})
            ctxmatchind = ctxidx_ep{ep};
            hpmatchind = hpidx_ep{ep};
            ctx_rowidx_temp = zeros(1,length(ctxmatchind(:,1)));
            hp_rowidx_temp = zeros(1,length(hpmatchind(:,1)));
            for i = 1:length(ctxmatchind(:,1))
                if ~isempty(ctxidx_temp)
                   if ~isempty(find(ismember(ctxidx_temp,ctxmatchind(i,:),'rows'),1))
                       ctx_rowidx_temp(i) = find(ismember(ctxidx_temp,ctxmatchind(i,:),'rows'),1);
                   else
                       ctx_rowidx_temp(i) = 0;
                   end
                else
                       ctx_rowidx_temp(i) = 0;
                end
            end
            ctxrowidx{ep} = ctx_rowidx_temp;
            for i = 1:length(hpmatchind(:,1))
                if ~isempty(hpidx_temp)
                   if ~isempty(find(ismember(hpidx_temp,hpmatchind(i,:),'rows'),1))
                       hp_rowidx_temp(i) = find(ismember(hpidx_temp,hpmatchind(i,:),'rows'),1);
                   else
                       hp_rowidx_temp(i) = 0;
                   end
                else
                       hp_rowidx_temp(i) = 0;
                end
            end
            hprowidx{ep} = hp_rowidx_temp;
        end
    end    
end
%%
ctxi1 = ctxrowidx{eps(2)};
hpi1 = hprowidx{eps(2)};
for i = 2:length(eps)-1
    ctxi1 = intersect(ctxi1,ctxrowidx{eps(i+1)});
    hpi1 = intersect(hpi1,hprowidx{eps(i+1)});
end
ctxindstart = min(find(ctxi1 > 0));
hpindstart = min(find(hpi1 > 0));
ctxidx = ctxidx_temp(ctxi1(ctxindstart:end),:);  %(tet, cell)
hpidx = hpidx_temp(hpi1(hpindstart:end),:);  %(tet, cell)

    
