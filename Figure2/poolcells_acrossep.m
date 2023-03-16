function [ctxidx, hpidx] = poolcells_acrossep(animaldir, animalprefix, day, eps,exclude_list)

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
            if  strcmp(tetinfo{day}{ep}{tet}.area,'CA1')
                hpneurons = find(~cellfun('isempty',spikes));
                nn = 0;
                for neuron = hpneurons
                    if ~isempty(spikes{neuron}.data)
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
    

    ctxidx_ep{ep} = cellind_ctx;
    hpidx_ep{ep} = cellind;
   
end
%%
% all cells
cellind_all = [];
cellind_all_ctx = [];
for ep = eps
    cellind_all = [cellind_all;hpidx_ep{ep}];
    cellind_all_ctx = [cellind_all_ctx;ctxidx_ep{ep}];
end
cellind_all = unique(cellind_all,'rows');
cellind_all_ctx = unique(cellind_all_ctx,'rows');
% use exclude list
if ~isempty(exclude_list)
    hpidx = setdiff(cellind_all,exclude_list,'rows');
    ctxidx = setdiff(cellind_all_ctx,exclude_list,'rows');
else
    hpidx = cellind_all;
    ctxidx = cellind_all_ctx;
end
    
