%---------------------------------------------------------------%
%  This is the main script for plotting the confusion matrix    %
%  -- Wenbo Tang (Jan 05, 2023)                                 %
%---------------------------------------------------------------%
clc;
clear all
close all;
%%
animalprefix_list = {'AM2','JS17','JS21','ZT2','JS34'};
day = 1;
eps = 4;
usingN2 = 1; % 1, decoding using N' rate maps; 0, decoding using rate maps of the current session
CTXHP = 0; % CTXHP = 1, PFC; CTXHP = 0, CA1
famnov_dir = '/Volumes/NovelFamiliar/'; % set directory

% directory for saved decoding results
filedir = '/Tang_CellRep_2023/Figure2/Decodepos/';
%%
for theAnimal = 1:length(animalprefix_list)
    animalprefix = animalprefix_list{theAnimal};
    % set animal directory
    dir = sprintf('%s/%s_direct/',famnov_dir,animalprefix);
   
    % load decoding results
    if CTXHP
        if ~usingN2
            load([filedir,animalprefix,'decodinginfo_PFC_cross_4trajs_200ms_sm.mat'])
        else
            load([filedir,animalprefix,'decodinginfo_PFC_cross_template6_200ms_sm.mat'])
        end
    else
        if ~usingN2
            load([filedir,animalprefix,'decodinginfo_CA1_cross_4trajs_200ms_sm.mat'])
        else
            load([filedir,animalprefix,'decodinginfo_CA1_cross_template6_200ms_sm.mat'])
        end
        decodinginfo_cross  = decodinginfo_cross_CA1;
    end
    %%
    linpos = loaddatastruct(dir, animalprefix, 'linpos', day); % get linpos
    trajinfo = loaddatastruct(dir, animalprefix, 'trajinfo', day); % get trajectory information
    %%
    for ep = eps
        end_point = 190; %cm, max position range 
         %----find trials-----%
        trajbound = trajinfo{day}{ep}.trajbound;
        trajtime = trajinfo{day}{ep}.trajtime;
        errorep = [];
        
        decodinginfo = decodinginfo_cross{ep}{1};
        binnum = length(decodinginfo);
        %%
        ii = 1;
        while ii <= binnum
            tmp = decodinginfo(ii).posrange;
            if ~isempty(tmp)
                break
            end
            ii = ii + 1;
        end
        maxpos = decodinginfo(ii).posrange;
        posvec = decodinginfo(ii).posvec;
        posrange_valid = [maxpos(:,1)+16,maxpos(:,2)-15];% exclude the 15cm start and end
        posvec_valid = [];
        trackid = [];
        for track = 1:4
            id = find((posvec >= maxpos(track,1)+16) & posvec <= maxpos(track,1)+end_point);
            posvec_valid = [posvec_valid;posvec(id)];
            trackid = [trackid;length(posvec_valid)];
        end
        % valid position range for each trajecory type
        posvec_valid = (posvec_valid(1:4:end-3)+posvec_valid(2:4:end-2)+...
            +posvec_valid(3:4:end-1)+posvec_valid(4:4:end))./4;% combine 4 bins; 8cm
        %%
        % reset matrices
        actual_matrix = zeros(length(posvec_valid),length(posvec_valid));
        decoded_matrix = zeros(length(posvec_valid),length(posvec_valid));

        for tr = 1:length(trajtime(:,1))
            decodinginfo = decodinginfo_cross{ep}{tr};
            if length(decodinginfo) > 1
                binnum = length(decodinginfo);
                for i = 1:binnum
                    actualpos = round(decodinginfo(i).actualpos/2)*2;
                    if ~isnan(actualpos)
                        if isExcluded(actualpos,posrange_valid) 
                            decodedpos = decodinginfo(i).decodedpos;
                            [~,id1] = min(abs(posvec_valid - actualpos));
                            [~,id2] = min(abs(posvec_valid - decodedpos));
                            errorep = [errorep;abs(decodinginfo(i).actualpos-decodedpos)];
                            if ~isempty(id1) &&  ~isempty(id2) 
                                actual_matrix(id1,:) = actual_matrix(id1,:) + 1;
                                decoded_matrix(id1,id2) =  decoded_matrix(id1,id2)+1;
                            end
                        end
                    end
                end
            end
        end
        
        if theAnimal == 1 && ep == eps(1) 
            trackid_raw = trackid;
            posvec_valid_raw = posvec_valid;
            actual_matrix_all = actual_matrix;
            decoded_matrix_all = decoded_matrix;
        else
            % scale positions for each animal, as the position tracking is slightly
            % different for different animals
            id_all = [];
            for track = 1:4
                if track == 1
                    id = 1:round(trackid(track)/4);
                    id_match = 1:round(trackid_raw(track)/4);                
                else
                    id = round(trackid(track-1)/4)+1:round(trackid(track)/4);
                    id_match = round(trackid_raw(track-1)/4)+1:round(trackid_raw(track)/4);
                end
                id_adj = (id-min(id))./(max(id)-min(id)).*(max(id_match)-min(id_match)) + min(id_match);
                id_all = [id_all,id_adj];
            end
            [X,Y] = meshgrid(id_all,id_all);
            [Xp,Yp] = meshgrid(1:round(trackid_raw(track)/4),1:round(trackid_raw(track)/4));
            try
                temp = interp2(X,Y,actual_matrix,Xp,Yp);
            catch
                temp = interp2(X(1:end-1,1:end-1),Y(1:end-1,1:end-1),actual_matrix,Xp,Yp);
            end
            temp(isnan(temp)) = 0;
            try
                actual_matrix_all =...
                            actual_matrix_all +temp;
            catch
                  actual_matrix_all =...
                            actual_matrix_all +temp(1:end-1,1:end-1);
            end
            try
                temp2 = interp2(X,Y,decoded_matrix,Xp,Yp);
            catch
                temp2 = interp2(X(1:end-1,1:end-1),Y(1:end-1,1:end-1),decoded_matrix,Xp,Yp);
            end
            temp2(isnan(temp2)) = 0;
            try
                decoded_matrix_all =...
                                decoded_matrix_all +temp2(1:end-1,1:end-1);
            catch
                decoded_matrix_all =...
                                decoded_matrix_all +temp2;
            end
        end   

    end
end
%%
% compute the confusion matrix
confusion_matrix = decoded_matrix_all./actual_matrix_all;
confusion_matrix(find(isnan(confusion_matrix))) = 0;
%%
% plot results
figure('position',[100,100,350,250])
imagesc(confusion_matrix)
axis square
caxis([0 0.2])
m=100;
cm_viridis= inferno(m); % use the matplotlib toolbox
colormap(cm_viridis)
hold on
for track = 1:3
    plot([1,length(posvec_valid_raw)],[trackid_raw(track)/4+0.5,trackid_raw(track)/4+0.5],'w:','linewidth',1.5)
    plot([trackid_raw(track)/4+0.5,trackid_raw(track)/4+0.5],[1,length(posvec_valid_raw)],'w:','linewidth',1.5)
end
axis off
