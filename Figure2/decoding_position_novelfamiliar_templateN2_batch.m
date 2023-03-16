%---------------------------------------------------------------%
%  This is the main script for decoding animal's current        %
%  position using rate maps from N' -- Wenbo Tang (Jan 05, 2023)%
%---------------------------------------------------------------%
clc;
clear all
close all;
%%
famnov_dir = '/Volumes/NovelFamiliar/'; % set directory
animalprefix = 'JS34'; %current animal
day = 1;
eps = 4; % 2 = N; 4 = F; 6 = N'
savedata = 1; % save results?
% set animal directory
dir = sprintf('%s/%s_direct/',famnov_dir,animalprefix);
%%
% get trial information
trajinfo = loaddatastruct(dir, animalprefix, 'trajinfo', day); 
% get linfields
load(sprintf('%s%slinfields0%d.mat',dir,animalprefix,day)); % get linearized place fields
linfields{day}{eps} = linfields{day}{6};% replace linfields in the current epoch with these in epoch 6 (N' session)

for ep = eps
    %----find trials-----%
    trajbound = trajinfo{day}{ep}.trajbound;
    trajtime = trajinfo{day}{ep}.trajtime;
    
    for tr = 1:length(trajtime(:,1))
        decoded_trajtime = trajtime(tr,:);
        exclude_trajtime = decoded_trajtime;
     
        decodinginfo_seg_CA1 = decoding_position_CA1_novelfamiliar_cross(animalprefix,day,ep,linfields,decoded_trajtime);
        decodinginfo_seg = decoding_position_PFC_novelfamiliar_cross(animalprefix,day,ep,linfields,decoded_trajtime);

        decodinginfo_ep{tr} = decodinginfo_seg;
        decodinginfo_ep_CA1{tr} = decodinginfo_seg_CA1;

        clear decodinginfo_seg;clear decodinginfo_seg_CA1
    end
    decodinginfo_cross{ep} = decodinginfo_ep;
    decodinginfo_cross_CA1{ep} = decodinginfo_ep_CA1;

    clear decodinginfo_ep
end      

%%
if savedata
    save([animalprefix, 'decodinginfo_PFC_cross_template6_200ms_sm.mat'],'decodinginfo_cross')
    save([animalprefix, 'decodinginfo_CA1_cross_template6_200ms_sm.mat'],'decodinginfo_cross_CA1')
end
