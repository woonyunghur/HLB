%% FRAP analysis for one embryo, one lif file
clear; tStart = tic;
date = '170317_FRAP_2';

pixelxy = 0.073; %um
cycle_time = 1.263; %sec

%% Load images
% Up to which one is each cycle
cc{11} = NaN;
cc{12} = 015;
cc{13} = 025;
cc{14} = NaN;

for num=11:14
    cycle_check(num)=num*~isnan(cc{num}(1));
end
start_cycle=min(cycle_check(cycle_check>0));
end_cycle=max(cycle_check(cycle_check>0));

dir_FRAP = '..';

index=0;
for n=0:cc{end_cycle}
    if ~isempty(dir([dir_FRAP filesep strcat('*FRAP_', sprintf('%03d', n), '*')]))
        index=index+1;
        for num=11:14
            if n==cc{num}
                cc_index(num)=index;
            end
        end
        files_pre_green(:,index) = dir([dir_FRAP filesep strcat('*FRAP_', sprintf('%03d', n), '*pre*t0_z0_ch00*')]);
        files_bleach(:,index) = dir([dir_FRAP filesep strcat('*FRAP_', sprintf('%03d', n), '*Bleach*t0_z0_ch00*')]);
        files_pb_green(:,index) = dir([dir_FRAP filesep strcat('*FRAP_', sprintf('%03d', n), '*Pb*z0_ch00*')]);
    end
end
index_length = index; % How many data sets in total
pb_length = size(files_pb_green,1);

% For each FRAP set, define which cell cycle it is (stored in mcc)
for m=1:index_length
    if m<=cc_index(11)
        mcc(m)=11;
    elseif m<=cc_index(12)
        mcc(m)=12;
    elseif m<=cc_index(13)
        mcc(m)=13;
    else
        mcc(m)=14;
    end
end

for m=1:index_length
    I_bleach{m}(:,:)=imread([dir_FRAP filesep files_bleach(:,m).name]);
    for t=1:pb_length+1
        if t==1
            I_pre=imread([dir_FRAP filesep files_pre_green(:,m).name]);
            I_comp{m}(:,:,t)=I_pre;
        else
            I_pb=imread([dir_FRAP filesep files_pb_green(t-1,m).name]);
            I_comp{m}(:,:,t)=I_pb;
        end
    end
    maxI{m}=max(max(I_comp{m}(:,:,1)));
    
    % Bleach mask
    I=I_bleach{m};
    for i=1:size(I,1)
        for j=1:size(I,2)
            threshold=10000;
            if(I(i,j)>threshold)
                I2(i,j)=1;
            else
                I2(i,j)=0;
            end
        end
    end
    I2 = bwmorph(I2, 'clean');
    for num_maj = 1:5
        I2 = bwmorph(I2, 'majority');
    end
    I2 = imfill(I2, 'holes');
    I_dilate=imdilate(I2,strel('disk',2,0));
    I_ring=I_dilate-I2;
    L = bwlabel(I2);
    stats = regionprops(L, 'Area','Centroid');
    
    for i = 1:max(max(L))
        area = stats(i).Area;
        [ind1, ind2]= find(L==i);
        if(area<10)
            I2(ind1,ind2)=0;
        end
    end

    I_bleach_mask{m} = I2;
    I_ring = 50000*uint16(I_ring);
    
    I_3D{m}(:,:,1,2)= I_ring;
    I_3D{m}(:,:,2,2)= 2*I_bleach{m};
    I_3D{m}(:,:,3,2)= 0;
    
    % HLB mask
    if m<=cc_index(11)
        threshold=0.13*maxI{m};
    elseif m<=cc_index(12)
        threshold=0.14*maxI{m};
    elseif m<=cc_index(13)
        threshold=0.15*maxI{m};
    else %cc14
        threshold=0.13*maxI{m};
    end
    for t=1:pb_length+1
        I=I_comp{m}(:,:,t);
        I = imgaussfilt(I, 0.5); % Blur to get smoother edge
        
        for i=1:size(I,1)
            for j=1:size(I,2)
                if(I(i,j)>threshold)
                    I2(i,j)=1;
                else
                    I2(i,j)=0;
                end
            end
        end
        
        I2 = bwmorph(I2, 'clean');
        for num_maj = 1:5
            I2 = bwmorph(I2, 'majority');
        end
        I2 = imfill(I2, 'holes');
        
        I2 = imdilate(I2,strel('disk',2,0));
        I2 = imerode(I2,strel('disk',2,0));
        
        %             [B,L2] = bwboundaries(I2, 'noholes');
        %             stats2 = regionprops(L2,'Area','Centroid');
        %             threshold = 0.85;
        %             for k=1:length(B)
        %                 boundary = B{k};
        %                 delta_sq = diff(boundary).^2;
        %                 perimeter = sum(sqrt(sum(delta_sq,2)));
        %                 area2 = stats2(k).Area;
        %                 metric = 4*pi*area2/perimeter^2;
        %                 [ind3, ind4] = find(L2==k);
        %                 if metric<threshold
        %                     I2(ind3,ind4)=0;
        %                 end
        %             end
        
                    L = bwlabel(I2);
                    stats = regionprops(L, 'Area','Centroid');
        
                    for i = 1:max(max(L))
                        area = stats(i).Area;
                        [ind1, ind2]= find(L==i);
                        if(area<10)
                            I2(ind1,ind2)=0;
                        end
                    end
        
        I_dilate=imdilate(I2,strel('disk',2,0));
        I_ring=I_dilate-I2;
        
        I_mask{m}(:,:,t) = I2;

        if t==1
            I_3D{m}(:,:,1,t)= 50000*uint16(I_ring);
            I_3D{m}(:,:,2,t)= imgaussfilt(3*I_comp{m}(:,:,t), 0.5);
            I_3D{m}(:,:,3,t)= 0;
        else
            I_3D{m}(:,:,1,t+1)= 50000*uint16(I_ring);
            I_3D{m}(:,:,2,t+1)= imgaussfilt(3*I_comp{m}(:,:,t), 0.5);
            I_3D{m}(:,:,3,t+1)= 0;
        end
    end

    mov{m} = squeeze(I_3D{m});
    mm=implay(mov{m}); set(findall(0,'tag','spcui_scope_framework'),'position',[150 150 820 550]);
    set(mm.Parent, 'Name', strcat('cc', num2str(mcc(m)), ' Set #', num2str(m)))
end



%%
save(strcat('getHLB_', date), 'pixelxy', 'cycle_time', 'cc', 'start_cycle', 'end_cycle',...
    'index_length', 'pb_length', 'mcc', 'I_bleach', 'I_comp', 'I_bleach_mask', 'I_mask', 'maxI')

%%
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n',floor(tEnd/60),rem(tEnd,60));



