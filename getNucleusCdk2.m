%% Nuclear mask and Cdk2 activity quantification for cycle 10-13
clearvars -except master_date master_series master_count isthisseriesconcatenated version
tStart = tic;
%% Nucleus mask and quantification of Cdk2
load('dateseries')
fprintf('getNucleusCdk2: %s / Series%s \n', date, series)

load(strcat('getHLB_', date, '_series', series), 'I_nondecon', 'I_4D_mask',...
    'index_begin', 'index_end', 'cc', 'M', 'ana', 'z_begin', 'z_end', 'zspan', 'z_middle', 'pixelxy', 'pixelz', 'cycle_time',...
    'series', 'series2', 'series3', 'conc11_12', 'conc12_13', 'start_cycle', 'end_cycle')

if exist(strcat('emask_', date, '_series', series, '.mat'), 'file')
    load(strcat('emask_', date, '_series', series), 'emask')
end

%% Load images
% Green channel
dir_nondecon = strcat('..\..\..\', date, '\Series', series); % Load non-deconvolved images
if ~isnan(series2), dir_nondecon2 = strcat('..\..\..\', date, '\Series', series2); end 
if ~isnan(series3), dir_nondecon3 = strcat('..\..\..\', date, '\Series', series3); end 

for n=1:12
    files_nondecon(:,n)= dir([dir_nondecon filesep strcat('*z', sprintf('%02d', n-1), '_ch00*')]);
    if ~isnan(series2)
        files_nondecon2(:,n)= dir([dir_nondecon2 filesep strcat('*z', sprintf('%02d', n-1), '_ch00*')]);
    end
    if ~isnan(series3)
        files_nondecon3(:,n)= dir([dir_nondecon3 filesep strcat('*z', sprintf('%02d', n-1), '_ch00*')]);
    end
end

if conc11_12 && ~conc12_13
    files_nondecon_comp(1:M{11}(2),:) = files_nondecon(1:M{11}(2),:);
    files_nondecon_comp(cc{12}(1):M{end_cycle}(2),:) = files_nondecon2(1:M{end_cycle}(2)-cc{12}(1)+1,:);
elseif ~conc11_12 && conc12_13
    files_nondecon_comp(1:M{12}(2),:) = files_nondecon(1:M{12}(2),:);
    files_nondecon_comp(cc{13}(1):M{13}(2),:) = files_nondecon2(1:M{13}(2)-cc{13}(1)+1,:);
elseif conc11_12 && conc12_13  % when you have 3 separate movies for each cycle
    files_nondecon_comp(1:M{11}(2),:) = files_nondecon(1:M{11}(2),:);
    files_nondecon_comp(cc{12}(1):M{12}(2),:) = files_nondecon2(1:M{12}(2)-cc{12}(1)+1,:);
    files_nondecon_comp(cc{13}(1):M{13}(2),:) = files_nondecon3(1:M{13}(2)-cc{13}(1)+1,:);
else
    files_nondecon_comp(1:M{end_cycle}(2),:) = files_nondecon(1:M{end_cycle}(2),:);
end

% Red channel
dir_red = strcat('..\..\..\', date, '\Series', series); 
if ~isnan(series2), dir_red2 = strcat('..\..\..\', date, '\Series', series2); end 
if ~isnan(series3), dir_red3 = strcat('..\..\..\', date, '\Series', series3); end 

for n=1:12
    files_red(:,n)= dir([dir_red filesep strcat('*z', sprintf('%02d', n-1), '_ch01*')]);
    if ~isnan(series2)
        files_red2(:,n)= dir([dir_red2 filesep strcat('*z', sprintf('%02d', n-1), '_ch01*')]);
    end
    if ~isnan(series3)
        files_red3(:,n)= dir([dir_red3 filesep strcat('*z', sprintf('%02d', n-1), '_ch01*')]);
    end
end

if conc11_12 && ~conc12_13
    files_red_comp(1:M{11}(2),:) = files_red(1:M{11}(2),:);
    files_red_comp(cc{12}(1):M{end_cycle}(2),:) = files_red2(1:M{end_cycle}(2)-cc{12}(1)+1,:);
elseif ~conc11_12 && conc12_13
    files_red_comp(1:M{12}(2),:) = files_red(1:M{12}(2),:);
    files_red_comp(cc{13}(1):M{13}(2),:) = files_red2(1:M{13}(2)-cc{13}(1)+1,:);
elseif conc11_12 && conc12_13 % when you have 3 separate movies for each cycle
    files_red_comp(1:M{11}(2),:) = files_red(1:M{11}(2),:);
    files_red_comp(cc{12}(1):M{12}(2),:) = files_red2(1:M{12}(2)-cc{12}(1)+1,:);
    files_red_comp(cc{13}(1):M{13}(2),:) = files_red3(1:M{13}(2)-cc{13}(1)+1,:);
else
    files_red_comp(1:M{end_cycle}(2),:) = files_red(1:M{end_cycle}(2),:);
end

z = z_middle; % Edit: just define in getHLB

for t=index_begin:index_end
    if conc11_12 && ~conc12_13
        if t<=M{11}(2)
            I=imread([dir_nondecon filesep files_nondecon_comp(t,z(t)).name]);
            J=imread([dir_red filesep files_red_comp(t,z(t)).name]);
        else
            I=imread([dir_nondecon2 filesep files_nondecon_comp(t,z(t)).name]);
            J=imread([dir_red2 filesep files_red_comp(t,z(t)).name]);
        end
    elseif ~conc11_12 && conc12_13
        if t<=M{12}(2)
            I=imread([dir_nondecon filesep files_nondecon_comp(t,z(t)).name]);
            J=imread([dir_red filesep files_red_comp(t,z(t)).name]);
        else
            I=imread([dir_nondecon2 filesep files_nondecon_comp(t,z(t)).name]);
            J=imread([dir_red2 filesep files_red_comp(t,z(t)).name]);
        end
    elseif conc11_12 && conc12_13 % when you have 3 separate movies for each cycle
        if t<=M{11}(2)
            I=imread([dir_nondecon filesep files_nondecon_comp(t,z(t)).name]);
            J=imread([dir_red filesep files_red_comp(t,z(t)).name]);
        elseif t<=M{12}(2)
            I=imread([dir_nondecon2 filesep files_nondecon_comp(t,z(t)).name]);
            J=imread([dir_red2 filesep files_red_comp(t,z(t)).name]);
        else
            I=imread([dir_nondecon3 filesep files_nondecon_comp(t,z(t)).name]);
            J=imread([dir_red3 filesep files_red_comp(t,z(t)).name]);
        end
    else
        I=imread([dir_nondecon filesep files_nondecon_comp(t,z(t)).name]);
        J=imread([dir_red filesep files_red_comp(t,z(t)).name]);
    end
    I_ori(:,:,t) = I;
    J_nondecon(:,:,t) = J;
    
%     maxI(t) = max(max(I));
    nucmeanI(t) = mean(mean(I.*uint16(1-(I_4D_mask(:,:,z_middle(t)-z_begin(t)+1,t))))); % Average Mxc level in nuclei, outside HLB
end

% for t=index_begin:index_end
%     for num=10:14
%         if ~isnan(cc{num}), if ismember(t,cc{num}(1):M{num}(2)), minmaxI(t) = min(min(maxI(cc{num}(1):M{num}(2)))); end, end
%         if ~isnan(cc{num}), if ismember(t,cc{num}(1):M{num}(2)), minmeanI(t) = min(min(meanI(cc{num}(1):M{num}(2)))); end, end
%     end
% end

for t=index_begin:index_end
    if conc11_12 && ~conc12_13
        if t<=M{11}(2)
            I=imread([dir_nondecon filesep files_nondecon_comp(t,z(t)).name]); % I = HLB
            J=imread([dir_red filesep files_red_comp(t,z(t)).name]); % J = Cdk2
        else
            I=imread([dir_nondecon2 filesep files_nondecon_comp(t,z(t)).name]);
            J=imread([dir_red2 filesep files_red_comp(t,z(t)).name]);
        end
    elseif ~conc11_12 && conc12_13
        if t<=M{12}(2)
            I=imread([dir_nondecon filesep files_nondecon_comp(t,z(t)).name]);
            J=imread([dir_red filesep files_red_comp(t,z(t)).name]);
        else
            I=imread([dir_nondecon2 filesep files_nondecon_comp(t,z(t)).name]);
            J=imread([dir_red2 filesep files_red_comp(t,z(t)).name]);
        end
    elseif conc11_12 && conc12_13 % when you have 3 separate movies for each cycle
        if t<=M{11}(2)
            I=imread([dir_nondecon filesep files_nondecon_comp(t,z(t)).name]);
            J=imread([dir_red filesep files_red_comp(t,z(t)).name]);
        elseif t<=M{12}(2)
            I=imread([dir_nondecon2 filesep files_nondecon_comp(t,z(t)).name]);
            J=imread([dir_red2 filesep files_red_comp(t,z(t)).name]);
        else
            I=imread([dir_nondecon3 filesep files_nondecon_comp(t,z(t)).name]);
            J=imread([dir_red3 filesep files_red_comp(t,z(t)).name]);
        end
    else
        I=imread([dir_nondecon filesep files_nondecon_comp(t,z(t)).name]);
        J=imread([dir_red filesep files_red_comp(t,z(t)).name]);
    end
    I = imgaussfilt(I, 4); % Blur to get smoother edge
    
    % I_nuc = nuclear mask
    if t>=cc{10}(1) && t<=M{10}(2)
        I_nuc = (I>1.00*nucmeanI(t));
    elseif t>=cc{11}(1) && t<=M{11}(2)
        I_nuc = (I>1.00*nucmeanI(t));
    elseif t>=cc{12}(1) && t<=M{12}(2)
        I_nuc = (I>1.00*nucmeanI(t));
    elseif t>=cc{13}(1) && t<=M{13}(2)
        I_nuc = (I>1.00*nucmeanI(t));
    elseif t>=cc{14}(1) && t<=M{14}(2)
        I_nuc = (I>0.80*nucmeanI(t));
    end
    if exist('emask', 'var'), I_nuc=I_nuc.*emask; end

    I_nuc = bwmorph(I_nuc, 'clean');
    
    for num_maj = 1:5
        I_nuc = bwmorph(I_nuc, 'majority');
    end
    
    I_nuc = imfill(I_nuc, 'holes');
    I_nuc = imerode(I_nuc,strel('disk',12,0)); % 15 for zoom of 10, 7 for zoom of 5
    I_nuc = imdilate(I_nuc,strel('disk',12,0)); % 15 for zoom of 10, 7 for zoom of 5
    
    I_center = imerode(I_nuc,strel('disk',12,0)); % 12 for zoom of 10, 6 for zoom of 5
    
    erosionnum = 20; % Pick a good number to make sure the nuclear masks to not overlap
    I_centroid = imerode(I_nuc,strel('disk',erosionnum,0)); 


    L = bwlabel(I_nuc);
    stats = regionprops(L, 'Area','Centroid');
    
    for i = 1:max(max(L))
        area = stats(i).Area;
        [ind1, ind2]= find(L==i);
        if(area<1000)
            I_nuc(ind1,ind2)=0;
        end
    end
    
    I_nucdilate=imdilate(I_nuc,strel('disk',4,0)); % 4 for zoom of 10, 2 for zoom of 5
    I_nucdilate2=imdilate(I_nuc,strel('disk',7,0)); % 7 for zoom of 10, 4 for zoom of 5
%     I_cytring=I_nucdilate2-I_nucdilate;
    if exist('emask', 'var')
        I_cytring=((I_nucdilate2-I_nucdilate)-emask)>0;
    else
        I_cytring=I_nucdilate2-I_nucdilate; % Use cytoplasmic ring of 4-7 pixels away from nuclear membrane
    end
    
    I_nucedge=imdilate(I_nuc,strel('disk',1,0))-I_nuc; % To display, not used for calculations

    I_centerdilate=imdilate(I_center,strel('disk',2,0));
    I_nucring=I_centerdilate-I_center;
    
    I_centroiddilate=imdilate(I_centroid,strel('disk',2,0));
    I_centroidring=I_centroiddilate-I_centroid;


    J_3D_centroid(:,:,t) = I_centroid;
    J_3D_mask(:,:,t) = I_nuc;
    
    area_nucleus = sum(sum(I_center));
    area_cytoplasm = sum(sum(I_cytring));
    cdk2_total_nucleus = sum(sum(I_center.*double(J)));
    cdk2_total_cytoplasm = sum(sum(I_cytring.*double(J)));
    
    cdk2_nucleus = cdk2_total_nucleus./area_nucleus;
    cdk2_cytoplasm = cdk2_total_cytoplasm./area_cytoplasm;
    
    cdk2_activity(t) = cdk2_cytoplasm./cdk2_nucleus;
    
    nucleus_wo_HLB = I_nuc-I_4D_mask(:,:,z_middle(t)-z_begin(t)+1,t);
    mxc_total = sum(sum(nucleus_wo_HLB.*double(I)));
    area_nucleus_wo_HLB = sum(sum(nucleus_wo_HLB));
    mxc_nuc_conc(t) = mxc_total./area_nucleus_wo_HLB;

    nucleus_whole = I_nuc;
    mxc_total_whole = sum(sum(nucleus_whole.*double(I_nondecon(:,:,z_middle(t)-z_begin(t)+1,t))));
    area_nucleus_whole = sum(sum(nucleus_whole));
    mxc_whole_conc(t) = mxc_total_whole./area_nucleus_whole;

    % Left = green channel, Right = red channel
    I_3D(:,:,1,t)= horzcat(50000*uint16(I_nucedge),50000*uint16(I_nucedge)+5*J_nondecon(:,:,t));
    I_3D(:,:,2,t)= horzcat(50000*uint16(I_nucedge)+3*I_ori(:,:,t),50000*uint16(I_nucedge));
    I_3D(:,:,3,t)= horzcat(50000*uint16(I_nucring+I_centroidring),50000*uint16(I_nucring+I_centroidring));
end

cdk2_activity(cdk2_activity==0)=NaN;
mov = I_3D(:,:,:,index_begin:index_end);

implay(mov)
set(findall(0,'tag','spcui_scope_framework'),'position',[150 150 1620 550]);

figure; clf;
hold on
for num=10:14
    if ~isnan(cc{num})
        cdk2_preM{num} = cdk2_activity(cc{num}(1):cc{num}(2));
        cdk2_postM{num} = cdk2_activity(M{num}(1):M{num}(2));
        fill([M{num}(1), M{num}(1), M{num}(2), M{num}(2)]*cycle_time/60, [0, 2.5, 2.5, 0], [.9 .9 .9], 'linestyle', 'none')
    end
end
plot((index_begin:index_end)*cycle_time/60, cdk2_activity(index_begin:index_end))
ylabel('Cdk2 Activity'), xlabel('Time (min)')
set(gca, 'ylim', [0.6 2.3])
set(gca,'layer','top')

for num=10:14
    if ~isnan(cc{num})
        mxc_nuc{num} = mxc_nuc_conc(cc{num}(1):M{num}(2));
        mxc_whole{num} = mxc_whole_conc(cc{num}(1):M{num}(2));
    else
        mxc_nuc{num} = NaN;
        mxc_whole{num} = NaN;
    end
end

% figure; clf;
% hold on
% for num=10:14
%     if ~isnan(cc{num})
%         fill([M{num}(1), M{num}(1), M{num}(2), M{num}(2)]*cycle_time/60, [0.9*min(mxc_nuc_conc(mxc_nuc_conc>0)), 1.1*max(mxc_nuc_conc), 1.1*max(mxc_nuc_conc), 0.9*min(mxc_nuc_conc(mxc_nuc_conc>0))], [.9 .9 .9], 'linestyle', 'none')
%     end
% end
% plot((index_begin:index_end)*cycle_time/60, mxc_nuc_conc(index_begin:index_end))
% ylabel('Mxc Concentration'), xlabel('Time (min)')
% set(gca, 'ylim', [0.9*min(mxc_nuc_conc(mxc_nuc_conc>0)) 1.1*max(mxc_nuc_conc)])
% set(gca,'layer','top')
% 
% figure; clf;
% hold on
% for num=10:14
%     if ~isnan(cc{num})
%         fill([M{num}(1), M{num}(1), M{num}(2), M{num}(2)]*cycle_time/60, [0.9*min(mxc_whole_conc(mxc_whole_conc>0)), 1.1*max(mxc_whole_conc), 1.1*max(mxc_whole_conc), 0.9*min(mxc_whole_conc(mxc_whole_conc>0))], [.9 .9 .9], 'linestyle', 'none')
%     end
% end
% plot((index_begin:index_end)*cycle_time/60, mxc_whole_conc(index_begin:index_end))
% ylabel('Total Mxc Concentration'), xlabel('Time (min)')
% set(gca, 'ylim', [0.9*min(mxc_whole_conc(mxc_whole_conc>0)) 1.1*max(mxc_whole_conc)])
% set(gca,'layer','top')

tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n',floor(tEnd/60),rem(tEnd,60));

%%  
save(strcat('getNucleusCdk2_', date, '_series', series), 'J_nondecon', 'cdk2_activity', 'cdk2_preM', 'cdk2_postM', 'J_3D_mask',...
    'J_3D_centroid', 'erosionnum', 'z_middle', 'mxc_nuc', 'mxc_nuc_conc', 'mxc_whole', 'mxc_whole_conc')