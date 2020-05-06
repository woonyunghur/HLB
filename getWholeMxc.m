%% Nuclear mask and Cdk2 activity quantification for cycle 10-13
clearvars -except master_date master_series master_count isthisseriesconcatenated version
tStart = tic;
%% Nucleus mask and quantification of Cdk2
load('dateseries')
fprintf('getWholeMxc: %s / Series%s \n', date, series)

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

for t=index_begin:index_end
    if z_middle(t)-3 >= z_begin(t)
        z_quantify_begin(t) = z_middle(t)-3;
    else
        z_quantify_begin(t) = z_begin(t);
    end
    if z_middle(t)+3 <= z_end(t)
        z_quantify_end(t) = z_middle(t)+3;
    else
        z_quantify_end(t) = z_end(t);
    end
end

for t=index_begin:index_end
    for z=z_quantify_begin(t):z_quantify_end(t)
        if conc11_12 && ~conc12_13
            if t<=M{11}(2)
                I=imread([dir_nondecon filesep files_nondecon_comp(t,z).name]); % I = HLB
                J=imread([dir_red filesep files_red_comp(t,z).name]); % J = Cdk2
            else
                I=imread([dir_nondecon2 filesep files_nondecon_comp(t,z).name]);
                J=imread([dir_red2 filesep files_red_comp(t,z).name]);
            end
        elseif ~conc11_12 && conc12_13
            if t<=M{12}(2)
                I=imread([dir_nondecon filesep files_nondecon_comp(t,z).name]);
                J=imread([dir_red filesep files_red_comp(t,z).name]);
            else
                I=imread([dir_nondecon2 filesep files_nondecon_comp(t,z).name]);
                J=imread([dir_red2 filesep files_red_comp(t,z).name]);
            end
        elseif conc11_12 && conc12_13 % when you have 3 separate movies for each cycle
            if t<=M{11}(2)
                I=imread([dir_nondecon filesep files_nondecon_comp(t,z).name]);
                J=imread([dir_red filesep files_red_comp(t,z).name]);
            elseif t<=M{12}(2)
                I=imread([dir_nondecon2 filesep files_nondecon_comp(t,z).name]);
                J=imread([dir_red2 filesep files_red_comp(t,z).name]);
            else
                I=imread([dir_nondecon3 filesep files_nondecon_comp(t,z).name]);
                J=imread([dir_red3 filesep files_red_comp(t,z).name]);
            end
        else
            I=imread([dir_nondecon filesep files_nondecon_comp(t,z).name]);
            J=imread([dir_red filesep files_red_comp(t,z).name]);
        end
        
        J_nondecon(:,:,z-z_quantify_begin(t)+1,t) = J;

        nucmeanI(z-z_quantify_begin(t)+1,t) = mean(mean(I.*uint16(1-(I_4D_mask(:,:,z-z_begin(t)+1,t))))); % Average Mxc level in nuclei, outside HLB
        
        I = imgaussfilt(I, 4); % Blur to get smoother edge
        
        % I_nuc = nuclear mask
        if t>=cc{10}(1) && t<=M{10}(2)
            I_nuc = (I>1.00*nucmeanI(z-z_quantify_begin(t)+1,t));
        elseif t>=cc{11}(1) && t<=M{11}(2)
            I_nuc = (I>1.00*nucmeanI(z-z_quantify_begin(t)+1,t));
        elseif t>=cc{12}(1) && t<=M{12}(2)
            I_nuc = (I>1.00*nucmeanI(z-z_quantify_begin(t)+1,t));
        elseif t>=cc{13}(1) && t<=M{13}(2)
            I_nuc = (I>1.00*nucmeanI(z-z_quantify_begin(t)+1,t));
        elseif t>=cc{14}(1) && t<=M{14}(2)
            I_nuc = (I>0.80*nucmeanI(z-z_quantify_begin(t)+1,t));
        end
%         threshold = 450;
%         if t>=cc{10}(1) && t<=M{10}(2)
%             I_nuc = (I>threshold);
%         elseif t>=cc{11}(1) && t<=M{11}(2)
%             I_nuc = (I>threshold);
%         elseif t>=cc{12}(1) && t<=M{12}(2)
%             I_nuc = (I>threshold);
%         elseif t>=cc{13}(1) && t<=M{13}(2)
%             I_nuc = (I>threshold);
%         elseif t>=cc{14}(1) && t<=M{14}(2)
%             I_nuc = (I>threshold);
%         end
        if exist('emask', 'var'), I_nuc=I_nuc.*emask; end
        
        I_nuc = bwmorph(I_nuc, 'clean');
        
        for num_maj = 1:5
            I_nuc = bwmorph(I_nuc, 'majority');
        end
        
        I_nuc = imfill(I_nuc, 'holes');
        I_nuc = imerode(I_nuc,strel('disk',12,0)); % 15 for zoom of 10, 7 for zoom of 5
        I_nuc = imdilate(I_nuc,strel('disk',12,0)); % 15 for zoom of 10, 7 for zoom of 5
        
        
        erosionnum = 20; % Pick a good number to make sure the nuclear masks to not overlap
        
        
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
        
        I_nucedge=imdilate(I_nuc,strel('disk',1,0))-I_nuc; % To display, not used for calculations
        
        J_4D_mask(:,:,z-z_quantify_begin(t)+1,t) = I_nuc;
        
        area_nucleus = sum(sum(I_nuc));
        
        nucleus_wo_HLB = I_nuc-I_4D_mask(:,:,z-z_begin(t)+1,t);
        mxc_total = sum(sum(nucleus_wo_HLB.*double(I)));
        area_nucleus_wo_HLB(z-z_quantify_begin(t)+1) = sum(sum(nucleus_wo_HLB));
        mxc_nuc_conc_perz(z-z_quantify_begin(t)+1) = mxc_total./area_nucleus_wo_HLB(z-z_quantify_begin(t)+1);
        
        nucleus_whole = I_nuc;
        mxc_total_whole = sum(sum(nucleus_whole.*double(I_nondecon(:,:,z-z_begin(t)+1,t))));
        area_nucleus_whole(z-z_quantify_begin(t)+1) = sum(sum(nucleus_whole));
        mxc_whole_conc_perz(z-z_quantify_begin(t)+1) = mxc_total_whole./area_nucleus_whole(z-z_quantify_begin(t)+1);
        
        % Left = green channel, Right = red channel
        I_3D(:,:,1,z-z_quantify_begin(t)+1,t)= horzcat(50000*uint16(I_nucedge),50000*uint16(I_nucedge)+5*J_nondecon(:,:,z-z_quantify_begin(t)+1,t));
        I_3D(:,:,2,z-z_quantify_begin(t)+1,t)= horzcat(50000*uint16(I_nucedge)+3*I_nondecon(:,:,z-z_begin(t)+1,t),50000*uint16(I_nucedge));
        I_3D(:,:,3,z-z_quantify_begin(t)+1,t)= 0;
    end
    mxc_3D_nuc(t) = sum(mxc_nuc_conc_perz)/sum(area_nucleus_wo_HLB);
    mxc_3D_whole(t) = sum(mxc_whole_conc_perz)/sum(area_nucleus_whole);
end


for n=1:size(I_3D,4)
    if ismember(n,1:zspan), mov{n} = squeeze(I_3D(:,:,:,n,index_begin:index_end)); end
end

for n=1:size(I_3D,4)
    if ismember(n,1:zspan), implay(mov{n}), set(findall(0,'tag','spcui_scope_framework'),'position',[150 150 1620 550]); end
end



tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n',floor(tEnd/60),rem(tEnd,60));

%%  
save(strcat('getWholeMxc_', date, '_series', series), 'mxc_3D_nuc', 'mxc_3D_whole')