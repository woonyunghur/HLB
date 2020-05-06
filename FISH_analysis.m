clear; 
tStart = tic;
date = '180715_H3_MPM2';
dir_folders = strcat('..\..\', date); 
files=dir([dir_folders filesep strcat('*Series*')]);

series_number = 13; % Pick a number from 1 to length(files)
% Thrown away: 9
% Bright: 
% No transcription: 
% No HLB but transcription: 

rr=regexp(files(series_number).name,'\d+','match'); % Extracts numbers from the file name
series = rr{1}; 

%% Set up series-specific z-stacks and threshold values
run_fresh = 0; % No retrieving if 1

% If you want to retrieve the saved threshold value for this specific series
if ~run_fresh && exist(strcat('FISH_', date, '_series', series, '.mat'), 'file')
    load(strcat('FISH_', date, '_series', series), 'z_begin', 'z_end', 'zspan', 'z_middle', 'cycle',...
    'strelnum', 'threshold_nuc', 'threshold_green', 'threshold_red', 'z_count_threshold')

% or if this is first time running this series or you want to override your values
else
% z_begin = 1; z_end = 16; zspan = z_end-z_begin+1;
% z_middle = 7; % Pick the best plane for nuclear mask
% 
% cycle = 12;
% 
% strelnum = 10; % 10
% 
% threshold_nuc = 6000; % 4500~7000
% threshold_green = 12000; % 10000~14000
% threshold_red = 14000; % 12000~15000
% 
% z_count_threshold = 2; % 3

end

if ~exist('z_count_threshold', 'var'), z_count_threshold = 3; end

pixelxy = 0.073; %um
pixelz = 0.500;

dir_nondecon = strcat('..\..\', date, '\Series', series); % Load green images

for n=1:length(dir([dir_nondecon filesep strcat('*_ch00*')]))
    files_green(:,n)= dir([dir_nondecon filesep strcat('*z', sprintf('%02d', n-1), '_ch00*')]);
    files_red(:,n)= dir([dir_nondecon filesep strcat('*z', sprintf('%02d', n-1), '_ch01*')]);
end

for z=z_begin:z_end
    I=imread([dir_nondecon filesep files_green(z).name]);
    I_nondecon(:,:,z-z_begin+1) = I;
    J=imread([dir_nondecon filesep files_red(z).name]);
    J_nondecon(:,:,z-z_begin+1) = J;
end

%% Nuclear mask

for z=z_begin:z_end
    % Red channel
    J=J_nondecon(:,:,z_middle-z_begin+1);
    J = imgaussfilt(J, 2); % Blur to get smoother edge
    J2 = (J>threshold_nuc);
    
    J2 = bwmorph(J2, 'clean');
    for num_maj = 1:5
        J2 = bwmorph(J2, 'majority');
    end
    J2 = imfill(J2, 'holes');
    
    J2 = imerode(J2,strel('disk',strelnum,0)); % 15 for zoom of 10, 7 for zoom of 5
    J2 = imdilate(J2,strel('disk',strelnum,0)); % 15 for zoom of 10, 7 for zoom of 5

    L = bwlabel(J2);
    stats = regionprops(L, 'Area','Centroid');
    
    for i = 1:max(max(L))
        area = stats(i).Area;
        [ind1, ind2]= find(L==i);
        if(area<200)
            J2(ind1,ind2)=0;
        end
    end

    J3=imdilate(J2,strel('disk',2,0));
    J4=J3-J2;

    J_nuclear_mask = J2;
    J_nuc = 50000*uint16(J4);

    J_nuc_3D(:,:,1,z-z_begin+1)= 2*J_nondecon(:,:,z-z_begin+1);
    J_nuc_3D(:,:,2,z-z_begin+1)= J_nuc; 
    J_nuc_3D(:,:,3,z-z_begin+1)= 0;
end

mov3 = squeeze(J_nuc_3D);
implay(mov3), set(findall(0,'tag','spcui_scope_framework'),'position',[150 150 820 550]);

% Calculate the N/C ratio
area_nucleus = sum(sum(J_nuclear_mask));

istherenucleus_dim1=any(J_nuclear_mask,2);
dim1_begin = 1; dim1_end = size(J_nuclear_mask,1);
while istherenucleus_dim1(dim1_begin)==0, dim1_begin = dim1_begin+1; end
while istherenucleus_dim1(dim1_end)==0, dim1_end = dim1_end-1; end

istherenucleus_dim2=any(J_nuclear_mask,1);
dim2_begin = 1; dim2_end = size(J_nuclear_mask,2);
while istherenucleus_dim2(dim2_begin)==0, dim2_begin = dim2_begin+1; end
while istherenucleus_dim2(dim2_end)==0, dim2_end = dim2_end-2; end

area_cytoplasm = ((dim1_end-dim1_begin+1)*(dim2_end-dim2_begin+1))-area_nucleus;
ncratio = area_nucleus/area_cytoplasm;
 

%% HLB mask for green channel - raw
for z=z_begin:z_end
    % Green channel
    I=I_nondecon(:,:,z-z_begin+1);
    I = imgaussfilt(I, 0.5); % Blur to get smoother edge
    I2 = (I>threshold_green).*J_nuclear_mask;
    
    I2 = bwmorph(I2, 'clean');
    for num_maj = 1:5
        I2 = bwmorph(I2, 'majority');
    end
    I2 = imfill(I2, 'holes');
    
    %     I2 = imerode(I2,strel('disk',1,0));
    %     I2 = imdilate(I2,strel('disk',1,0));
    %
    %     [B,L2] = bwboundaries(I2, 'noholes');
    %     stats2 = regionprops(L2,'Area','Centroid');
    %     threshold = 0.75;
    %     for k=1:length(B)
    %         boundary = B{k};
    %         delta_sq = diff(boundary).^2;
    %         perimeter = sum(sqrt(sum(delta_sq,2)));
    %         area2 = stats2(k).Area;
    %         metric = 4*pi*area2/perimeter^2;
    %         [ind3, ind4] = find(L2==k);
    %         if metric<threshold
    %             I2(ind3,ind4)=0;
    %         end
    %     end
    %
    L = bwlabel(I2);
    stats = regionprops(L, 'Area','Centroid');
    
    for i = 1:max(max(L))
        area = stats(i).Area;
        [ind1, ind2]= find(L==i);
        if(area<10)
            I2(ind1,ind2)=0;
        end
    end
    
    I3=imdilate(I2,strel('disk',2,0));
    I4=I3-I2;
    
    I_3D_mask(:,:,z-z_begin+1) = I2;
    I_mask = 50000*uint16(I4);
    
    I_3D(:,:,1,z-z_begin+1)= I_mask;
    I_3D(:,:,2,z-z_begin+1)= I_nondecon(:,:,z-z_begin+1);
    I_3D(:,:,3,z-z_begin+1)= 0;
end

mov = squeeze(I_3D);
implay(mov), set(findall(0,'tag','spcui_scope_framework'),'position',[150 150 820 550]);

I_HLB_any = zeros(500,800);
for z=z_begin:z_end
    I_HLB_any=I_HLB_any+I_3D_mask(:,:,z-z_begin+1);
end
I_HLB_any = (I_HLB_any>0);

%% HLB mask for red channel
for z=z_begin:z_end
    % Red channel
    J=J_nondecon(:,:,z-z_begin+1);
    J = imgaussfilt(J, 0.5); % Blur to get smoother edge
    J2 = (J>threshold_red).*imdilate(I_HLB_any,strel('disk',5,0)); % Do it only when there is a transcription signal in any z at that location
    
    J2 = bwmorph(J2, 'clean');
    for num_maj = 1:5
        J2 = bwmorph(J2, 'majority');
    end
    J2 = imfill(J2, 'holes');
    
    L = bwlabel(J2);
    stats = regionprops(L, 'Area','Centroid');
    
    for i = 1:max(max(L))
        area = stats(i).Area;
        [ind1, ind2]= find(L==i);
        if(area>300)
            J2(ind1,ind2)=0;
        end
    end
    
    J3=imdilate(J2,strel('disk',2,0));
    J4=J3-J2;
    
    J_3D_mask(:,:,z-z_begin+1) = J2;
    J_mask = 50000*uint16(J4);
    
    J_3D(:,:,1,z-z_begin+1)= 2*J_nondecon(:,:,z-z_begin+1);
    J_3D(:,:,2,z-z_begin+1)= J_mask;
    J_3D(:,:,3,z-z_begin+1)= 0;
end

mov2 = squeeze(J_3D);
implay(mov2), set(findall(0,'tag','spcui_scope_framework'),'position',[150 150 820 550]);

%% Refine green mask based on red mask
for z=z_begin:z_end
    % Green channel
    I2 = I_3D_mask(:,:,z-z_begin+1);
    if z==z_begin
        J_3D_span = J_3D_mask(:,:,z-z_begin+1)+J_3D_mask(:,:,z-z_begin+2);
    elseif z==z_end
        J_3D_span = J_3D_mask(:,:,z-z_begin+1)+J_3D_mask(:,:,z-z_begin);
    else
        J_3D_span = J_3D_mask(:,:,z-z_begin+1)+J_3D_mask(:,:,z-z_begin)+J_3D_mask(:,:,z-z_begin+2);
    end
    J_3D_span = (J_3D_span>0);
    
    I2 = I2.*imdilate(J_3D_span,strel('disk',2,0));
    
    I3=imdilate(I2,strel('disk',2,0));
    I4=I3-I2;
    
    I_3D_mask(:,:,z-z_begin+1) = I2;
    I_mask = 50000*uint16(I4);
    
    I_3D_re(:,:,1,z-z_begin+1)= I_mask;
    I_3D_re(:,:,2,z-z_begin+1)= I_nondecon(:,:,z-z_begin+1);
    I_3D_re(:,:,3,z-z_begin+1)= 0;
end

mov4 = squeeze(I_3D_re);
implay(mov4), set(findall(0,'tag','spcui_scope_framework'),'position',[150 150 820 550]);

%% Refine red mask based on green mask 
for z=z_begin:z_end
    % Red channel
    J2 = J_3D_mask(:,:,z-z_begin+1);
    if z==z_begin
        I_3D_span = I_3D_mask(:,:,z-z_begin+1)+I_3D_mask(:,:,z-z_begin+2);
    elseif z==z_end
        I_3D_span = I_3D_mask(:,:,z-z_begin+1)+I_3D_mask(:,:,z-z_begin);
    else
        I_3D_span = I_3D_mask(:,:,z-z_begin+1)+I_3D_mask(:,:,z-z_begin)+I_3D_mask(:,:,z-z_begin+2);
    end
    I_3D_span = (I_3D_span>0);
    
    J2 = J2.*imdilate(I_3D_span,strel('disk',5,0));
    
    J3=imdilate(J2,strel('disk',2,0));
    J4=J3-J2;
    
    J_3D_mask(:,:,z-z_begin+1) = J2;
    J_mask = 50000*uint16(J4);
    
    J_3D_re(:,:,1,z-z_begin+1)= 2*J_nondecon(:,:,z-z_begin+1);
    J_3D_re(:,:,2,z-z_begin+1)= J_mask;
    J_3D_re(:,:,3,z-z_begin+1)= 0;
end    

mov5 = squeeze(J_3D_re);
implay(mov5), set(findall(0,'tag','spcui_scope_framework'),'position',[150 150 820 550]);

%% Calculate noise
% Calculate the noise level based on the fluorescence
noise_nuc=(J_nuclear_mask-I_3D_mask(:,:,z_middle-z_begin+1)).*double(I_nondecon(:,:,z_middle-z_begin+1)); % Nucleus only
noise_cyt=(1-J_nuclear_mask).*double(I_nondecon(:,:,z_middle-z_begin+1)); % Cytoplasm only
noise_nuc=reshape(noise_nuc,[],1); noise_cyt=reshape(noise_cyt,[],1);
noise_nuc(noise_nuc==0)=[]; noise_cyt(noise_cyt==0)=[];
noise_level_nuc = mean(noise_nuc); noise_level_cyt = mean(noise_cyt);


%% H3 transcription
% Green channel
CC_green = bwlabeln(I_3D_mask);
num_obj_CC_green = max(max(max(CC_green)));

for i=1:num_obj_CC_green
    clear xx yy zz
    clear CCt_logical CCt_logical_dilated CCt_logical_edge
    clear area V_trapcone
    z_count = 0;
    
    for z=z_begin:z_end
        CC_green_logical(:,:,z-z_begin+1)=(CC_green(:,:,z-z_begin+1)==i);
    end
    
    for z=z_begin:z_end
        area(z-z_begin+1) = sum(sum(CC_green_logical(:,:,z-z_begin+1)));
        area_GFP(z-z_begin+1) = sum(sum(double(CC_green_logical(:,:,z-z_begin+1)).*double(I_nondecon(:,:,z-z_begin+1))));
        if area(z-z_begin+1)>0
            z_count = z_count+1;
        end
    end
    if z_count>=z_count_threshold
        for z=z_begin:z_end-1
            V_trapcone(z-z_begin+1) = (area(z-z_begin+1)+area(z-z_begin+1+1)+sqrt(area(z-z_begin+1)*area(z-z_begin+1+1)))/3*(pixelxy^2)*pixelz;
        end
        V_green_HLB(i) = sum(V_trapcone);
        F_green_HLB(i) = sum(area_GFP)/sum(area)*V_green_HLB(i);
    else
        V_green_HLB(i) = NaN;
        F_green_HLB(i) = NaN;
    end
end

%% MPM-2
% Red channel
CC_red = bwlabeln(J_3D_mask);
num_obj_CC_red = max(max(max(CC_red)));

for i=1:num_obj_CC_red
    clear xx yy zz
    clear CCt_logical CCt_logical_dilated CCt_logical_edge
    clear area V_trapcone
    z_count = 0;
    
    for z=z_begin:z_end
        CC_red_logical(:,:,z-z_begin+1)=(CC_red(:,:,z-z_begin+1)==i);
    end
    
    for z=z_begin:z_end
        area(z-z_begin+1) = sum(sum(CC_red_logical(:,:,z-z_begin+1)));
        area_GFP(z-z_begin+1) = sum(sum(double(CC_red_logical(:,:,z-z_begin+1)).*double(I_nondecon(:,:,z-z_begin+1))));
        if area(z-z_begin+1)>0
            z_count = z_count+1;
        end
    end
    if z_count>=z_count_threshold
        for z=z_begin:z_end-1
            V_trapcone(z-z_begin+1) = (area(z-z_begin+1)+area(z-z_begin+1+1)+sqrt(area(z-z_begin+1)*area(z-z_begin+1+1)))/3*(pixelxy^2)*pixelz;
        end
        V_red_HLB(i) = sum(V_trapcone);
        F_red_HLB(i) = sum(area_GFP)/sum(area)*V_red_HLB(i);
    else
        V_red_HLB(i) = NaN;
        F_red_HLB(i) = NaN;
    end
end

%% Quantify single HLB

for i=1:num_obj_CC_green
    red_index_green_raw = unique((CC_green==i).*CC_red);
    red_index_green_raw(red_index_green_raw==0) = [];
    if ~isempty(red_index_green_raw) && length(red_index_green_raw)==1
        red_index_green(i) = red_index_green_raw(1);
        V_red_HLB_single(i) = V_red_HLB(red_index_green_raw(1));
        F_red_HLB_single(i) = F_red_HLB(red_index_green_raw(1));
    else
        red_index_green(i) = NaN;
        V_red_HLB_single(i) = NaN;
        F_red_HLB_single(i) = NaN;
    end
end

figure('Position', [0 500 560 420]); clf; clf
hold on
for i=1:num_obj_CC_green
    plot(F_red_HLB_single(i), F_green_HLB(i), 'linestyle', 'none', 'color', 'k', 'marker', '.', 'markersize', 15)
end
xlabel('F_{HLB} (\mum^3)'), ylabel('Transcription (A.U.)')

F_green_HLB_single = F_green_HLB;

%% Quantify HLB per nuclei
DD = bwlabel(J_nuclear_mask);
num_nuc = max(max(DD)); % Label nuclei

% Green channel
for n=1:num_nuc
    HLB_green_label = unique(repmat(imdilate((DD==n),strel('disk',10,0)), [1 1 zspan]).*CC_green);
    HLB_green_label(HLB_green_label==0)=[];
    if ~isempty(HLB_green_label)
        for k=1:length(HLB_green_label)
            index_HLB_green(n,k)=HLB_green_label(k);
        end
    else
        index_HLB_green(n,1)=0;
    end
    clear HLB_green_label
end

for n=1:num_nuc
    countnum = 0;
    if any(index_HLB_green(n,:))
        for k=1:length(index_HLB_green(n,:))
            if index_HLB_green(n,k)~=0
                V_green_each(k)=V_green_HLB(index_HLB_green(n,k));
                F_green_each(k)=F_green_HLB(index_HLB_green(n,k));
                countnum=countnum+1;
            else
                V_green_each(k)=0;
                F_green_each(k)=0;
            end
            V_green_pernuc(n)=sum(V_green_each);
            F_green_pernuc(n)=sum(F_green_each);
        end
    else
        V_green_pernuc(n)=0;
        F_green_pernuc(n)=0;
    end
    countHLBpernuc(n) = countnum;
    clear V_green_each F_green_each
end
V_green_pernuc(V_green_pernuc==0)=NaN;
F_green_pernuc(F_green_pernuc==0)=NaN;

% Red channel
for n=1:num_nuc
    HLB_red_label = unique(repmat((DD==n), [1 1 zspan]).*CC_red);
    HLB_red_label(HLB_red_label==0)=[];
    if ~isempty(HLB_red_label)
        for k=1:length(HLB_red_label)
            index_HLB_red(n,k)=HLB_red_label(k);
        end
    else
        index_HLB_red(n,1)=0;
    end
    clear HLB_red_label
end

for n=1:num_nuc
    if any(index_HLB_red(n,:))
        for k=1:length(index_HLB_red(n,:))
            if index_HLB_red(n,k)~=0
                V_red_each(k)=V_red_HLB(index_HLB_red(n,k));
                F_red_each(k)=F_red_HLB(index_HLB_red(n,k));
            else
                V_red_each(k)=0;
                F_red_each(k)=0;
            end
            V_red_pernuc(n)=sum(V_red_each);
            F_red_pernuc(n)=sum(F_red_each);
        end
    else
        V_red_pernuc(n)=0;
        F_red_pernuc(n)=0;
    end
    clear V_red_each F_red_each
end
V_red_pernuc(V_red_pernuc==0)=NaN;
F_red_pernuc(F_red_pernuc==0)=NaN;

F_green_pernuc_nanremoved = F_green_pernuc;
F_green_pernuc_nanremoved(isnan(F_green_pernuc_nanremoved))=[];
transcription = trimmean(F_green_pernuc_nanremoved,25);
F_red_pernuc_nanremoved = F_red_pernuc;
F_red_pernuc_nanremoved(isnan(F_red_pernuc_nanremoved))=[];
F_HLB = trimmean(F_red_pernuc_nanremoved,25);
V_red_pernuc_nanremoved = V_red_pernuc;
V_red_pernuc_nanremoved(isnan(V_red_pernuc_nanremoved))=[];
V_HLB = trimmean(V_red_pernuc_nanremoved,25);

%% Sort single HLB to paired/unpaired
for i=1:num_obj_CC_green
    whichnuc = unique((CC_green==i).*repmat(imdilate(DD,strel('disk',5,0)), [1 1 zspan]));
    whichnuc(whichnuc==0) = [];
    if countHLBpernuc(whichnuc)==1
        F_green_HLB_single_paired(i) = F_green_HLB_single(i);
        F_green_HLB_single_unpaired(i) = NaN;
        F_red_HLB_single_paired(i) = F_red_HLB_single(i);
        F_red_HLB_single_unpaired(i) = NaN;
        V_red_HLB_single_paired(i) = V_red_HLB_single(i);
        V_red_HLB_single_unpaired(i) = NaN;
    elseif countHLBpernuc(whichnuc)==2
        F_green_HLB_single_paired(i) = NaN;
        F_green_HLB_single_unpaired(i) = F_green_HLB_single(i);
        F_red_HLB_single_paired(i) = NaN;
        F_red_HLB_single_unpaired(i) = F_red_HLB_single(i);
        V_red_HLB_single_paired(i) = NaN;
        V_red_HLB_single_unpaired(i) = V_red_HLB_single(i);
    end
end
F_green_HLB_single_paired(F_green_HLB_single_paired==0)=NaN;
F_green_HLB_single_unpaired(F_green_HLB_single_unpaired==0)=NaN;
F_red_HLB_single_paired(F_red_HLB_single_paired==0)=NaN;
F_red_HLB_single_unpaired(F_red_HLB_single_unpaired==0)=NaN;
V_red_HLB_single_paired(V_red_HLB_single_paired==0)=NaN;
V_red_HLB_single_unpaired(V_red_HLB_single_unpaired==0)=NaN;

%% Plot and save

figure('Position', [600 500 560 420]); clf
plot(F_red_pernuc, F_green_pernuc, 'linestyle', 'none', 'color', 'k', 'marker', '.', 'markersize', 15)
xlabel('F_{HLB}/nucleus (\mum^3)'), ylabel('Transcription (A.U.)')

figure('Position', [1200 500 560 420]); clf
hold on
plot(F_red_HLB_single_paired, F_green_HLB_single_paired, 'linestyle', 'none', 'color', 'r', 'marker', '.', 'markersize', 15)
plot(F_red_HLB_single_unpaired, F_green_HLB_single_unpaired, 'linestyle', 'none', 'color', 'b', 'marker', '.', 'markersize', 15)
xlabel('F_{HLB} (\mum^3)'), ylabel('Transcription (A.U.)')

tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n',floor(tEnd/60),rem(tEnd,60));

save(strcat('FISH_', date, '_series', series), 'z_begin', 'z_end', 'zspan', 'z_middle', 'cycle',...
    'strelnum', 'threshold_nuc', 'threshold_green', 'threshold_red', 'z_count_threshold', 'ncratio',...
    'F_green_HLB_single', 'V_red_HLB_single', 'F_red_HLB_single', ...
    'F_green_HLB_single_paired', 'F_green_HLB_single_unpaired', 'F_red_HLB_single_paired', 'F_red_HLB_single_unpaired',...
    'V_red_HLB_single_paired', 'V_red_HLB_single_unpaired',...
    'transcription', 'F_HLB', 'V_HLB', 'F_green_pernuc', 'F_red_pernuc', 'V_red_pernuc',...
    'noise_level_nuc', 'noise_level_cyt')