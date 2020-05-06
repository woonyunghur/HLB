%% HLB mask
clearvars -except master_date master_series master_count isthisseriesconcatenated version
tStart = tic;
current_directory = strsplit(pwd, '\');
date = current_directory{end-1};
series = char(regexp(current_directory{end},'\d+','match'));
save('dateseries', 'date', 'series')
fprintf('getHLB: %s / Series%s \n', date, series)

%% Set up series-specific z-stacks and threshold values
run_fresh = 1; % No retrieving if 1

% If you want to retrieve the saved threshold value for this specific series
if ~run_fresh && exist(strcat('getHLB_', date, '_series', series, '.mat'), 'file')
    load(strcat('getHLB_', date, '_series', series), 'series2', 'series3', 'conc11_12', 'conc12_13', 'pixelxy', 'pixelz', 'cycle_time',...
        'cc', 'start_cycle', 'end_cycle', 'index_begin', 'index_end', 'z_begin', 'z_middle', 'z_end', 'zspan',...
        'threshold_S', 'threshold_M')
        
% or if this is first time running this series or you want to override your values
else

series2 = NaN; series3 = NaN; % Concatenate up to 3 movies
conc11_12 = 0; conc12_13 = 0; % If cycle 11-12 or 12-13 are concatenated

pixelxy = 0.073; %um
pixelz = 1.000;
cycle_time = 16.24;

threshold_S(10) = 6000; threshold_M(10) = threshold_S(10)+1000;
threshold_S(11) = 6000; threshold_M(11) = threshold_S(11)+1000;
threshold_S(12) = 6000; threshold_M(12) = threshold_S(12)+1000;
threshold_S(13) = 6000; threshold_M(13) = threshold_S(13)+1000;
threshold_S(14) = 6000; threshold_M(14) = threshold_S(14)+1000;
threshold_S(threshold_S==0)=NaN; threshold_M(threshold_M==0)=NaN;

%% Define time points manually
for num=10:14, cc{num}=NaN; ana{num}=NaN; M{num}=NaN; end
cc{13}=[1 57]; M{13}=[58 81]; ana{13}=71; 
% cc{10}=[]; M{10}=[]; ana{10}=;
% cc{11}=[]; M{11}=[]; ana{11}=; 
% cc{12}=[]; M{12}=[]; ana{12}=; 
% cc{13}=[]; M{13}=[]; ana{13}=; 
% cc{14}=[]; M{14}=[]; ana{14}=; 

if conc11_12 && ~conc12_13
    cc{12}=cc{12}+conc11_12*M{11}(2); M{12}=M{12}+conc11_12*M{11}(2); ana{12}=ana{12}+conc11_12*M{11}(2); 
    cc{13}=cc{13}+conc11_12*M{11}(2); M{13}=M{13}+conc11_12*M{11}(2); ana{13}=ana{13}+conc11_12*M{11}(2); 
elseif conc12_13 && ~conc11_12
    cc{13}=cc{13}+conc12_13*M{12}(2); M{13}=M{13}+conc12_13*M{12}(2); ana{13}=ana{13}+conc12_13*M{12}(2); 
elseif conc11_12 && conc12_13
    cc{13}=cc{13}+conc11_12*M{11}(2)+conc12_13*M{12}(2); M{13}=M{13}+conc11_12*M{11}(2)+conc12_13*M{12}(2); ana{13}=ana{13}+conc11_12*M{11}(2)+conc12_13*M{12}(2); 
    cc{12}=cc{12}+conc11_12*M{11}(2); M{12}=M{12}+conc11_12*M{11}(2); ana{12}=ana{12}+conc11_12*M{11}(2); 
end

for num=10:14, cycle_check(num)=num*~isnan(cc{num}(1)); end
start_cycle=min(cycle_check(cycle_check>0));
end_cycle=max(cycle_check(cycle_check>0));

index_begin = cc{start_cycle}(1);
index_end = M{end_cycle}(2);

%% Define z stack ranges
% if ~isnan(cc{10}), z_begin(cc{10}(1):M{10}(2)) = 1; z_end(cc{10}(1):M{10}(2)) = 12; end % Define z for each cycle
% if ~isnan(cc{11}), z_begin(cc{11}(1):M{11}(2)) = 1; z_end(cc{11}(1):M{11}(2)) = 12; end
% if ~isnan(cc{12}), z_begin(cc{12}(1):M{12}(2)) = 1; z_end(cc{12}(1):M{12}(2)) = 12; end
if ~isnan(cc{13}), z_begin(cc{13}(1):M{13}(2)) = 1; z_end(cc{13}(1):M{13}(2)) = 12; end 
% if ~isnan(cc{14}), z_begin(cc{14}(1):M{14}(2)) = 1; z_end(cc{14}(1):M{14}(2)) = 12; end 
zspan = 12; % How many z's are covered

% if ~isnan(cc{10}), z_middle(cc{10}(1):M{10}(2)) = 6; end % Choose which z you will use for Cdk2 quantification
% if ~isnan(cc{11}), z_middle(cc{11}(1):M{11}(2)) = 8; end
% if ~isnan(cc{12}), z_middle(cc{12}(1):M{12}(2)) = 6; end
if ~isnan(cc{13}), z_middle(cc{13}(1):M{13}(2)) = 5; end 
% if ~isnan(cc{14}), z_middle(cc{14}(1):M{14}(2)) = 8; end 

end

%% Load images
% Green channel, nondeconvolved
dir_nondecon = strcat('..\..\..\', date, '\Series', series);
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

%% Optional: Select area to exclude, do it only once and save it
% if exist(strcat('emask_', date, '_series', series, '.mat'), 'file')
%     load(strcat('emask_', date, '_series', series), 'emask')
% else
%     figure;
%     t=; % Select which time point
%     if conc11_12 && ~conc12_13
%         if t<=M{11}(2)
%             I=imread([dir_nondecon filesep files_nondecon_comp(t,z_middle(t)).name]);
%         else
%             I=imread([dir_nondecon2 filesep files_nondecon_comp(t,z_middle(t)).name]);
%         end
%     elseif ~conc11_12 && conc12_13
%         if t<=M{12}(2)
%             I=imread([dir_nondecon filesep files_nondecon_comp(t,z_middle(t)).name]);
%         else
%             I=imread([dir_nondecon2 filesep files_nondecon_comp(t,z_middle(t)).name]);
%         end
%     elseif conc11_12 && conc12_13 % when you have 3 separate movies for each cycle
%         if t<=M{11}(2)
%             I=imread([dir_nondecon filesep files_nondecon_comp(t,z_middle(t)).name]);
%         elseif t<=M{12}(2)
%             I=imread([dir_nondecon2 filesep files_nondecon_comp(t,z_middle(t)).name]);
%         else
%             I=imread([dir_nondecon3 filesep files_nondecon_comp(t,z_middle(t)).name]);
%         end
%     else
%         I=imread([dir_nondecon filesep files_nondecon_comp(t,z_middle(t)).name]);
%     end
%     imshow(3*I)
%     h=impoly;
%     emask=1-createmask(h);
%     save(strcat('emask_', date, '_series', series), 'emask')
% end

%% Load images - continued
maxI=NaN(max(z_end), max([cc{:}, M{:}])); meanI=NaN(max(z_end), max([cc{:}, M{:}]));
for t=index_begin:index_end
    for z=z_begin(t):z_end(t)
        if conc11_12 && ~conc12_13
            if t<=M{11}(2)
                I=imread([dir_nondecon filesep files_nondecon_comp(t,z).name]);
            else
                I=imread([dir_nondecon2 filesep files_nondecon_comp(t,z).name]);
            end
        elseif ~conc11_12 && conc12_13
            if t<=M{12}(2)
                I=imread([dir_nondecon filesep files_nondecon_comp(t,z).name]);
            else
                I=imread([dir_nondecon2 filesep files_nondecon_comp(t,z).name]);
            end
        elseif conc11_12 && conc12_13 % when you have 3 separate movies for each cycle
            if t<=M{11}(2)
                I=imread([dir_nondecon filesep files_nondecon_comp(t,z).name]);
            elseif t<=M{12}(2)
                I=imread([dir_nondecon2 filesep files_nondecon_comp(t,z).name]);
            else
                I=imread([dir_nondecon3 filesep files_nondecon_comp(t,z).name]);
            end
        else
            I=imread([dir_nondecon filesep files_nondecon_comp(t,z).name]);
        end
        I_nondecon(:,:,z-z_begin(t)+1,t) = I;
        maxI(z,t) = max(max(I));
        meanI(z,t) = mean(mean(I));
    end
end

for t=index_begin:index_end
    for num=10:14
        if ~isnan(cc{num}), if ismember(t,cc{num}(1):M{num}(2)), minmaxI(t) = min(min(maxI(z_begin(t):z_end(t),cc{num}(1):M{num}(2)))); end, end
        if ~isnan(cc{num}), if ismember(t,cc{num}(1):M{num}(2)), minmeanI(t) = min(min(meanI(z_begin(t):z_end(t),cc{num}(1):M{num}(2)))); end, end
    end
end

%% Create mask
for t=index_begin:index_end
    for z=z_begin(t):z_end(t)
        I=I_nondecon(:,:,z-z_begin(t)+1,t);
        % I = imgaussfilt(I, 0.5); % Blur to get smoother edge
        
        % I2 = HLB mask
        for num=10:14
            % Interphase
            if t>=cc{num}(1) && t<=cc{num}(2)
                I2 = (I>threshold_S(num));
            % Mitosis
            elseif t>=M{num}(1) && t<ana{num}
                I2 = (I>threshold_M(num));
            elseif t>=ana{num} && t<=M{num}(2)
                I2=0;
            end
        end
        
        if exist('emask', 'var')
            I2=I2.*emask;
        end
        
        I2 = bwmorph(I2, 'clean');
        for num_maj = 1:5
            I2 = bwmorph(I2, 'majority');
        end
        I2 = imfill(I2, 'holes');
        
        %             I2 = imerode(I2,strel('disk',1,0));
        %             I2 = imdilate(I2,strel('disk',1,0));
        
        %             [B,L2] = bwboundaries(I2, 'noholes');
        %             stats2 = regionprops(L2,'Area','Centroid');
        %             threshold = 0.75;
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
        %
        L = bwlabel(I2);
        stats = regionprops(L, 'Area','Centroid');
        
        for i = 1:max(max(L))
            area = stats(i).Area;
            [ind1, ind2]= find(L==i);
%             if t>=cc{11}(1) && t<=cc{11}(1)+3
%                 I2(ind1,ind2)=0;
%             end
%             if t>=cc{14}(1) && t<=cc{14}(1)+3
%                 I2(ind1,ind2)=0;
%             end
            if(area<10)
                I2(ind1,ind2)=0;
            end
        end
        
        I3=imdilate(I2,strel('disk',2,0));
        I4=I3-I2;
        
        I_4D_mask(:,:,z-z_begin(t)+1,t) = I2;
        I_mask = 50000*uint16(I4);
        
        I_3D(:,:,1,z-z_begin(t)+1,t)= I_mask;
        I_3D(:,:,2,z-z_begin(t)+1,t)= 3*I_nondecon(:,:,z-z_begin(t)+1,t); %65535/2/max(minmaxI)*I_ori(:,:,z-z_begin(t)+1,t);
        I_3D(:,:,3,z-z_begin(t)+1,t)= 0;
    end
end

for n=1:12
    if ismember(n,1:zspan), mov{n} = squeeze(I_3D(:,:,:,n,index_begin:index_end)); end
end

for n=5:8
    if ismember(n,1:zspan), implay(mov{n}), set(findall(0,'tag','spcui_scope_framework'),'position',[150 150 820 550]); end
end
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n',floor(tEnd/60),rem(tEnd,60));

%%
save(strcat('getHLB_', date, '_series', series), 'I_nondecon', 'I_4D_mask', 'maxI', 'minmaxI',...
    'index_begin', 'index_end', 'cc', 'M', 'ana', 'z_begin', 'z_end', 'zspan', 'z_middle', 'pixelxy', 'pixelz', 'cycle_time',...
    'series', 'series2', 'series3', 'conc11_12', 'conc12_13', 'start_cycle', 'end_cycle', 'threshold_S', 'threshold_M', '-v7.3')
