%% Use this script along with simpleTrackHLB
clearvars -except master_date master_series master_count isthisseriesconcatenated version
tStart = tic;
load('dateseries')
fprintf('quantifyMxcCdk2HLBpernuc: %s / Series%s \n', date, series)

load(strcat('getHLB_', date, '_series', series))
load(strcat('getNucleusCdk2_', date, '_series', series))

for num=10:14
    if ~isnan(cc{num})
        load(strcat('trackHLB_', date, '_series', series, '_cc', num2str(num)))
        CCt_total(:,:,:,cc{num}(1):M{num}(2)) = CCt(:,:,:,cc{num}(1):M{num}(2));
        V_HLB{num}=cat(2,V_preM,V_postM); F_HLB{num}=cat(2,F_preM,F_postM);
    end
end
clear CCt V_HLB_final V_GFP_final

CCt = CCt_total; clear CCt_total

%% Mark objects (nuclei)
for t=index_begin:index_end
    J_2D_mask = J_3D_centroid(:,:,t);
    DDD = bwlabel(J_2D_mask);
    DD(:,:,t) = DDD;
end
cent_total = zeros(max(max(max(DD))),2,index_end); 
cent_total(cent_total == 0) = NaN;

for t=index_begin:index_end
    J_2D_mask = J_3D_centroid(:,:,t);
    DDD = bwlabel(J_2D_mask);
    DD(:,:,t) = DDD;
    num_obj(t) = max(max(DDD));
    stats=regionprops(DDD, 'Centroid');
    
    for i=1:max(max(DDD))
        if i==1
            cent_initial = stats(1).Centroid;
             cent_one = cent_initial;
        else
            cent = stats(i).Centroid;
            cent_one = vertcat(cent_one, cent);
        end
    end
    
    cent_total(1:size(cent_one,1),:,t) = cent_one;
end

%% Track nuclei
for num=11:14
    if ~isnan(cc{num})
        mindist = zeros(max(num_obj), index_end);
        for t = cc{num}(1):M{num}(2)
            for k = 1:num_obj(t)
                if t==cc{num}(1) % First t
                    DDt(:,:,t) = DD(:,:,t);
                    m_new(k,t) = k;
                    num_label = num_obj(t);
                elseif t==cc{num}(1)+1
                    if num_obj(t) > num_label % If there are more objects
                        num_label = num_obj(t);
                    end
                    cent_curr = cent_total(k,:,t);
                    cent_prev = cent_total(:,:,t-1); % All objects of previous t
                    dist_cent = repmat(cent_curr, size(cent_prev,1),1) - cent_prev;
                    % Calculate distance from kth object at current t to all objects at previous t
                    dist = sqrt(dist_cent(:,1).^2+dist_cent(:,2).^2);
                    new_label = m_new(dist==min(dist),t-1);
                    if isempty(new_label)
                        m_new(k,t) = 0;
                    else
                        m_new(k,t) = new_label(1);
                    end
                    mindist(k,t) = min(dist);
                elseif t==cc{num}(1)+2
                    if num_obj(t) > num_label % If there are more objects
                        num_label = num_obj(t);
                    end
                    cent_curr = cent_total(k,:,t);
                    cent_prev = cent_total(:,:,t-1); % All objects of previous t
                    cent_prev2 = cent_total(:,:,t-2);
                    dist_cent = repmat(cent_curr, size(cent_prev,1),1) - cent_prev;
                    dist_cent2 = repmat(cent_curr, size(cent_prev2,1),1) - cent_prev2;
                    % Calculate distance from kth object at current t to all
                    % objects at previous t AND previous previous t
                    dist = sqrt(dist_cent(:,1).^2+dist_cent(:,2).^2);
                    dist2 = sqrt(dist_cent2(:,1).^2+dist_cent2(:,2).^2);
                    if min(dist2)<0.8*min(dist)
                        new_label = m_new(dist2==min(dist2),t-2);
                    else
                        new_label = m_new(dist==min(dist),t-1);
                    end
                    if isempty(new_label)
                        m_new(k,t) = 0;
                    else
                        m_new(k,t) = new_label(1);
                    end
                    if min(dist2)<0.8*min(dist)
                        mindist(k,t) = min(dist2);
                    else
                        mindist(k,t) = min(dist);
                    end
                else
                    if num_obj(t) > num_label % If there are more objects
                        num_label = num_obj(t);
                    end
                    cent_curr = cent_total(k,:,t);
                    cent_prev = cent_total(:,:,t-1); % All objects of previous t
                    cent_prev2 = cent_total(:,:,t-2);
                    cent_prev3 = cent_total(:,:,t-3);
                    
                    dist_cent = repmat(cent_curr, size(cent_prev,1),1) - cent_prev;
                    dist_cent2 = repmat(cent_curr, size(cent_prev2,1),1) - cent_prev2;
                    dist_cent3 = repmat(cent_curr, size(cent_prev3,1),1) - cent_prev3;
                    % Calculate distance from kth object at current t to all
                    % objects at previous t AND previous previous t AND previous previous previous t
                    dist = sqrt(dist_cent(:,1).^2+dist_cent(:,2).^2);
                    dist2 = sqrt(dist_cent2(:,1).^2+dist_cent2(:,2).^2);
                    dist3 = sqrt(dist_cent3(:,1).^2+dist_cent3(:,2).^2);
                    if min(dist3)<0.9*min(dist2) && min(dist3)<0.9*min(dist)
                        new_label = m_new(dist3==min(dist3),t-3);
                    elseif min(dist2)<0.9*min(dist)
                        new_label = m_new(dist2==min(dist2),t-2);
                    else
                        new_label = m_new(dist==min(dist),t-1);
                    end
                    if isempty(new_label)
                        m_new(k,t) = 0;
                    else
                        m_new(k,t) = new_label(1);
                    end
                    if min(dist3)<0.9*min(dist2) && min(dist3)<0.9*min(dist)
                        mindist(k,t) = min(dist3);
                    elseif min(dist2)<0.9*min(dist)
                        mindist(k,t) = min(dist2);
                    else
                        mindist(k,t) = min(dist);
                    end
                end
                clear cent_prev cent_curr dist_cent dist new_label
            end
            
            for k = 1:num_obj(t)
                if mindist(k,t)<20 && length(find(m_new(:,t)==m_new(k,t)))==1
                    [y_label, x_label] = find(DD(:,:,t)==k);
                    for m = 1:length(x_label)
                        DDt(y_label(m), x_label(m), t) = m_new(k,t); % Replace the target index to new label
                    end
                elseif mindist(k,t)<20 && length(find(m_new(:,t)==m_new(k,t)))>1
                    overlap_index = find(m_new(:,t)==m_new(k,t));
                    mindist_overlap = mindist(overlap_index,t);
                    survivor = find(mindist_overlap == min(mindist_overlap));
                    if k ~= overlap_index(survivor) % If k isn't the closest one, assign a new label
                        num_label = num_label+1;
                        [y_label, x_label] = find(DD(:,:,t)==k);
                        for m = 1:length(x_label)
                            DDt(y_label(m), x_label(m), t) = num_label;
                        end
                        m_new(k,t) = num_label;
                    end
                else
                    num_label = num_label+1;
                    [y_label, x_label] = find(DD(:,:,t)==k);
                    for m = 1:length(x_label)
                        DDt(y_label(m), x_label(m), t) = num_label;
                    end
                    m_new(k,t) = num_label; % Reassign new label
                end
            end
        end
        
        for p=1:max(max(m_new)) % for each object, count how many times it appears
            m_count_{num}(p) = sum(sum(m_new(:,cc{num}(1):cc{num}(2))==p));
        end
        
        m_label_{num} = find(m_count_{num}>0.8*(cc{num}(2)-cc{num}(1)+1)); % This has the object numbers that span over the whole time
        if length(m_label_{num})<5 % If 0.8 threshold is too high so that there aren't enough nuclei, lower the threshold
            m_label_{num} = find(m_count_{num}>0.6*(cc{num}(2)-cc{num}(1)+1)); 
        end
        if length(m_label_{num})<5 
            m_label_{num} = find(m_count_{num}>0.4*(cc{num}(2)-cc{num}(1)+1)); 
        end
        if length(m_label_{num})<5 
            m_label_{num} = find(m_count_{num}>0.2*(cc{num}(2)-cc{num}(1)+1)); 
        end
    end
end

%% Whole Tracking
% To visualize whole elements that are being tracked

for num=10:14
    if ~isnan(cc{num})
        for t = cc{num}(1):M{num}(2)
            DDt_whole(:,:,t) = imdilate(DDt(:,:,t),strel('disk',erosionnum,0)).*ismember(imdilate(DDt(:,:,t),strel('disk',erosionnum,0)), m_label_{num});
            DDt_whole_dilate(:,:,t) = imdilate(DDt_whole(:,:,t),strel('disk',2,0));
            DDt_whole_edge(:,:,t) = DDt_whole_dilate(:,:,t)-DDt_whole(:,:,t);
            CC_whole(:,:,t) = I_4D_mask(:,:,z_middle(t)-z_begin(t)+1,t);
            CC_whole_dilate(:,:,t) = imdilate(CC_whole(:,:,t),strel('disk',2,0));
            CC_whole_edge(:,:,t) = CC_whole_dilate(:,:,t)-CC_whole(:,:,t);

            J_3D_wh(:,:,1,t) = 50000*uint16(DDt_whole_edge(:,:,t));
            J_3D_wh(:,:,2,t) = 3*I_nondecon(:,:,z_middle(t)-z_begin(t)+1,t);
            J_3D_wh(:,:,3,t) = 50000*uint16(CC_whole_edge(:,:,t));
        end
        mov = J_3D_wh(:,:,:,cc{num}(1):M{num}(2));
        mm = implay(mov);
        set(findall(0,'tag','spcui_scope_framework'),'position',[150 150 820 550]);
        set(mm.Parent, 'Name', strcat('Whole Tracking: cc', num2str(num)))
    end
end

%% Individual Tracking
% To visualize tracking for one element

% for num=11:14
%     if ~isnan(cc{num})
%         target = m_label_{num}(1); % Pick a nucleus to track and put in ( )
%         for t = cc{num}(1):M{num}(2)
%             DDt_dilate(:,:,t) = imdilate(DDt(:,:,t)==target,strel('disk',erosionnum+2,0));
%             DDt_edge(:,:,t) = DDt_dilate(:,:,t)-imdilate(DDt(:,:,t)==target,strel('disk',erosionnum,0));
%             
%             J_3D(:,:,1,t) = 50000*uint16(DDt_edge(:,:,t));
%             J_3D(:,:,2,t) = 3*I_nondecon(:,:,z_middle(t)-z_begin(t)+1,t);
%             J_3D(:,:,3,t) = 50000*uint16(CC_whole_edge(:,:,t));
%         end
%         
%         mov2 = J_3D(:,:,:,cc{num}(1):M{num}(2));
%         mm = implay(mov2);
%         set(findall(0,'tag','spcui_scope_framework'),'position',[150 150 820 550]);
%         set(mm.Parent, 'Name', strcat('Individual Tracking: cc', num2str(num)))
%     end
% end

%% Count # of HLBs in nuclei and measure Mxc/Cdk2 per nucleus
for num=10:14
    mxc_nuc_single{num} = NaN; mxc_whole_single{num} = NaN; cdk2_activity_single{num} = NaN;
end
for num=10:14
    if ~isnan(cc{num})
        for t = cc{num}(1):M{num}(2)
            for n=1:length(m_label_{num})
                nucleus_whole_single = imdilate((DDt(:,:,t)==m_label_{num}(n)),strel('disk',erosionnum,0));
                HLB_label = unique(repmat(nucleus_whole_single, [1 1 zspan]).*CCt(:,:,:,t));
                HLB_label(HLB_label==0)=[];
                if ~isempty(HLB_label)
                    for k=1:length(HLB_label)
                        index_HLB(num,t,n,k)=HLB_label(k);
                    end
                else
                    index_HLB(num,t,n,1)=0;
                end
                nucleus_wo_HLB_single = (nucleus_whole_single-I_4D_mask(:,:,z_middle(t)-z_begin(t)+1,t))>0;
                
                % Single nucleoplasmic Mxc
                mxc_total_single = sum(sum(nucleus_wo_HLB_single.*double(I_nondecon(:,:,z_middle(t)-z_begin(t)+1,t))));
                area_nucleus_wo_HLB_single = sum(sum(nucleus_wo_HLB_single));
                mxc_nuc_single{num}(n,t-cc{num}(1)+1) = mxc_total_single./area_nucleus_wo_HLB_single;
                
                % Single whole nuclear Mxc
                mxc_total_whole_single = sum(sum(nucleus_whole_single.*double(I_nondecon(:,:,z_middle(t)-z_begin(t)+1,t))));
                area_nucleus_whole_single = sum(sum(nucleus_whole_single));
                mxc_whole_single{num}(n,t-cc{num}(1)+1) = mxc_total_whole_single./area_nucleus_whole_single;
                
                % Single nuclear Cdk2
                nucleus_whole_single_center = imerode(nucleus_whole_single,strel('disk',12,0)); % Use center part of nucleus for Cdk2 quantification
                nucleus_whole_single_dilate = imdilate(nucleus_whole_single,strel('disk',4,0));
                nucleus_whole_single_dilate2 = imdilate(nucleus_whole_single,strel('disk',7,0));
                nucleus_cytring = nucleus_whole_single_dilate2-nucleus_whole_single_dilate;
                
                area_nucleus_center = sum(sum(nucleus_whole_single_center));
                area_cytoplasm = sum(sum(nucleus_cytring));
                cdk2_total_nucleus = sum(sum(nucleus_whole_single_center.*double(J_nondecon(:,:,t))));
                cdk2_total_cytoplasm = sum(sum(nucleus_cytring.*double(J_nondecon(:,:,t))));
                
                cdk2_nucleus = cdk2_total_nucleus./area_nucleus_center;
                cdk2_cytoplasm = cdk2_total_cytoplasm./area_cytoplasm;
                
                cdk2_activity_single{num}(n,t-cc{num}(1)+1) = cdk2_cytoplasm./cdk2_nucleus;
                
                clear HLB_label nucleus_whole_single nucleus_wo_HLB_single mxc_total_single area_nucleus_wo_HLB_single...
                    mxc_total_whole_single area_nucleus_whole_single nucleus_whole_single_center nucleus_whole_single_dilate...
                    nucleus_whole_single_dilate2 nucleus_cytring area_nucleus_center area_cytoplasm cdk2_total_nucleus cdk2_total_cytoplasm...
                    cdk2_nucleus cdk2_cytoplasm
            end
        end
    end
end

%% Calculate V_HLB/F_HLB per nuclei
for num=10:14
    V_pernuc{num} = NaN; F_pernuc{num} = NaN;
end

for num=10:14
    if ~isnan(cc{num})
        for t = cc{num}(1):M{num}(2)
            for n=1:length(m_label_{num})
                if any(index_HLB(num,t,n,:))
                    for k=1:length(index_HLB(num,t,n,:))
                        if index_HLB(num,t,n,k)~=0
                            V_each(k)=V_HLB{num}(index_HLB(num,t,n,k),t-cc{num}(1)+1);
                            F_each(k)=F_HLB{num}(index_HLB(num,t,n,k),t-cc{num}(1)+1);
                        else
                            V_each(k)=0;
                            F_each(k)=0;
                        end
                    end
                    V_pernuc{num}(n,t-cc{num}(1)+1)=sum(V_each);
                    F_pernuc{num}(n,t-cc{num}(1)+1)=sum(F_each);
                else
                    V_pernuc{num}(n,t-cc{num}(1)+1)=NaN;
                    F_pernuc{num}(n,t-cc{num}(1)+1)=NaN;
                end
                clear V_each F_each
            end
        end
        
        V_pernuc{num}(V_pernuc{num}==0)=NaN; F_pernuc{num}(F_pernuc{num}==0)=NaN;
    end
end

upperbound1 = 6; %1.1*max(max(([V_pernuc{num}])));
upperbound2 = 6e4; %1.1*max(max(([F_pernuc{num}])));

figure; clf;
hold on
for num=10:14
    if ~isnan(cc{num})
        fill([M{num}(1), M{num}(1), M{num}(2), M{num}(2)]*cycle_time/60, [0, upperbound1, upperbound1, 0], [.9 .9 .9], 'linestyle', 'none')
        for t = cc{num}(1):M{num}(2)
            for n=1:length(m_label_{num})
                plot(t*cycle_time/60, V_pernuc{num}(n,t-cc{num}(1)+1), 'k*')
            end
        end
        plot((cc{num}(1):M{num}(2))*cycle_time/60, nanmean(V_pernuc{num}), 'k-')
    end
end
xlabel('Time (min)'), ylabel('V_{HLB}'), title('V_{HLB} per nuclei')
set(gca, 'ylim', [0 upperbound1]), set(gca,'layer','top'), set(gca, 'box', 'off')

figure; clf;
hold on
for num=10:14
    if ~isnan(cc{num})
        fill([M{num}(1), M{num}(1), M{num}(2), M{num}(2)]*cycle_time/60, [0, upperbound2, upperbound2, 0], [.9 .9 .9], 'linestyle', 'none')
        for t = cc{num}(1):M{num}(2)
            for n=1:length(m_label_{num})
                plot(t*cycle_time/60, F_pernuc{num}(n,t-cc{num}(1)+1), 'k*')
            end
        end
        plot((cc{num}(1):M{num}(2))*cycle_time/60, nanmean(F_pernuc{num}), 'k-')
    end
end
xlabel('Time (min)'), ylabel('F_{HLB}'), title('F_{HLB} per nuclei')
set(gca, 'ylim', [0 upperbound2]), set(gca,'layer','top'), set(gca, 'box', 'off')

figure; clf;
hold on
for num=10:14
    if ~isnan(cc{num})
        fill([M{num}(1), M{num}(1), M{num}(2), M{num}(2)]*cycle_time/60, [0, upperbound2, upperbound2, 0], [.9 .9 .9], 'linestyle', 'none')
        for n=1:length(m_label_{num})
            plot((cc{num}(1):M{num}(2))*cycle_time/60, F_pernuc{num}(n,:), 'k-')
        end

        plot((cc{num}(1):M{num}(2))*cycle_time/60, nanmean(F_pernuc{num}), 'r-')
    end
end
xlabel('Time (min)'), ylabel('F_{HLB}'), title('F_{HLB} per nuclei')
set(gca, 'ylim', [0 upperbound2]), set(gca,'layer','top'), set(gca, 'box', 'off')

figure; clf;
hold on
for num=10:14
    if ~isnan(cc{num})
        fill([M{num}(1), M{num}(1), M{num}(2), M{num}(2)]*cycle_time/60, [0, 2.2, 2.2, 0], [.9 .9 .9], 'linestyle', 'none')
        for n=1:length(m_label_{num})
            plot((cc{num}(1):M{num}(2))*cycle_time/60, cdk2_activity_single{num}(n,:), 'k-')
        end

        plot((cc{num}(1):M{num}(2))*cycle_time/60, nanmean(cdk2_activity_single{num}), 'r-')
        plot((cc{num}(1):cc{num}(2))*cycle_time/60, cdk2_preM{num}, 'g-')
    end
end
xlabel('Time (min)'), ylabel('Cdk2 Activity'), title('Cdk2 Activity per nuclei')
set(gca, 'ylim', [0.7 2.5]), set(gca,'layer','top'), set(gca, 'box', 'off')

% for num=10:14
%     if ~isnan(cc{num})
%         figure; clf;
%         hold on
%         for t = cc{num}(2):cc{num}(2)
%             for n=1:length(m_label_{num})
%                 plot(mxc_whole_single{num}(n,t-cc{num}(1)+1), F_pernuc{num}(n,t-cc{num}(1)+1),  'k*')
%             end
%         end
%         xlabel('Mxc\_whole'), ylabel('F_{HLB}'), title(strcat('F_{HLB} vs Mxc\_whole per nucleus, cycle ', num2str(num)))
%         set(gca, 'ylim', [0 upperbound2]), set(gca,'layer','top'), set(gca, 'box', 'off')
%     end
% end

save(strcat('quantifyMxcCdk2HLBpernuc_', date, '_series', series), 'V_pernuc', 'F_pernuc', 'mxc_nuc_single', 'mxc_whole_single', 'cdk2_activity_single')

%%
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n',floor(tEnd/60),rem(tEnd,60));