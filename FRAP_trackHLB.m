%% Calculate fluorescence for bleached/non-bleached HLBs
clear; tStart = tic;
date = '170317_FRAP_2';
load(strcat('getHLB_', date))

%% Mark objects (HLB)
for m=1:index_length
    for t=1:pb_length+1
        I_2D_mask = I_mask{m}(:,:,t);
        CCC = bwlabel(I_2D_mask);
        CC{m}(:,:,t) = CCC;
    end
    cent_total = zeros(max(max(max(CC{m}))),2,pb_length+1);
    cent_total(cent_total == 0) = NaN;
    
    for t=1:pb_length+1
        I_2D_mask = I_mask{m}(:,:,t);
        CCC = bwlabel(I_2D_mask);
        CC{m}(:,:,t) = CCC;
        num_obj(t) = max(max(CCC));
        stats=regionprops(CCC, 'Centroid');
        
        for i=1:max(max(CCC))
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
    
    % Find center of bleaching point
    CCC_bleach = bwlabel(I_bleach_mask{m});
    stats=regionprops(CCC_bleach, 'Centroid');
    cent_bleach{m} = stats.Centroid;
    
    %% Track HLB
    mindist = zeros(max(num_obj), pb_length+1);
    for t = 1:pb_length+1
        for k = 1:num_obj(t)
            if t==1 % First t
                CCt{m}(:,:,t) = CC{m}(:,:,t);
                m_new(k,t) = k;
                num_label = num_obj(t);
                cent_initall = cent_total(:,:,t);
                dist_bleach_cent = repmat(cent_bleach{m}, size(cent_initall,1),1) - cent_initall;
                dist_bleach = sqrt(dist_bleach_cent(:,1).^2+dist_bleach_cent(:,2).^2);
                bleach_index=find(dist_bleach==min(dist_bleach));
            elseif t==2
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
            elseif t==3
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
            clear cent_curr dist new_label cent_prev cent_prev2 cent_prev3 dist_cent dist_cent2 dist_cent3
        end
        
        for k = 1:num_obj(t)
            if mindist(k,t)<20 && length(find(m_new(:,t)==m_new(k,t)))==1
                [y_label, x_label] = find(CC{m}(:,:,t)==k);
                for q = 1:length(x_label)
                    CCt{m}(y_label(q), x_label(q), t) = m_new(k,t); % Replace the target index to new label
                end
            elseif mindist(k,t)<20 && length(find(m_new(:,t)==m_new(k,t)))>1
                overlap_index = find(m_new(:,t)==m_new(k,t));
                mindist_overlap = mindist(overlap_index,t);
                survivor = find(mindist_overlap == min(mindist_overlap));
                if k ~= overlap_index(survivor) % If k isn't the closest one, assign a new label
                    num_label = num_label+1;
                    [y_label, x_label] = find(CC{m}(:,:,t)==k);
                    for q = 1:length(x_label)
                        CCt{m}(y_label(q), x_label(q), t) = num_label;
                    end
                    m_new(k,t) = num_label;
                end
            else
                num_label = num_label+1;
                [y_label, x_label] = find(CC{m}(:,:,t)==k);
                for q = 1:length(x_label)
                    CCt{m}(y_label(q), x_label(q), t) = num_label;
                end
                m_new(k,t) = num_label; % Reassign new label
            end
        end
    end
    clear mindist num_label new_label overlap_index mindist_overlap survivor  
    
    for p=1:max(max(m_new)) % for each object, count how many times it appears
        m_count(p) = sum(sum(m_new(:,1:pb_length+1)==p));
    end
    
    clear m_new
    m_label = find(m_count>(pb_length-2)); % This has the object numbers that span over the whole time
    clear m_count
    
    %% Whole Tracking
    % To visualize whole elements that are being tracked
    for t = 1:pb_length+1
        CCt_whole(:,:,t) = CCt{m}(:,:,t).*ismember(CCt{m}(:,:,t), m_label);
    end
%     for t = 1:pb_length+1
%         CCt_whole_dilate(:,:,t) = imdilate(CCt_whole(:,:,t),strel('disk',2,0));
%         CCt_whole_edge(:,:,t) = CCt_whole_dilate(:,:,t)-CCt_whole(:,:,t);
%         
%         if t==1
%             I_3D_wh(:,:,1,t) = 50000*uint16(CCt_whole_edge(:,:,t));
%             I_3D_wh(:,:,2,t) = 2*I_comp{m}(:,:,t);
%             I_3D_wh(:,:,3,t) = 0;
%         else
%             I_3D_wh(:,:,1,t+1) = 50000*uint16(CCt_whole_edge(:,:,t));
%             I_3D_wh(:,:,2,t+1) = 2*I_comp{m}(:,:,t);
%             I_3D_wh(:,:,3,t+1) = 0;
%         end
%     end
%     I_3D_wh(:,:,1,2)= 0;
%     I_3D_wh(:,:,2,2)= 2*I_bleach{m};
%     I_3D_wh(:,:,3,2)= 0;
% 
%     mov = I_3D_wh(:,:,:,1:pb_length+2);
%     mm = implay(mov);
%     set(findall(0,'tag','spcui_scope_framework'),'position',[150 150 820 550]);
%     set(mm.Parent, 'Name', strcat('Whole Tracking: cc', num2str(mcc(m)), ' Set #', num2str(m)))
    
    %% Individual Tracking
    % To visualize tracking for one element
%     target = m_label(1); % Pick a number
%     for t = 1:pb_length+1
%         CCt_dilate(:,:,t) = imdilate(CCt{m}(:,:,t)==target,strel('disk',2,0));
%         CCt_edge(:,:,t) = CCt_dilate(:,:,t)-(CCt{m}(:,:,t)==target);
%     
%         if t==1
%         I_3D_ind(:,:,1,t) = 50000*uint16(CCt_edge(:,:,t));
%         I_3D_ind(:,:,2,t) = 2*I_comp{m}(:,:,t);
%         I_3D_ind(:,:,3,t) = 0;
%         else
%         I_3D_ind(:,:,1,t+1) = 50000*uint16(CCt_edge(:,:,t));
%         I_3D_ind(:,:,2,t+1) = 2*I_comp{m}(:,:,t);
%         I_3D_ind(:,:,3,t+1) = 0;
%         end
%     end
%     I_3D_ind(:,:,1,2)= 0;
%     I_3D_ind(:,:,2,2)= 2*I_bleach{m};
%     I_3D_ind(:,:,3,2)= 0;
%     
%     mov2 = I_3D_ind(:,:,:,1:pb_length+2);
%     mm = implay(mov2);
%     set(findall(0,'tag','spcui_scope_framework'),'position',[150 150 820 550]);
%     set(mm.Parent, 'Name', strcat('Individual Tracking: cc', num2str(mcc(m)), ' Set #', num2str(m)))
    
    %% Bleached HLB Tracking
    
    for t = 1:pb_length+1
        CCt_bleach(:,:,t) = CCt{m}(:,:,t).*double(CCt{m}(:,:,t)==bleach_index);
        CCt_bleach_dilate(:,:,t) = imdilate(CCt_bleach(:,:,t),strel('disk',2,0));
        CCt_bleach_edge(:,:,t) = CCt_bleach_dilate(:,:,t)-(CCt_bleach(:,:,t));
        
        CCt_nonbleach(:,:,t) = CCt_whole(:,:,t)-CCt_bleach(:,:,t);
        CCt_nonbleach_dilate(:,:,t) = imdilate(CCt_nonbleach(:,:,t),strel('disk',2,0));
        CCt_nonbleach_edge(:,:,t) = CCt_nonbleach_dilate(:,:,t)-(CCt_nonbleach(:,:,t));
        
        if t==1
            I_3D_bl(:,:,1,t) = 50000*uint16(CCt_bleach_edge(:,:,t));
            I_3D_bl(:,:,2,t) = 2*I_comp{m}(:,:,t);
            I_3D_bl(:,:,3,t) = 50000*uint16(CCt_nonbleach_edge(:,:,t));
        else
            I_3D_bl(:,:,1,t+1) = 50000*uint16(CCt_bleach_edge(:,:,t));
            I_3D_bl(:,:,2,t+1) = 2*I_comp{m}(:,:,t);
            I_3D_bl(:,:,3,t+1) = 50000*uint16(CCt_nonbleach_edge(:,:,t));
        end
    end
    I_3D_bl(:,:,1,2)= 0;
    I_3D_bl(:,:,2,2)= 2*I_bleach{m};
    I_3D_bl(:,:,3,2)= 0;
    
%     mov3 = I_3D_bl(:,:,:,1:pb_length+2);
%     mm = implay(mov3);
%     set(findall(0,'tag','spcui_scope_framework'),'position',[150 150 820 550]);
%     set(mm.Parent, 'Name', strcat('Bleached HLB Tracking: cc', num2str(mcc(m)), ' Set #', num2str(m)))
    
    %% Calculate and record fluorescence
    
    for t = 1:pb_length+1
        F_bleached{m}(t) = sum(sum(double(I_comp{m}(:,:,t)).*(CCt_bleach(:,:,t)>0))) / sum(sum(CCt_bleach(:,:,t)>0)); % total fluorescence / area
        F_rest{m}(t) = sum(sum(double(I_comp{m}(:,:,t)).*(CCt_nonbleach(:,:,t)>0))) / sum(sum(CCt_nonbleach(:,:,t)>0));
        if t>3 && isnan(F_bleached{m}(t-1)) % connect missing points
            F_bleached{m}(t-1)= (F_bleached{m}(t-2)+F_bleached{m}(t))/2;
        end
    end
    
    % Smoothen the bleached data
    F_bleached_smooth{m}=sgolayfilt(F_bleached{m}(2:end),3,5);
    
    % Nonbleached & bleached
    x{m}=(1:1:pb_length)*cycle_time;
    y1=(F_rest{m}(2:end)-F_bleached_smooth{m})./F_rest{m}(2:end); y1(y1<=0)=NaN; y{m}=log(y1);
    x{m}(isnan(y{m}))=[]; y{m}(isnan(y{m}))=[];
    [f,gof]=fit(x{m}',y{m}','poly1');
    R2_pre{m} = gof.rsquare;
    koff_pre{m} = -f.p1;
    const_pre{m} = f.p2;
    
    diffxy{m} = abs(f.p1*x{m}+f.p2-y{m});
    index_extreme{m}= find(diffxy{m}>=2*std(diffxy{m}));
    index_normal{m} = find(diffxy{m}<2*std(diffxy{m}));
    
    [f2,gof2]=fit((x{m}(index_normal{m}))',(y{m}(index_normal{m}))','poly1');
    R2{m} = gof2.rsquare;
    koff{m} = -f2.p1;
    const{m} = f2.p2;
    
%     if R2{m}>= 0.6
        fig = figure; fig.Name = strcat('#', num2str(m), ' cc', num2str(mcc(m))); set(gcf, 'position', [10, 400, 1000, 400]);
        subplot(1,2,1), hold on
        plot((2:1:pb_length+1)*cycle_time, F_bleached{m}(2:end), 'color', [1 0.8 0.8]);
        pk=plot((1:1:pb_length+1)*cycle_time, F_rest{m}, 'k-');
        pr=plot((2:1:pb_length+1)*cycle_time, F_bleached_smooth{m}, 'r-');
        set(findobj(gca, 'Type', 'Line'), 'LineWidth', 2)
        set(gca, 'ylim', [0 inf]), set(gca, 'xlim', [0 pb_length*cycle_time])
        set(gca,'layer','top'), set(gca, 'box', 'off'), set(gca, 'fontweight','bold','fontsize', 12)
        legend([pr pk], {'Bleached', 'Non-bleached'}, 'location', 'southeast')
        xlabel('Time (sec)')
        
        if ~isempty(y)
            subplot(1,2,2), hold on
            plot(x{m}(index_normal{m}), y{m}(index_normal{m}), 'linestyle', 'none', 'marker', '*', 'color', [0 0 1])
            plot(x{m}(index_extreme{m}), y{m}(index_extreme{m}), 'linestyle', 'none', 'marker', '*', 'color', [0.7 0.7 1])
            plot(x{m}, f.p1*x{m}+f.p2, 'color', [0.7 1 0.7])
            plot(x{m}(index_normal{m}), f2.p1*x{m}(index_normal{m})+f2.p2, 'color', [0.2 1 0.2])
            set(gca, 'xlim', [0 pb_length*cycle_time])
            xlim = get(gca, 'xlim'); ylim = get(gca, 'ylim');
            text(xlim(2),ylim(2), strcat('R^2=', num2str(R2{m}, '%.4f')), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'fontsize', 12, 'fontweight', 'bold')
            text(xlim(2),ylim(2)-(ylim(2)-ylim(1))*0.1, strcat('k_{off}=', num2str(koff{m}, '%.3e')), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'fontsize', 12, 'fontweight', 'bold')
            set(gca,'layer','top'), set(gca, 'box', 'off'), set(gca, 'fontweight','bold','fontsize', 12)
            set(findobj(gca, 'Type', 'Line'), 'LineWidth', 2)
            xlabel('Time (sec)')
            %         savefig(strcat('#', num2str(m), ' cc',num2str(mcc(m))))
        end
%     end
    
end

%%
save(strcat('trackHLB_', date), 'F_bleached', 'F_rest', 'F_bleached_smooth', 'R2_pre', 'koff_pre', 'const_pre',...
'R2', 'koff', 'const', 'x', 'y', 'index_normal', 'index_extreme')

%%
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n',floor(tEnd/60),rem(tEnd,60));
