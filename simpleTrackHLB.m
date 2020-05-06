%% Measure volume and total fluorescence of HLBs for all cycles, no individual tracking
clearvars -except master_date master_series master_count isthisseriesconcatenated version
tStart = tic;
load('dateseries')
fprintf('simpleTrackHLB: %s / Series%s \n', date, series)

load(strcat('getHLB_', date, '_series', series))

for num=10:14
    if ~isnan(cc{num})
%% Assign which cell cycle this is 
iscc(10:14) = 0;
iscc(num) = 1; % assign

%% Mark objects
for t=index_begin:index_end
    I_3D_mask = squeeze(I_4D_mask(:,:,:,t));
    CCC = bwlabeln(I_3D_mask);
    CCt(:,:,:,t) = CCC;
    num_obj(t) = max(max(max(CCC)));
end

%% Track HLBs
for num=10:14
    if iscc(num)
        index_begin_track = cc{num}(1); index_end_track = M{num}(2); 
        index_end_cycle = cc{num}(2); index_begin_mitosis = M{num}(1);
    end
end

index_begin_track_new = index_begin_track;
while num_obj(index_begin_track_new)==0
    index_begin_track_new = index_begin_track_new+1;
end


%% Calculate volume
% CCt_logical = zeros(size(CCt));
V_HLB = zeros(max(num_obj), index_end_track);
F_HLB = zeros(max(num_obj), index_end_track);

if iscc(10)%||iscc(11) % Do not use trapezoidal cone with these cycles
    for t=index_begin_track:index_end_track
        for i=1:num_obj(t)
            clear xx yy zz
            clear CCt_logical CCt_logical_dilated CCt_logical_edge
            clear area V_trapcone
            z_count = 0;
            
            for z=z_begin(t):z_end(t)
                CCt_logical(:,:,z-z_begin(t)+1)=(CCt(:,:,z-z_begin(t)+1,t)==i);
            end
            
            for z=z_begin(t):z_end(t)
                area(z-z_begin(t)+1) = sum(sum(CCt_logical(:,:,z-z_begin(t)+1)));
                area_GFP(z-z_begin(t)+1) = sum(sum(double(CCt_logical(:,:,z-z_begin(t)+1)).*double(I_nondecon(:,:,z-z_begin(t)+1,t))));
            end
            for z=z_begin(t):z_end(t)
                V_trapcone(z-z_begin(t)+1) = area(z-z_begin(t)+1)*(pixelxy^2)*pixelz;
            end
            V_HLB(i,t) = sum(V_trapcone);
            F_HLB(i,t) = sum(area_GFP)/sum(area)*V_HLB(i,t); 
        end
    end
else
    for t=index_begin_track:index_end_track
        for i=1:num_obj(t)
            clear xx yy zz
            clear CCt_logical CCt_logical_dilated CCt_logical_edge
            clear area V_trapcone
            z_count = 0;
            
            for z=z_begin(t):z_end(t)
                CCt_logical(:,:,z-z_begin(t)+1)=(CCt(:,:,z-z_begin(t)+1,t)==i);
            end
            
            for z=z_begin(t):z_end(t)
                area(z-z_begin(t)+1) = sum(sum(CCt_logical(:,:,z-z_begin(t)+1)));
                area_GFP(z-z_begin(t)+1) = sum(sum(double(CCt_logical(:,:,z-z_begin(t)+1)).*double(I_nondecon(:,:,z-z_begin(t)+1,t))));
                if area(z-z_begin(t)+1)>0
                    z_count = z_count+1;
                end
            end
            if iscc(11), z_count_threshold = 2; else z_count_threshold = 2; end % 3 for later cycles
            if z_count>=z_count_threshold
                for z=z_begin(t):z_end(t)-1
                    V_trapcone(z-z_begin(t)+1) = (area(z-z_begin(t)+1)+area(z-z_begin(t)+1+1)+sqrt(area(z-z_begin(t)+1)*area(z-z_begin(t)+1+1)))/3*(pixelxy^2)*pixelz;
                end
                V_HLB(i,t) = sum(V_trapcone);
                F_HLB(i,t) = sum(area_GFP)/sum(area)*V_HLB(i,t); 
            else
                V_HLB(i,t) = 0;
                F_HLB(i,t) = 0;
            end
        end
    end
end
if ~any(num_obj(index_begin_track:index_end_track))
    V_HLB = zeros(1, index_end_track-index_begin_track+1);
    F_HLB = zeros(1, index_end_track-index_begin_track+1);
end
V_HLB(V_HLB<0) = NaN;
V_HLB(V_HLB>10) = NaN;
V_HLB(V_HLB==0) = NaN;

F_HLB(isnan(V_HLB)) = NaN;

upperbound_V = 10;
upperbound_F = 1.1*max(max(F_HLB));

% Plot all objects & nanmeans
figure; clf
hold on
for num=10:14
    if iscc(num), fill([M{num}(1), M{num}(1), M{num}(2), M{num}(2)]*cycle_time/60, [0, upperbound_V, upperbound_V, 0], [.9 .9 .9], 'linestyle', 'none'), end
end
for t=index_begin_track:index_end_track
    for i=1:num_obj(t)
        plot(t*cycle_time/60, V_HLB(i,t), 'k*'); 
    end
end
plot((index_begin_track:index_end_track)*cycle_time/60, nanmean(V_HLB(:,index_begin_track:index_end_track)))
xlabel('Time (min)')
ylabel('V_{HLB}')
if ~all(isnan(V_HLB)), set(gca, 'ylim', [0 upperbound_V]), end
set(gca, 'xlim', [index_begin_track index_end_track]*cycle_time/60)
set(gca,'layer','top')

figure; clf
hold on
for num=10:14
if iscc(num), fill([M{num}(1), M{num}(1), M{num}(2), M{num}(2)]*cycle_time/60, [0, upperbound_F, upperbound_F, 0], [.9 .9 .9], 'linestyle', 'none'), end
end
for t=index_begin_track:index_end_track
    for i=1:num_obj(t)
        plot(t*cycle_time/60, F_HLB(i,t), 'k*');
    end
end
plot((index_begin_track:index_end_track)*cycle_time/60, nanmean(F_HLB(:,index_begin_track:index_end_track)))
xlabel('Time (min)')
ylabel('Total Fluorescence (arbitrary units)')
if ~all(isnan(F_HLB)), set(gca, 'ylim', [0 upperbound_F]), end
set(gca, 'xlim', [index_begin_track index_end_track]*cycle_time/60)
set(gca,'layer','top')

V_preM = V_HLB(:,index_begin_track:index_end_cycle);
V_postM = V_HLB(:,index_begin_mitosis:index_end_track);
F_preM = F_HLB(:,index_begin_track:index_end_cycle);
F_postM = F_HLB(:,index_begin_mitosis:index_end_track);

for num=10:14
    if iscc(num), save(strcat('trackHLB_', date, '_series', series, '_cc', num2str(num)), ...
            'V_preM', 'V_postM', 'F_preM', 'F_postM', 'CCt', 'num_obj', '-v7.3'), end
end

%% Whole Tracking
% To visualize whole elements that are being tracked
% CCt_whole = zeros(size(CCt));
% CCt_whole_dilate = zeros(size(CCt));
% CCt_whole_edge = zeros(size(CCt));
% SE=strel('disk',2,0);
% for t = index_begin_track:index_end_track
%     for z=z_begin(t):z_end(t)
%         CCt_whole(:,:,z-z_begin(t)+1,t) = CCt(:,:,z-z_begin(t)+1,t)>0;
%         CCt_whole_dilate(:,:,z-z_begin(t)+1,t) = imdilate(CCt_whole(:,:,z-z_begin(t)+1,t)>0,SE);
%         CCt_whole_edge(:,:,z-z_begin(t)+1,t) = CCt_whole_dilate(:,:,z-z_begin(t)+1,t)-(CCt_whole(:,:,z-z_begin(t)+1,t));
%         
%         I_3D_wh(:,:,1,z-z_begin(t)+1,t) = 50000*uint16(CCt_whole_edge(:,:,z-z_begin(t)+1,t));
%         I_3D_wh(:,:,2,z-z_begin(t)+1,t) = 3*I_nondecon(:,:,z-z_begin(t)+1,t);
%         I_3D_wh(:,:,3,z-z_begin(t)+1,t) = 0;
%     end
% end

% for n=round(zspan/2)%1:12
% if ismember(n,1:zspan), mov = squeeze(I_3D_wh(:,:,:,n,index_begin_track:index_end_track)); end
% end
% 
% for n=round(zspan/2)%1:12
% if ismember(n,1:zspan)
%     mm = implay(mov);
%     set(findall(0,'tag','spcui_scope_framework'),'position',[150 150 820 550])
%     set(mm.Parent, 'Name', strcat('cc ', num2str(find(iscc==1))))
% end
% end

    end
end

%%
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n',floor(tEnd/60),rem(tEnd,60));

% send_msg({'919-452-2983'}, 'Script finished running', strcat('trackHLB_', date, '_series', series), 'T-Mobile')