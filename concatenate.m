clearvars -except master_date master_series master_count isthisseriesconcatenated version
%% Series to be concatenated
load dateseries
q=1;
% Pre-injection
date_new{q}=date;   series_new{q}=series; q=q+1;

qq=q-1;

% Move to parental folders, grab what series numbers there are, and detect what series has to be concatenated to this
cd ..
parental_series_list = dir('*Series*');
for n=1:length(parental_series_list)
    if all(series==char(regexp(parental_series_list(n).name,'\d+','match')))
        series_post = char(regexp(parental_series_list(n+1).name,'\d+','match'));
    end
end
cd(strcat('Series', series))

% Post-injection
date_new{q}=date;   series_new{q}=series_post; q=q+1;

qqq=q-1;

bef = [1 qq]; % begin end
after = [qq+1 qqq];
block = 0;

%% Load data

for x=1:length(date_new)/2
    % Post-injection
    load(strcat('..\..\', date_new{x+qq}, '\Series', series_new{x+qq}, '\compile_', date_new{x+qq}, '_series', series_new{x+qq}))
    index_end_2=index_end;
    cc2=cc;
    cdk2_activity2=cdk2_activity;
    F2=F;
    F_pernuc2=F_pernuc;
    index_begin2=index_begin;
    index_end2=index_end;
    M2=M;
    V2=V;
    V_pernuc2=V_pernuc;
    mxc_nuc2=mxc_nuc;
    mxc_whole2=mxc_whole;
    mxc_nuc_single2=mxc_nuc_single;
    mxc_whole_single2=mxc_whole_single;
    cdk2_activity_single2=cdk2_activity_single;
    mxc_3D_nuc2=mxc_3D_nuc;
    mxc_3D_whole2=mxc_3D_whole;
    %     F_hist2=F_hist;
    %     F_mean2=F_mean;
    %     F_std2=F_std;
    %     F_stdmean2=F_stdmean;
    %     F_kstest_H2=F_kstest_H;
    %     F_kstest_P2=F_kstest_P;
    
    % Pre-injection
    load(strcat('..\..\', date_new{x}, '\Series', series_new{x}, '\compile_', date_new{x}, '_series', series_new{x}))
    adjustnum = cc{start_cycle}(1)-1; % If the first series didn't start at t=1, otherwise set to 0

    if adjustnum~=0
        cdk2_activity(1:adjustnum)=[];
    end
    for num=11:14
        if ~isnan(cc{num}(1))
            cc{num}(1)=cc{num}(1)-adjustnum; cc{num}(2)=cc{num}(2)-adjustnum;
            M{num}(1)=M{num}(1)-adjustnum; M{num}(2)=M{num}(2)-adjustnum;
        end
    end
    index_end=index_end-adjustnum+block+index_end2;
    index_begin=index_begin-adjustnum+block+index_begin2;
    for num=11:14
        if ~isnan(cc{num}(1)) && ~isnan(cc2{num}(1))
            conc_frame(num)=cc{num}(2)+1; % Record where it gets concatenated
            
            cc{num}(2)=cc{num}(2)+block+cc2{num}(2);
            M{num}(1)=M{num}(1)+block+M2{num}(1);
            M{num}(2)=M{num}(2)+block+M2{num}(2);
            F{num}=F{num}(:,1:end-1); % Remove the extra column
            V{num}=V{num}(:,1:end-1);
            
            if size(F2{num},1)>size(F{num},1)
                F{num}=cat(1,F{num},NaN(size(F2{num},1)-size(F{num},1),size(F{num},2)));
                V{num}=cat(1,V{num},NaN(size(V2{num},1)-size(V{num},1),size(V{num},2)));
            elseif size(F{num},1)>size(F2{num},1)
                F2{num}=cat(1,F2{num},NaN(size(F{num},1)-size(F2{num},1),size(F2{num},2)));
                V2{num}=cat(1,V2{num},NaN(size(V{num},1)-size(V2{num},1),size(V2{num},2)));
            end
            F{num}=cat(2,F{num},NaN(size(F{num},1),block),F2{num});
            V{num}=cat(2,V{num},NaN(size(V{num},1),block),V2{num});
            
            if size(F_pernuc2{num},1)>size(F_pernuc{num},1)
                F_pernuc{num}=cat(1,F_pernuc{num},NaN(size(F_pernuc2{num},1)-size(F_pernuc{num},1),size(F_pernuc{num},2)));
                V_pernuc{num}=cat(1,V_pernuc{num},NaN(size(V_pernuc2{num},1)-size(V_pernuc{num},1),size(V_pernuc{num},2)));
                mxc_nuc_single{num}=cat(1,mxc_nuc_single{num},NaN(size(mxc_nuc_single2{num},1)-size(mxc_nuc_single{num},1),size(mxc_nuc_single{num},2)));
                mxc_whole_single{num}=cat(1,mxc_whole_single{num},NaN(size(mxc_whole_single2{num},1)-size(mxc_whole_single{num},1),size(mxc_whole_single{num},2)));
                cdk2_activity_single{num}=cat(1,cdk2_activity_single{num},NaN(size(cdk2_activity_single2{num},1)-size(cdk2_activity_single{num},1),size(cdk2_activity_single{num},2)));
            elseif size(F_pernuc{num},1)>size(F_pernuc2{num},1)
                F_pernuc2{num}=cat(1,F_pernuc2{num},NaN(size(F_pernuc{num},1)-size(F_pernuc2{num},1),size(F_pernuc2{num},2)));
                V_pernuc2{num}=cat(1,V_pernuc2{num},NaN(size(V_pernuc{num},1)-size(V_pernuc2{num},1),size(V_pernuc2{num},2)));
                mxc_nuc_single2{num}=cat(1,mxc_nuc_single2{num},NaN(size(mxc_nuc_single{num},1)-size(mxc_nuc_single2{num},1),size(mxc_nuc_single2{num},2)));
                mxc_whole_single2{num}=cat(1,mxc_whole_single2{num},NaN(size(mxc_whole_single{num},1)-size(mxc_whole_single2{num},1),size(mxc_whole_single2{num},2)));
                cdk2_activity_single2{num}=cat(1,cdk2_activity_single2{num},NaN(size(cdk2_activity_single{num},1)-size(cdk2_activity_single2{num},1),size(cdk2_activity_single2{num},2)));
            end
            F_pernuc{num}=cat(2,F_pernuc{num},NaN(size(F_pernuc{num},1),block),F_pernuc2{num});
            V_pernuc{num}=cat(2,V_pernuc{num},NaN(size(V_pernuc{num},1),block),V_pernuc2{num});
            mxc_nuc_single{num}=cat(2,mxc_nuc_single{num},NaN(size(mxc_nuc_single{num},1),block),mxc_nuc_single2{num});
            mxc_whole_single{num}=cat(2,mxc_whole_single{num},NaN(size(mxc_whole_single{num},1),block),mxc_whole_single2{num});
            cdk2_activity_single{num}=cat(2,cdk2_activity_single{num},NaN(size(cdk2_activity_single{num},1),block),cdk2_activity_single2{num});
            
            mxc_nuc{num}=cat(2,mxc_nuc{num},NaN(1,block),mxc_nuc2{num});
            mxc_whole{num}=cat(2,mxc_whole{num},NaN(1,block),mxc_whole2{num});
            F_nuc{num} = F_pernuc{num}/mean(mxc_nuc{num});
        end
    end
    cdk2_adjustment = cdk2_activity2(1)-cdk2_activity(end);
    cdk2_activity=cat(2,cdk2_activity,NaN(1,block),cdk2_activity2-cdk2_adjustment);
    mxc_3D_nuc=cat(2,mxc_3D_nuc,NaN(1,block),mxc_3D_nuc2);
    mxc_3D_whole=cat(2,mxc_3D_whole,NaN(1,block),mxc_3D_whole2);
    
    mkdir(strcat('..\..\', date_new{x}, '\Series', series_new{x}, '_conc'))
    save(strcat('..\..\', date_new{x}, '\Series', series_new{x}, '_conc', '\compile_', date_new{x}, '_series', series_new{x}, '_conc'),...
        'cc', 'cdk2_activity', 'cycle_time', 'F', 'F_nuc', 'F_pernuc', 'index_begin', 'index_end',...
        'M', 'V', 'V_pernuc', 'mxc_nuc', 'mxc_whole', 'mxc_nuc_single', 'mxc_whole_single', 'cdk2_activity_single',...
        'start_cycle', 'end_cycle', 'ccHLB', 'conc_frame', 'mxc_3D_nuc', 'mxc_3D_whole')
end

%% Figures
figure; clf;
set(gcf, 'position', [10, 400, 1500, 500]);
hold on
upperbound = 1.1*max(max([F_pernuc{end_cycle}]));
yyaxis left
for num=10:14
    if ~isnan(cc{num})
        fill([M{num}(1), M{num}(1), M{num}(2), M{num}(2)]*cycle_time/60, [0, upperbound, upperbound, 0], [.9 .9 .9], 'linestyle', 'none')
    end
end
for num=10:14
    if ccHLB(num) && ~isnan(cc{num}(1))
        if num~=14
            p1=plot((cc{num}(1):M{num}(2))*cycle_time/60, nanmean(F_pernuc{num}), 'color', [0 0 0], 'linestyle', '-', 'marker', 'none');
        else
            p1=plot((cc{num}(1):cc{num}(2))*cycle_time/60, nanmean(F_pernuc{num}), 'color', [0 0 0], 'linestyle', '-', 'marker', 'none');
        end
    end
end
set(gca, 'ylim', [0 upperbound])
ylabel('Total Fluorescence/nuclei')
yyaxis right
pc=plot((index_begin:index_end)*cycle_time/60, cdk2_activity(index_begin:index_end), 'g-');
set(findobj(gca, 'Type', 'Line'), 'LineWidth', 2)
set(gca, 'FontWeight', 'bold')
ylabel('Cdk2 Activity')
set(gca, 'ylim', [0.8 3.0])
xlabel('Time (min)'), set(gca,'layer','top')
legend([pc p1], {'Cdk2 Activity', 'HLB'})
set(gca, 'box', 'off')
title('Total Fluorescence per nuclei & Cdk2 activity')

% figure; clf;
% set(gcf, 'position', [10, 400, 1500, 500]);
% hold on
% upperbound = 1.1*max(max([F_nuc{end_cycle}]));
% yyaxis left
% for num=10:14
%     if ~isnan(cc{num})
%         fill([M{num}(1), M{num}(1), M{num}(2), M{num}(2)]*cycle_time/60, [0, upperbound, upperbound, 0], [.9 .9 .9], 'linestyle', 'none')
%     end
% end
% for num=10:14
%     if ccHLB(num)&&~isnan(cc{num}(1))
%         p3=plot((cc{num}(1):M{num}(2))*cycle_time/60, nanmean(F_nuc{num}),...
%             'color', [0 0 0], 'linestyle', '-', 'marker', 'none');
%     end
% end
% set(gca, 'ylim', [0 upperbound])
% ylabel('Adjusted Total Fluorescence/nuclei')
% yyaxis right
% pc=plot((index_begin:index_end)*cycle_time/60, cdk2_activity(index_begin:index_end), 'g-');
% set(findobj(gca, 'Type', 'Line'), 'LineWidth', 2)
% set(gca, 'FontWeight', 'bold')
% ylabel('Cdk2 Activity')
% set(gca, 'ylim', [0.8 3.0])
% xlabel('Time (min)'), set(gca,'layer','top')
% legend([pc p3], {'Cdk2 Activity', 'Adjusted Total Fluorescence'})
% set(gca, 'box', 'off')
% title('Adjusted Total Fluorescence per nuclei & Cdk2 activity')
% 
% figure; clf;
% set(gcf, 'position', [10, 400, 1500, 500]);
% hold on
% upperbound2 = 1.1*max(max([V_pernuc{end_cycle}]));
% yyaxis left
% for num=10:14
%     if ~isnan(cc{num})
%         fill([M{num}(1), M{num}(1), M{num}(2), M{num}(2)]*cycle_time/60, [0, upperbound2, upperbound2, 0], [.9 .9 .9], 'linestyle', 'none')
%     end
% end
% for num=10:14
%     if ccHLB(num)&&~isnan(cc{num}(1))
%         p3=plot((cc{num}(1):M{num}(2))*cycle_time/60, nanmean(V_pernuc{num}),...
%             'color', [0 0 0], 'linestyle', '-', 'marker', 'none');
%     end
% end
% set(gca, 'ylim', [0 upperbound2])
% ylabel('V_{HLB}')
% yyaxis right
% pc=plot((index_begin:index_end)*cycle_time/60, cdk2_activity(index_begin:index_end), 'g-');
% set(findobj(gca, 'Type', 'Line'), 'LineWidth', 2)
% set(gca, 'FontWeight', 'bold')
% ylabel('Cdk2 Activity')
% set(gca, 'ylim', [0.8 3.0])
% xlabel('Time (min)'), set(gca,'layer','top')
% legend([pc p3], {'Cdk2 Activity', 'V_{HLB}'})
% set(gca, 'box', 'off')
% title('V_{HLB} & Cdk2 activity')
% 
% figure; clf
% set(gcf, 'position', [10, 400, 1500, 500]);
% hold on
% upperbound_F = 1.1*max(max([F{end_cycle}]));
% for num=10:14
%     if ~isnan(cc{num})
%         fill([M{num}(1), M{num}(1), M{num}(2), M{num}(2)]*cycle_time/60, [0, upperbound_F, upperbound_F, 0], [.9 .9 .9], 'linestyle', 'none')
%     end
% end
% for num=10:14
%     if ccHLB(num)&&~isnan(cc{num}(1))
%         for t=cc{num}(1):M{num}(2)
%             for i=1:size(F{num},1)
%                 plot(t*cycle_time/60, F{num}(i,t), 'k*')
%             end
%         end
%     end
% end
% xlabel('Time (min)')
% ylabel('Total Fluorescence (A.U.)')
% set(gca, 'ylim', [0 upperbound_F])
% set(gca,'layer','top')


