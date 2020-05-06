clearvars -except master_date master_series master_count isthisseriesconcatenated version
load('dateseries')
ccHLB(10)=1; ccHLB(11)=1; ccHLB(12)=1; ccHLB(13)=1; ccHLB(14)=1;

load(strcat('getHLB_', date, '_series', series),...
    'index_begin', 'index_end', 'cc', 'M', 'cycle_time', 'start_cycle', 'end_cycle')

for num=10:14
    if ccHLB(num)
        if ~isnan(cc{num}(1))
            load(strcat('trackHLB_', date, '_series', series, '_cc', num2str(num)), 'V_preM', 'V_postM', 'F_preM', 'F_postM')
            V1{num}=V_preM; V2{num}=V_postM; F1{num}=F_preM; F2{num}=F_postM;
            V{num}=cat(2,V_preM,V_postM); F{num}=cat(2,F_preM,F_postM);
        end
    end
end

load(strcat('getNucleusCdk2_', date, '_series', series), 'cdk2_activity', 'mxc_nuc', 'mxc_nuc_conc', 'mxc_whole', 'mxc_whole_conc')
load(strcat('quantifyMxcCdk2HLBpernuc_', date, '_series', series))
load(strcat('getWholeMxc_', date, '_series', series))


for num=10:14
    F_nuc{num}=NaN;
    if ~isnan(cc{num})
        F_nuc{num} = F_pernuc{num}/mean(mxc_nuc{num});
    end
end

%% 3D histogram of HLB / Calculation of STD
n=1;
for num=10:14
    if ccHLB(num) && ~isnan(cc{num}(1))
        for t=cc{num}(1):M{num}(2)
            for x=1:size(F_pernuc{num},1)
                if ~isnan(F_pernuc{num}(x,t-cc{num}(1)+1))
                    F_hist(n,1)=t*cycle_time/60;
                    F_hist(n,2)=F_pernuc{num}(x,t-cc{num}(1)+1);
                    n=n+1;
                end
            end
            F_mean(t) = nanmean(F_pernuc{num}(:,t-cc{num}(1)+1));
            F_std(t) = nanstd(F_pernuc{num}(:,t-cc{num}(1)+1));
            F_stdmean(t) = F_std(t)/F_mean(t);
            if sum(~isnan(F_pernuc{num}(:,t-cc{num}(1)+1)))>1 && any(~isnan((F_pernuc{num}(:,t-cc{num}(1)+1)-F_mean(t))/F_std(t)))
            [H, P] = kstest((F_pernuc{num}(:,t-cc{num}(1)+1)-F_mean(t))/F_std(t));
            F_kstest_H(t) = H;
            F_kstest_P(t) = P;
            else
            F_kstest_H(t) = 0;
            F_kstest_P(t) = 0;
            end
        end
    end
end
nbins = [round(max(F_hist(:,1)/cycle_time*60)-min(F_hist(:,1)/cycle_time*60)+1), 200];
% figure; clf;
% set(gcf, 'position', [10, 400, 1500, 500]);
% hist3(F_hist, nbins)
% colormap jet
% % set(gca, 'xlim', [cc{start_cycle}(1) M{end_cycle}(2)]*cycle_time/60)
% % xticks([0 20 40 60 80]+((cc{start_cycle}(1)-1)*cycle_time/60)), xticklabels({'0', '20', '40', '60', '80'})
% % set(gca, 'ylim', [0 2.5e5])
% % yticks([0 0.5 1 1.5 2 2.5]*10^5)
% % set(get(gca,'child'), 'FaceColor', 'flat', 'CDataMode', 'auto', 'linestyle', 'none');
% set(get(gca,'child'), 'FaceColor', 'flat', 'CDataMode', 'auto');
% xlabel('Time (min)'), ylabel('Total Fluorescence (A.U.)')
% view(2)
% colorbar;
% % axis off

%%
% figure; clf; hold on
% set(gcf, 'position', [10, 400, 1500, 500]);
% for num=11:14
%     if ccHLB(num) && ~isnan(cc{num}(1))
%         plot((cc{num}(1):M{num}(2))*cycle_time/60, F_stdmean(cc{num}(1):M{num}(2)), 'linestyle', 'none', 'marker', '.')
%         fill([M{num}(1), M{num}(1), M{num}(2), M{num}(2)]*cycle_time/60, [0, 1, 1, 0], [0 0 0], 'linestyle', 'none')
%     end
% end
% xlabel('Time (min)'), ylabel('Coefficient of Variation')
% set(gca, 'xlim', [0 80]+((cc{start_cycle}(1)-1)*cycle_time/60)), set(gca, 'ylim', [0 1]), set(gca,'layer','top'), set(gca, 'box', 'off')
% set(gca, 'fontweight','bold','fontsize', 12)
% yticks([0.2 0.4 0.6 0.8 1]), yticklabels({'0.2', '0.4', '0.6', '0.8', '1.0'})
% xticks([0 20 40 60 80]+((cc{start_cycle}(1)-1)*cycle_time/60)), xticklabels({'0', '20', '40', '60', '80'})
% set(gca,'layer','top'), set(gca, 'box', 'off')
%%
% figure; clf; hold on
% for num=11:14
%     if ccHLB(num) && ~isnan(cc{num}(1))
%         plot([0 M{num}(2)+1]*cycle_time/60, [0.05 0.05], 'linestyle', '-', 'color', 'k')
%         for t=cc{num}(1):M{num}(2)
%             plot(t*cycle_time/60, F_kstest_P(t), 'k.')
%         end
%         fill([M{num}(1), M{num}(1), M{num}(2), M{num}(2)]*cycle_time/60, [0, 1.1*max(F_kstest_H), 1.1*max(F_kstest_H), 0], [.9 .9 .9], 'linestyle', 'none')
%     end
% end
% xlabel('Time (min)'), ylabel('P value'), title('Kolmogorov-Smirnov Test')
% set(gca, 'ylim', [0 1]), set(gca, 'xlim', [0 (M{end_cycle}(2)+1)*cycle_time/60]), set(gca,'layer','top'), set(gca, 'box', 'off')

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
        p1=plot((cc{num}(1):M{num}(2))*cycle_time/60, nanmean(F_pernuc{num}), 'color', [0 0 0], 'linestyle', '-', 'marker', 'none');
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

figure; clf;
set(gcf, 'position', [10, 400, 1500, 500]);
hold on
upperbound2 = 1.1*max(max([V_pernuc{end_cycle}]));
yyaxis left
for num=10:14
    if ~isnan(cc{num})
        fill([M{num}(1), M{num}(1), M{num}(2), M{num}(2)]*cycle_time/60, [0, upperbound2, upperbound2, 0], [.9 .9 .9], 'linestyle', 'none')
    end
end
for num=10:14
    if ccHLB(num)&&~isnan(cc{num}(1))
        p3=plot((cc{num}(1):M{num}(2))*cycle_time/60, nanmean(V_pernuc{num}),...
            'color', [0 0 0], 'linestyle', '-', 'marker', 'none');
    end
end
set(gca, 'ylim', [0 upperbound2])
ylabel('V_{HLB}')
yyaxis right
pc=plot((index_begin:index_end)*cycle_time/60, cdk2_activity(index_begin:index_end), 'g-');
set(findobj(gca, 'Type', 'Line'), 'LineWidth', 2)
set(gca, 'FontWeight', 'bold')
ylabel('Cdk2 Activity')
set(gca, 'ylim', [0.8 3.0])
xlabel('Time (min)'), set(gca,'layer','top')
legend([pc p3], {'Cdk2 Activity', 'V_{HLB}'})
set(gca, 'box', 'off')
title('V_{HLB} & Cdk2 activity')


figure; clf; hold on
set(gcf, 'position', [10, 400, 1500, 500]);
hold on
upperbound = 1.1*max(max([F_pernuc{end_cycle}])); upperbound2 = 3*max([mxc_nuc{:}]);
yyaxis left
for num=10:14
    if ~isnan(cc{num})
        fill([M{num}(1), M{num}(1), M{num}(2), M{num}(2)]*cycle_time/60, [0, upperbound, upperbound, 0], [.9 .9 .9], 'linestyle', 'none')
    end
end
for num=10:14
    if ccHLB(num)&&~isnan(cc{num}(1))
        p3=plot((cc{num}(1):M{num}(2))*cycle_time/60, nanmean(F_pernuc{num}),'color', [0 0 0], 'linestyle', '-', 'marker', 'none');
    end
end
set(gca, 'ylim', [0 upperbound])
ylabel('Total Fluorescence/nuclei (arbitrary unit)')
yyaxis right
p4=plot((index_begin:index_end)*cycle_time/60, mxc_nuc_conc(index_begin:index_end),'color', [0 1 0], 'linestyle', '-', 'marker', 'none');
set(findobj(gca, 'Type', 'Line'), 'LineWidth', 2)
set(gca, 'FontWeight', 'bold')
ylabel('Mxc Concentration (arbitrary unit)')
set(gca, 'ylim', [0 upperbound2])
xlabel('Time (min)'), set(gca,'layer','top')
legend([p4 p3], {'Mxc Concentration', 'Total Fluorescence'})
set(gca, 'box', 'off')
title('Total Fluorescence per nuclei & Mxc Concentration')

%%
save(strcat('compile_', date, '_series', series), 'cc', 'cdk2_activity', 'cycle_time', 'F', 'F_pernuc',...
    'index_begin', 'index_end',...
    'M', 'V', 'V_pernuc', 'mxc_nuc', 'mxc_whole', 'mxc_nuc_single', 'mxc_whole_single', 'cdk2_activity_single',...
    'F_hist', 'F_mean', 'F_std', 'F_stdmean', 'F_kstest_H', 'F_kstest_P', 'start_cycle', 'end_cycle', 'ccHLB',...
    'mxc_3D_nuc', 'mxc_3D_whole')
