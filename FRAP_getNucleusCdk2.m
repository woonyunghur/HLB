%% Nuclear mask and Cdk2 activity quantification for FRAP experiments
clear; tStart = tic;
date = '170317_FRAP_2';

load(strcat('getHLB_', date)), load(strcat('trackHLB_', date))

%% Load images (red channel)
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
        files_pre_red(:,index) = dir([dir_FRAP filesep strcat('*FRAP_', sprintf('%03d', n), '*pre*t0_z0_ch01*')]);
        files_pb_red(:,index) = dir([dir_FRAP filesep strcat('*FRAP_', sprintf('%03d', n), '*Pb*z0_ch01*')]);
    end
end

for m=1:index_length
    for t=1:pb_length+1
        if t==1
            J_pre=imread([dir_FRAP filesep files_pre_red(:,m).name]);
            J_comp{m}(:,:,t)=J_pre;
        else
            J_pb=imread([dir_FRAP filesep files_pb_red(t-1,m).name]);
            J_comp{m}(:,:,t)=J_pb;
        end
    end
    
    % Nuclear mask
    for t=1:pb_length+1
        I=I_comp{m}(:,:,t);
        I = imgaussfilt(I, 0.5); % Blur to get smoother edge

        J=J_comp{m}(:,:,t);
        
        for i=1:size(I,1)
            for j=1:size(I,2)
                if m<=cc_index(11)
                    threshold=maxI{m}*0.03;
                elseif m<=cc_index(12)
                    threshold=maxI{m}*0.03; 
                elseif m<=cc_index(13)
                    threshold=maxI{m}*0.03;
                else %cc14
                    threshold=maxI{m}*0.01;
                end
                if(I(i,j)>threshold)
                    I_nuc(i,j)=1;
                else
                    I_nuc(i,j)=0;
                end
            end
        end
        
        I_nuc = bwmorph(I_nuc, 'clean');
        for num_maj = 1:5
            I_nuc = bwmorph(I_nuc, 'majority');
        end
        I_nuc = imfill(I_nuc, 'holes');
        
        I_nuc = imerode(I_nuc,strel('disk',12,0));
        I_nuc = imdilate(I_nuc,strel('disk',12,0));
        I_center = imerode(I_nuc,strel('disk',12,0)); % Use smaller nuclear circle
        
        I_nucdilate=imdilate(I_nuc,strel('disk',4,0));
        I_nucdilate2=imdilate(I_nuc,strel('disk',9,0));
        I_cytring=I_nucdilate2-I_nucdilate;

        I_nucedge=imdilate(I_nuc,strel('disk',1,0))-I_nuc; % To display, not used for calculations
        
        I_centerdilate=imdilate(I_center,strel('disk',2,0));
        I_nucring=I_centerdilate-I_center;

        area_nucleus = sum(sum(I_center));
        area_cytoplasm = sum(sum(I_cytring));
        cdk2_total_nucleus = sum(sum(I_nuc.*double(J)));
        cdk2_total_cytoplasm = sum(sum(I_cytring.*double(J)));
        
        cdk2_nucleus = cdk2_total_nucleus./area_nucleus;
        cdk2_cytoplasm = cdk2_total_cytoplasm./area_cytoplasm;
        
        cdk2_activity{m}(t) = cdk2_cytoplasm./cdk2_nucleus;

        if t==1 % Left = green channel, Right = red channel
            I_3D{m}(:,:,1,t)= horzcat(50000*uint16(I_nucedge),50000*uint16(I_nucedge)+5*J_comp{m}(:,:,t));
            I_3D{m}(:,:,2,t)= horzcat(50000*uint16(I_nucedge)+3*I_comp{m}(:,:,t),50000*uint16(I_nucedge));
            I_3D{m}(:,:,3,t)= horzcat(50000*uint16(I_nucring),50000*uint16(I_nucring));
        else
            I_3D{m}(:,:,1,t+1)= horzcat(50000*uint16(I_nucedge),50000*uint16(I_nucedge)+5*J_comp{m}(:,:,t));
            I_3D{m}(:,:,2,t+1)= horzcat(50000*uint16(I_nucedge)+3*I_comp{m}(:,:,t),50000*uint16(I_nucedge));
            I_3D{m}(:,:,3,t+1)= horzcat(50000*uint16(I_nucring),50000*uint16(I_nucring));
        end
    end
    I_3D{m}(:,:,1,2)= 0;
    I_3D{m}(:,:,2,2)= horzcat(I_bleach{m},zeros(500,800));
    I_3D{m}(:,:,3,2)= 0;

    mov{m} = squeeze(I_3D{m});

    mm=implay(mov{m}); set(findall(0,'tag','spcui_scope_framework'),'position',[150 150 1620 550]);
    set(mm.Parent, 'Name', strcat('cc', num2str(mcc(m)), ' Set #', num2str(m)))
    
    fig = figure; fig.Name = strcat('#', num2str(m), ' cc', num2str(mcc(m))); set(gcf, 'position', [10, 400, 1500, 400]);
    % Nonbleached & bleached
    subplot(1,3,1), hold on
    plot((2:1:pb_length+1)*cycle_time, F_bleached{m}(2:end), 'color', [1 0.8 0.8]);
    pk=plot((1:1:pb_length+1)*cycle_time, F_rest{m}, 'k-');
    pr=plot((2:1:pb_length+1)*cycle_time, F_bleached_smooth{m}, 'r-');
    set(findobj(gca, 'Type', 'Line'), 'LineWidth', 2)
    set(gca, 'ylim', [0 inf]), set(gca, 'xlim', [0 pb_length*cycle_time])
    set(gca,'layer','top'), set(gca, 'box', 'off'), set(gca, 'fontweight','bold','fontsize', 12)
    legend([pr pk], {'Bleached', 'Non-bleached'}, 'location', 'southeast')
    xlabel('Time (sec)')

    subplot(1,3,2), hold on
    plot(x{m}(index_normal{m}), y{m}(index_normal{m}), 'linestyle', 'none', 'marker', '*', 'color', [0 0 1])
    plot(x{m}(index_extreme{m}), y{m}(index_extreme{m}), 'linestyle', 'none', 'marker', '*', 'color', [0.7 0.7 1])
    plot(x{m}, -koff_pre{m}*x{m}+const_pre{m}, 'color', [0.7 1 0.7])
    plot(x{m}(index_normal{m}), -koff{m}*x{m}(index_normal{m})+const{m}, 'color', [0.2 1 0.2])
    set(gca, 'xlim', [0 pb_length*cycle_time])
    xlim = get(gca, 'xlim'); ylim = get(gca, 'ylim');
    text(xlim(2),ylim(2), strcat('R^2=', num2str(R2{m}, '%.4f')), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'fontsize', 12, 'fontweight', 'bold')
    text(xlim(2),ylim(2)-(ylim(2)-ylim(1))*0.1, strcat('k_{off}=', num2str(koff{m}, '%.3e')), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'fontsize', 12, 'fontweight', 'bold')
    set(gca,'layer','top'), set(gca, 'box', 'off'), set(gca, 'fontweight','bold','fontsize', 12)
    set(findobj(gca, 'Type', 'Line'), 'LineWidth', 2)
    xlabel('Time (sec)')

    if R2{m}>0.6 && (max(cdk2_activity{m})-min(cdk2_activity{m}))<0.1
        isusefulset(m)=1;
    else
        isusefulset(m)=0;
    end
    subplot(1,3,3), hold on
    plot((1:1:pb_length+1)*cycle_time, cdk2_activity{m}, 'r-');
    set(findobj(gca, 'Type', 'Line'), 'LineWidth', 2)
    set(gca, 'ylim', [0.3 1.2]), set(gca, 'xlim', [0 (pb_length+1)*cycle_time])
    xlim = get(gca, 'xlim'); ylim = get(gca, 'ylim'); 
    if (max(cdk2_activity{m})-min(cdk2_activity{m}))<0.1
        text(xlim(2),ylim(2), strcat('Stable Cdk2'), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'fontsize', 12, 'fontweight', 'bold')
    end
    set(gca,'layer','top'), set(gca, 'box', 'off'), set(gca, 'fontweight','bold','fontsize', 12)
    xlabel('Time (sec)')
%     savefig(strcat('#', num2str(m), ' cc',num2str(mcc(m))))

end

%%
linecolor{11}=[1 0 0]; linecolor{12}=[1 0.5 0]; linecolor{13}=[0 0.8 0]; linecolor{14}=[0 0 1];
figure; clf; hold on
for m=1:index_length
    if isusefulset(m) % Plot only if R2>0.6 and Cdk2 is constant
        usefulcdk2(m)=mean(cdk2_activity{m});
        usefulkoff(m)=koff{m};
        scatter(usefulcdk2(m), usefulkoff(m), [], linecolor{mcc(m)}, 'filled')
    elseif (max(cdk2_activity{m})-min(cdk2_activity{m}))<0.1 % Plot sets with smaller R2 values
        
    end
end
set(gca, 'xlim', [0.3 1.2]), set(gca, 'ylim', [0 1.1*max(usefulkoff)])%, set(gca, 'yscale', 'log')
yticks([0 1e-3 5e-3 1e-2 2e-2 4e-2])
yticklabels({'0', '1.0e-3','5.0e-3', '1.0e-2','2.0e-2', '4.0e-2'})
xlim = get(gca, 'xlim'); ylim = get(gca, 'ylim'); 
for num=11:14
    if ismember(num, mcc(isusefulset>0))
        text(xlim(2)-(xlim(2)-xlim(1))*0.01,ylim(2)-(ylim(2)-ylim(1))*0.02-(num-min(mcc(isusefulset>0)))*(ylim(2)-ylim(1))*0.07, strcat('Cycle ', num2str(num)), 'color', linecolor{num}, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'fontsize', 12, 'fontweight', 'bold')
    end
end
xlabel('Cdk2 Activity (A.U.)')
ylabel('k_{off}')
savefig(strcat('k_{off} vs Cdk2, ', date))

%%
save(strcat('getNucleusCdk2_', date), 'cdk2_activity', 'isusefulset')

%%
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n',floor(tEnd/60),rem(tEnd,60));