%% Compile selected dataset
clear;
files=dir('*mat*'); 

for k=1:length(files)
    load(files(k).name)
    F_green_HLB_single_comp{k} = F_green_HLB_single;
    V_red_HLB_single_comp{k} = V_red_HLB_single;
    F_red_HLB_single_comp{k} = F_red_HLB_single;
    transcription_comp(k) = transcription;
    V_HLB_comp(k) = V_HLB;
    F_green_HLB_single_paired_comp{k}=F_green_HLB_single_paired;
    F_green_HLB_single_unpaired_comp{k}=F_green_HLB_single_unpaired;
    F_red_HLB_single_paired_comp{k}=F_red_HLB_single_paired;
    F_red_HLB_single_unpaired_comp{k}=F_red_HLB_single_unpaired;
    F_green_pernuc_comp{k} = F_green_pernuc;
    F_red_pernuc_comp{k} = F_red_pernuc;
    cycle_comp{k} = cycle;
    ncratio_comp{k} = ncratio;
end

% figure; clf
% plot(V_HLB_comp, transcription_comp, 'linestyle', 'none', 'color', 'k', 'marker', '.', 'markersize', 15)
% xlabel('V_HLB (\mum^3)'), ylabel('Transcription (A.U.)')

%% Plot by the cell cycle
color{11} = [1 0 0]; color{12} = [1 0 0]; color{13} = 'g'; color{14} = 'b'; %color{14} = 'b';
% color{11} = [1 0 0]; color{12} = 'none'; color{13} = 'none'; color{14} = 'b'; %color{14} = 'b';

% If you want to quantify single HLB instead of HLB/nuclei, activate these two lines
F_green_pernuc_comp = F_green_HLB_single_comp;
V_red_pernuc_comp = V_red_HLB_single_comp;

figure; clf
hold on
for k=1:length(files)
    color_rand{k} = rand(1,3);
    for p=1:length(F_green_pernuc_comp{k})
        if ~isnan(F_green_pernuc_comp{k}(p)) && ~isnan(V_red_pernuc_comp{k}(p))
            % Single HLB/nuc
            plot(V_red_pernuc_comp{k}(p), F_green_pernuc_comp{k}(p), 'linestyle', 'none', 'color', color{cycle_comp{k}},...
                'marker', '.', 'markersize', 15)
%             if cycle_comp{k}~=14, plot(V_red_pernuc_comp{k}(p), F_green_pernuc_comp{k}(p), 'linestyle', 'none', 'color', color_rand{k},...
%                     'marker', '.', 'markersize', 15), end
        end
    end
end

q=0; q2=0;
for k=1:length(files)
    qq=0;
    for p=1:length(F_green_pernuc_comp{k})
        if ~isnan(F_green_pernuc_comp{k}(p)) && ~isnan(V_red_pernuc_comp{k}(p))
            if cycle_comp{k} ~= 14
                q2=q2+1;
                V_HLB_single_pre14(q2)=V_red_pernuc_comp{k}(p);
                transcription_single_pre14(q2)=F_green_pernuc_comp{k}(p);
            end
            q=q+1;
            V_HLB_single(q)=V_red_pernuc_comp{k}(p);
            transcription_single(q)=F_green_pernuc_comp{k}(p);
            qq=qq+1;
            V_HLB_set{k}(qq)=V_red_pernuc_comp{k}(p);
            transcription_set{k}(qq)=F_green_pernuc_comp{k}(p);
        end
    end
    % Trimmed mean of each embryo
%     plot(trimmean(V_HLB_set{k},25), trimmean(transcription_set{k}, 25), 'linestyle', 'none', 'color', color_rand{k}, 'marker', '*', 'markersize', 15)
    V_set_mean(k) = trimmean(V_HLB_set{k},25);
    transcription_set_mean(k) = trimmean(transcription_set{k}, 25);
    
    if length(V_HLB_set{k})>=5
        [f_set,gof_set]=fit((V_HLB_set{k})',(transcription_set{k})','poly1');
        R2_set(k) = gof_set.rsquare;
        slope_set{k}=f_set.p1;
        const_set{k}=f_set.p2;
        x_set{k}=linspace(prctile(V_HLB_set{k}, 20), prctile(V_HLB_set{k}, 80));
        if R2_set(k)>=0.3
            setstyle='-'; 
        else
            setstyle='--'; 
        end
        plot(x_set{k}, f_set.p1*x_set{k}+f_set.p2, 'linestyle', setstyle, 'linewidth', 1, 'color', color{cycle_comp{k}})
%         if cycle_comp{k}~=14, plot(x_set{k}, f_set.p1*x_set{k}+f_set.p2, 'linestyle', setstyle, 'linewidth', 2, 'color', color_rand{k}), end
    else
        slope_set{k}=NaN;
        const_set{k}=NaN;
        x_set{k}=NaN;
    end
end

[f,gof]=fit((V_HLB_single)',(transcription_single)','poly1');
R2 = gof.rsquare
slope = f.p1;
const = f.p2;
x=linspace(0,5);

[f2,gof2]=fit((V_HLB_single_pre14)',(transcription_single_pre14)','poly1');
R2_2 = gof2.rsquare
slope_2 = f2.p1;
const_2 = f2.p2;
x=linspace(0,5);


plot(x, f.p1*x+f.p2, 'linestyle', '-', 'linewidth', 2, 'color', [0 0 0])
% plot(x, f2.p1*x+f2.p2, 'linestyle', '-', 'linewidth', 2, 'color', [0 0 0])

xlabel('V_{HLB} (\mum^3)'), ylabel('Transcription (A.U.)')

xlim = get(gca, 'xlim'); ylim = get(gca, 'ylim'); 
for num=11:14
text(xlim(1)+(xlim(2)-xlim(1))*0.05,ylim(2)*0.95-(ylim(2)-ylim(1))*0.05*(num-11), strcat('Cycle ', num2str(num)), 'color', color{num}, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'fontsize', 12, 'fontweight', 'bold')
end

%% Plot by the cell cycle - Total Fluorescence
color{11} = [1 0 0]; color{12} = [1 0 0]; color{13} = 'g'; color{14} = 'b'; %color{14} = 'b';
% color{11} = [1 0 0]; color{12} = 'none'; color{13} = 'none'; color{14} = 'b'; %color{14} = 'b';

% If you want to quantify single HLB instead of HLB/nuclei, activate these two lines
F_green_pernuc_comp = F_green_HLB_single_comp;
F_red_pernuc_comp = F_red_HLB_single_comp; % Use total MPM-2 fluorescence instead of volume

figure; clf
hold on
for k=1:length(files)
    color_rand{k} = rand(1,3);
    for p=1:length(F_green_pernuc_comp{k})
        if ~isnan(F_green_pernuc_comp{k}(p)) && ~isnan(F_red_pernuc_comp{k}(p))
            % Single HLB/nuc
            plot(F_red_pernuc_comp{k}(p), F_green_pernuc_comp{k}(p), 'linestyle', 'none', 'color', color{cycle_comp{k}},...
                'marker', '.', 'markersize', 15)
%             if cycle_comp{k}~=14, plot(V_red_pernuc_comp{k}(p), F_green_pernuc_comp{k}(p), 'linestyle', 'none', 'color', color_rand{k},...
%                     'marker', '.', 'markersize', 15), end
        end
    end
end

q=0; q2=0;
for k=1:length(files)
    qq=0;
    for p=1:length(F_green_pernuc_comp{k})
        if ~isnan(F_green_pernuc_comp{k}(p)) && ~isnan(F_red_pernuc_comp{k}(p))
            if cycle_comp{k} ~= 14
                q2=q2+1;
                V_HLB_single_pre14(q2)=F_red_pernuc_comp{k}(p);
                transcription_single_pre14(q2)=F_green_pernuc_comp{k}(p);
            end
            q=q+1;
            V_HLB_single(q)=F_red_pernuc_comp{k}(p);
            transcription_single(q)=F_green_pernuc_comp{k}(p);
            qq=qq+1;
            V_HLB_set{k}(qq)=F_red_pernuc_comp{k}(p);
            transcription_set{k}(qq)=F_green_pernuc_comp{k}(p);
        end
    end
    % Trimmed mean of each embryo
%     plot(trimmean(V_HLB_set{k},25), trimmean(transcription_set{k}, 25), 'linestyle', 'none', 'color', color_rand{k}, 'marker', '*', 'markersize', 15)
    V_set_mean(k) = trimmean(V_HLB_set{k},25);
    transcription_set_mean(k) = trimmean(transcription_set{k}, 25);
    
    if length(V_HLB_set{k})>=5
        [f_set,gof_set]=fit((V_HLB_set{k})',(transcription_set{k})','poly1');
        R2_set(k) = gof_set.rsquare;
        slope_set{k}=f_set.p1;
        const_set{k}=f_set.p2;
        x_set{k}=linspace(prctile(V_HLB_set{k}, 20), prctile(V_HLB_set{k}, 80));
        if R2_set(k)>=0.3
            setstyle='-'; 
        else
            setstyle='--'; 
        end
        plot(x_set{k}, f_set.p1*x_set{k}+f_set.p2, 'linestyle', setstyle, 'linewidth', 1, 'color', color{cycle_comp{k}})
%         if cycle_comp{k}~=14, plot(x_set{k}, f_set.p1*x_set{k}+f_set.p2, 'linestyle', setstyle, 'linewidth', 2, 'color', color_rand{k}), end
    else
        slope_set{k}=NaN;
        const_set{k}=NaN;
        x_set{k}=NaN;
    end
end

[f,gof]=fit((V_HLB_single)',(transcription_single)','poly1');
R2 = gof.rsquare
slope = f.p1;
const = f.p2;
x=linspace(0,5);

[f2,gof2]=fit((V_HLB_single_pre14)',(transcription_single_pre14)','poly1');
R2_2 = gof2.rsquare
slope_2 = f2.p1;
const_2 = f2.p2;
x=linspace(0,5);


plot(x, f.p1*x+f.p2, 'linestyle', '-', 'linewidth', 2, 'color', [0 0 0])
% plot(x, f2.p1*x+f2.p2, 'linestyle', '-', 'linewidth', 2, 'color', [0 0 0])

xlabel('F_{HLB} (A.U.)'), ylabel('Transcription (A.U.)')

xlim = get(gca, 'xlim'); ylim = get(gca, 'ylim'); 
for num=11:14
text(xlim(1)+(xlim(2)-xlim(1))*0.05,ylim(2)*0.95-(ylim(2)-ylim(1))*0.05*(num-11), strcat('Cycle ', num2str(num)), 'color', color{num}, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'fontsize', 12, 'fontweight', 'bold')
end

% figure; hold on
% for k=1:length(files)
%     plot(cycle_comp{k}, ncratio_comp{k}, 'linestyle', 'none', 'color', color{cycle_comp{k}}, 'marker', '.', 'markersize', 15)
%     if cycle_comp{k}==14 && ncratio_comp{k}<0.4
%         files(k).name
%     end
% end

%% Paired vs unpaired

figure; clf
hold on

for k=1:length(files)
    color_rand{k} = rand(1,3);
    for p=1:length(F_green_HLB_single_paired_comp{k})
        if ~isnan(F_green_HLB_single_paired_comp{k}(p)) && ~isnan(F_red_HLB_single_paired_comp{k}(p))
            % Single HLB
            if cycle_comp{k}~=14
                plot(F_red_HLB_single_paired_comp{k}(p), F_green_HLB_single_paired_comp{k}(p), 'linestyle', 'none', 'color', 'r',...
                    'marker', '.', 'markersize', 15)
            end
        elseif ~isnan(F_green_HLB_single_unpaired_comp{k}(p)) && ~isnan(F_red_HLB_single_unpaired_comp{k}(p))
            % Single HLB
            if cycle_comp{k}~=14
                plot(F_red_HLB_single_unpaired_comp{k}(p), F_green_HLB_single_unpaired_comp{k}(p), 'linestyle', 'none', 'color', 'b',...
                    'marker', '.', 'markersize', 15), 
            end
        end
    end
end

% Paired
q=0; q2=0;
for k=1:length(files)
    qq=0; 
    for p=1:length(F_green_HLB_single_paired_comp{k})
        if ~isnan(F_green_HLB_single_paired_comp{k}(p)) && ~isnan(F_red_HLB_single_paired_comp{k}(p))
            if cycle_comp{k} ~= 14
                q2=q2+1;
                V_HLB_single_paired_pre14(q2)=F_red_HLB_single_paired_comp{k}(p);
                transcription_single_paired_pre14(q2)=F_green_HLB_single_paired_comp{k}(p);
            end
            q=q+1;
            V_HLB_single_paired(q)=F_red_HLB_single_paired_comp{k}(p);
            transcription_single_paired(q)=F_green_HLB_single_paired_comp{k}(p);
            qq=qq+1;
            V_HLB_paired_set{k}(qq)=F_red_HLB_single_paired_comp{k}(p);
            transcription_paired_set{k}(qq)=F_green_HLB_single_paired_comp{k}(p);
        end
    end
    % Trimmed mean of each embryo
%     plot(trimmean(V_HLB_paired_set{k},25), trimmean(transcription_paired_set{k}, 25), 'linestyle', 'none', 'color', color_rand{k}, 'marker', '*', 'markersize', 15)
%     V_set_mean(k) = trimmean(V_HLB_paired_set{k},25);
%     transcription_paired_set_mean(k) = trimmean(transcription_paired_set{k}, 25);
    
    if qq>=5 && length(V_HLB_paired_set{k})>=5
        [f_paired_set,gof_set]=fit((V_HLB_paired_set{k})',(transcription_paired_set{k})','poly1');
        R2_paired_set(k) = gof_set.rsquare;
        slope_paired_set{k}=f_paired_set.p1;
        const_paired_set{k}=f_paired_set.p2;
        x_paired_set{k}=linspace(prctile(V_HLB_paired_set{k}, 20), prctile(V_HLB_paired_set{k}, 80));
        if R2_paired_set(k)>=0.3
            setstyle='-'; 
        else
            setstyle='--'; 
        end
%         plot(x_set{k}, f_set.p1*x_set{k}+f_set.p2, 'linestyle', setstyle, 'linewidth', 1, 'color', color{cycle_comp{k}})
        if cycle_comp{k}~=14, plot(x_paired_set{k}, f_paired_set.p1*x_paired_set{k}+f_paired_set.p2, 'linestyle', setstyle, 'linewidth', 2, 'color', [1 0 0]), end
    else
        slope_paired_set{k}=NaN;
        const_paired_set{k}=NaN;
        x_paired_set{k}=NaN;
    end
end

[f3,gof3]=fit((V_HLB_single_paired)',(transcription_single_paired)','poly1');
R2_paired = gof3.rsquare
slope_paired = f3.p1;
const_paired = f3.p2;
x=linspace(0,5);

[f4,gof4]=fit((V_HLB_single_paired_pre14)',(transcription_single_paired_pre14)','poly1');
R2_2_paired = gof4.rsquare
slope_2 = f4.p1;
const_2 = f4.p2;
x=linspace(0,5);


% plot(x, f3.p1*x+f3.p2, 'linestyle', '-', 'linewidth', 2, 'color', [1 0 0]) % all cycles
plot(x, f4.p1*x+f4.p2, 'linestyle', '-', 'linewidth', 2, 'color', [1 0 0]) % exclude cycle 13

% Unpaired
q=0; q2=0;
for k=1:length(files)
    qq=0; 
    for p=1:length(F_green_HLB_single_unpaired_comp{k})
        if ~isnan(F_green_HLB_single_unpaired_comp{k}(p)) && ~isnan(F_red_HLB_single_unpaired_comp{k}(p))
            if cycle_comp{k} ~= 14
                q2=q2+1;
                V_HLB_single_unpaired_pre14(q2)=F_red_HLB_single_unpaired_comp{k}(p);
                transcription_single_unpaired_pre14(q2)=F_green_HLB_single_unpaired_comp{k}(p);
            end
            q=q+1;
            V_HLB_single_unpaired(q)=F_red_HLB_single_unpaired_comp{k}(p);
            transcription_single_unpaired(q)=F_green_HLB_single_unpaired_comp{k}(p);
            qq=qq+1;
            V_HLB_unpaired_set{k}(qq)=F_red_HLB_single_unpaired_comp{k}(p);
            transcription_unpaired_set{k}(qq)=F_green_HLB_single_unpaired_comp{k}(p);
        end
    end
    % Trimmed mean of each embryo
%     plot(trimmean(V_HLB_unpaired_set{k},25), trimmean(transcription_unpaired_set{k}, 25), 'linestyle', 'none', 'color', color_rand{k}, 'marker', '*', 'markersize', 15)
%     V_set_mean(k) = trimmean(V_HLB_unpaired_set{k},25);
%     transcription_unpaired_set_mean(k) = trimmean(transcription_unpaired_set{k}, 25);
    
    if qq>=5 && length(V_HLB_unpaired_set{k})>=5
        [f_unpaired_set,gof_set]=fit((V_HLB_unpaired_set{k})',(transcription_unpaired_set{k})','poly1');
        R2_unpaired_set(k) = gof_set.rsquare;
        slope_unpaired_set{k}=f_unpaired_set.p1;
        const_unpaired_set{k}=f_unpaired_set.p2;
        x_unpaired_set{k}=linspace(prctile(V_HLB_unpaired_set{k}, 20), prctile(V_HLB_unpaired_set{k}, 80));
        if R2_unpaired_set(k)>=0.3
            setstyle='-'; 
        else
            setstyle='--'; 
        end
%         plot(x_set{k}, f_set.p1*x_set{k}+f_set.p2, 'linestyle', setstyle, 'linewidth', 1, 'color', color{cycle_comp{k}})
        if cycle_comp{k}~=14, plot(x_unpaired_set{k}, f_unpaired_set.p1*x_unpaired_set{k}+f_unpaired_set.p2, 'linestyle', setstyle, 'linewidth', 2, 'color', [0 0 1]), end
    else
        slope_unpaired_set{k}=NaN;
        const_unpaired_set{k}=NaN;
        x_unpaired_set{k}=NaN;
    end
end

[f5,gof5]=fit((V_HLB_single_unpaired)',(transcription_single_unpaired)','poly1');
R2_unpaired = gof5.rsquare
slope_unpaired = f5.p1;
const_unpaired = f5.p2;
x=linspace(0,6);

[f6,gof6]=fit((V_HLB_single_unpaired_pre14)',(transcription_single_unpaired_pre14)','poly1');
R2_2_unpaired = gof6.rsquare
slope_2 = f6.p1;
const_2 = f6.p2;
x=linspace(0,6);


% plot(x, f3.p1*x+f5.p2, 'linestyle', '-', 'linewidth', 2, 'color', [0 0 1]) % all cycles
plot(x, f4.p1*x+f6.p2, 'linestyle', '-', 'linewidth', 2, 'color', [0 0 1]) % exclude cycle 14




xlabel('F_{HLB} (A.U.)'), ylabel('Transcription (A.U.)')



% 
% 
% %%
% for kkk=1:12
%     bin(kkk)=0;
%     transbin{kkk}=0;
% end
% for kk=1:length(V_HLB_single_pre14)
% %     if V_HLB_single_pre14(kk)<=0.5
% %         bin(1)=bin(1)+1;
% %         transbin{1}(kk)=transcription_single_pre14(kk);
% %     elseif V_HLB_single_pre14(kk)>0.5 && V_HLB_single_pre14(kk)<=1
% %         bin(2)=bin(2)+1;
% %         transbin{2}(kk)=transcription_single_pre14(kk);
% %     elseif V_HLB_single_pre14(kk)>1 && V_HLB_single_pre14(kk)<=1.5
% %         bin(3)=bin(3)+1;
% %         transbin{3}(kk)=transcription_single_pre14(kk);
% %     elseif V_HLB_single_pre14(kk)>1.5 && V_HLB_single_pre14(kk)<=2
% %         bin(4)=bin(4)+1;
% %         transbin{4}(kk)=transcription_single_pre14(kk);
% %     elseif V_HLB_single_pre14(kk)>2 && V_HLB_single_pre14(kk)<=2.5
% %         bin(5)=bin(5)+1;
% %         transbin{5}(kk)=transcription_single_pre14(kk);
% %     elseif V_HLB_single_pre14(kk)>2.5 && V_HLB_single_pre14(kk)<=3
% %         bin(6)=bin(6)+1;
% %         transbin{6}(kk)=transcription_single_pre14(kk);
% %     elseif V_HLB_single_pre14(kk)>3 && V_HLB_single_pre14(kk)<=3.5
% %         bin(7)=bin(7)+1;
% %         transbin{7}(kk)=transcription_single_pre14(kk);
% %     elseif V_HLB_single_pre14(kk)>3.5 && V_HLB_single_pre14(kk)<=4
% %         bin(8)=bin(8)+1;
% %         transbin{8}(kk)=transcription_single_pre14(kk);
% %     elseif V_HLB_single_pre14(kk)>4 && V_HLB_single_pre14(kk)<=4.5
% %         bin(9)=bin(9)+1;
% %         transbin{9}(kk)=transcription_single_pre14(kk);
% %     elseif V_HLB_single_pre14(kk)>4.5 && V_HLB_single_pre14(kk)<=5
% %         bin(10)=bin(10)+1;
% %         transbin{10}(kk)=transcription_single_pre14(kk);
% %     elseif V_HLB_single_pre14(kk)>5 && V_HLB_single_pre14(kk)<=5.5
% %         bin(11)=bin(11)+1;
% %         transbin{11}(kk)=transcription_single_pre14(kk);
% %     else
% %         bin(12)=bin(12)+1;
% %         transbin{12}(kk)=transcription_single_pre14(kk);
% %     end
%     
%     if V_HLB_single_pre14(kk)<=1
%         bin(1)=bin(1)+1;
%         transbin{1}(kk)=transcription_single_pre14(kk);
%     elseif V_HLB_single_pre14(kk)>1 && V_HLB_single_pre14(kk)<=2
%         bin(2)=bin(2)+1;
%         transbin{2}(kk)=transcription_single_pre14(kk);
%     elseif V_HLB_single_pre14(kk)>2 && V_HLB_single_pre14(kk)<=3
%         bin(3)=bin(3)+1;
%         transbin{3}(kk)=transcription_single_pre14(kk);
%     elseif V_HLB_single_pre14(kk)>3 && V_HLB_single_pre14(kk)<=5
%         bin(4)=bin(4)+1;
%         transbin{4}(kk)=transcription_single_pre14(kk);
%     end
% 
% end
% 
% for kkk=1:12
%     transbin{kkk}(transbin{kkk}==0)=[];
% end
% 
% clear xlim
% % figure; clf
% % hold on
% % for kkk=1:9
% %     subplot(3,4,kkk)
% %     hist(transbin{kkk}, 20)
% % %     xlim([0, 15])
% % end
%         
% 
% %%
% k13=0;
% for k=1:length(files)
%     if cycle_comp{k}==13
%         k13=k13+1;
%         transcription_13(k) = nanmean(transcription_set{k});
%     else
%         transcription_13(k) = 0;
%     end
% end


%% Plot by the cell cycle
% color{11} = [1 0 0]; color{12} = [1 0.7 0]; color{13} = 'g'; color{14} = 'b';
% figure; clf
% hold on
% for k=1:length(files)
%     for p=1:length(F_green_pernuc_comp{k})
%         if ~isnan(F_green_pernuc_comp{k}(p)) && ~isnan(V_red_pernuc_comp{k}(p))
% %             plot(V_red_pernuc_comp{k}(p), F_green_pernuc_comp{k}(p), 'linestyle', 'none', 'color', color{cycle_comp{k}}, 'marker', '.', 'markersize', 15)
%             plot(F_green_pernuc_comp{k}(p), V_red_pernuc_comp{k}(p), 'linestyle', 'none', 'color', color{cycle_comp{k}}, 'marker', '.', 'markersize', 15)
%         end
%     end
% end
% 
% q=0;
% for k=1:length(files)
%     for p=1:length(F_green_pernuc_comp{k})
%         if ~isnan(F_green_pernuc_comp{k}(p)) && ~isnan(V_red_pernuc_comp{k}(p))
%             q=q+1;
%             V_HLB_single(q)=V_red_pernuc_comp{k}(p);
%             transcription_single(q)=F_green_pernuc_comp{k}(p);
%         end
%     end
% end
% 
% % [f,gof]=fit((V_HLB_single)',(transcription_single)','poly1');
% [f,gof]=fit((transcription_single)',(V_HLB_single)','poly1');
% R2 = gof.rsquare;
% koff = -f.p1;
% const = f.p2;
% x=linspace(0,15);
% 
% plot(x, f.p1*x+f.p2, 'k-')
% xlabel('Transcription (A.U.)'), ylabel('V_HLB (\mum^3)')