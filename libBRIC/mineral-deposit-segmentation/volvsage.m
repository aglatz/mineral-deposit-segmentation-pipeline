addpath ./NIFTI
addpath ./libBRIC
addpath ./LIBRA

close all; clear all

basepath = '/media/LP3TBdisk/Andreas_PhD/mineral-deposit-segmentation-pipeline/BRICpipe/';
pathfile = [basepath 'subjects_asps_318_checked_both_t2w.csv'];
pipefile = [basepath 'subjects_asps_371.csv'];
outfile = [basepath 'subjects_asps_371_lt65.csv'];
anotfile = '/media/LP3TBdisk/Andreas_PhD/mineral/annot.csv';
matfile = 'subjects_asps_318_volvsage_t1wgre_2_q0_6_both_t2w_all.mat';


%% Uncomment below for the first time
% %
% fd = fopen(pathfile);
% Subj = textscan(fd, '%s%*[^\n]', ...
%     'delimiter',',', ...
%     'treatAsEmpty',{'NA','na'}, ...
%     'commentStyle', '#');
% fclose(fd);
% %
% fd = fopen(anotfile);
% SubjAnot = textscan(fd, '%s%s%s%s%s%s%*[^\n]', ...
%     'delimiter',',',...
%     'treatAsEmpty',{'NA','na'}, ...
%     'commentStyle', '#');
% fclose(fd);
% %
% fd = fopen(pipefile);
% SubjPipe = textscan(fd, '%s%s%s%*[^\n]', ...
%     'delimiter',',',...
%     'treatAsEmpty',{'NA','na'}, ...
%     'commentStyle', '#');
% fclose(fd);
% %
% % fd = fopen(outfile, 'w');
% %
% N_sub = length(Subj{1});
% In = NaN(N_sub, 11+8+2);
% Paths = cell(N_sub, 1);
% Chi = SubjAnot{1};
% Chi_pipe = SubjPipe{1};
% idx_in = 1;
% for idx_subj = 1:N_sub
%     Path_in = [basepath Subj{1}{idx_subj}];
%     % Find corresponding annotation (chi number that matches path)
%     f=@(in) (~isempty(strfind(Path_in, in)));
%     M_subj = cellfun(f, Chi);
%     if sum(M_subj) == 1
%         In(idx_in, 7) = str2double(SubjAnot{2}{M_subj}); % Age
%         In(idx_in, 8) = str2double(SubjAnot{3}{M_subj}); % mmse
%         In(idx_in, 9) = str2double(SubjAnot{4}{M_subj}); % WMLvol
%         In(idx_in, 10) = str2double(SubjAnot{5}{M_subj}); % Brainvol
%         if strcmp(SubjAnot{6}{M_subj}, 'm')
%             In(idx_in, 11) = 0;
%         end
%         if strcmp(SubjAnot{6}{M_subj}, 'w')
%             In(idx_in, 11) = 1;
%         end
%         Paths{idx_in} = Path_in;
%         Path_in
%         SubjAnot{1}{M_subj}
%         %
%         fname = 'T2swHypoT1wiso_mask';
%         S_t2swhypo = load_series(fullfile(Path_in, fname), []);
% %         fname = 'T2swHypoT1whypo_mask';
% %         S_tmp = load_series(fullfile(Path_in, fname), []);
% %         S_t2swhypo(logical(S_tmp)) = 0;
% %         fname = 'T2swHypoT1whyper_mask';
% %         S_tmp = load_series(fullfile(Path_in, fname), []);
% %         S_t2swhypo(logical(S_tmp)) = 0;
% 
%         NII = load_series(fullfile(Path_in, fname), 0);
%         F = NII.hdr.dime.pixdim(2:4);
%         In(idx_in, 1) = get_volume(S_t2swhypo == 13 | S_t2swhypo == 52, F);
%         In(idx_in, 2) = get_volume(S_t2swhypo == 11 | S_t2swhypo == 50, F);
%         In(idx_in, 3) = get_volume(S_t2swhypo == 12 | S_t2swhypo == 51, F);
%         In(idx_in, 4) = get_volume(S_t2swhypo == 14 | S_t2swhypo == 55, F);
%         In(idx_in, 5) = get_volume(logical(S_t2swhypo), F);
%         In(idx_in, 6) = sum(In(idx_in, [1 3 4]));
%         %
% %         S_r2s = load_series(fullfile(Path_in, 'R2s'), []);
% %         S_ntis = load_series(fullfile(Path_in, 'NormTis_mask'), []);
% %         SM = S_ntis == 13 | S_ntis == 52;
% %         In(idx_in, 12) = median(S_r2s(SM));
% %         SM = S_ntis == 11 | S_ntis == 50;
% %         In(idx_in, 13) = median(S_r2s(SM));
% %         SM = S_ntis == 12 | S_ntis == 51;
% %         In(idx_in, 14) = median(S_r2s(SM));
% %         SM = S_ntis == 14 | S_ntis == 55;
% %         In(idx_in, 15) = median(S_r2s(SM));
% %         SM = S_ntis == 13 | S_ntis == 52 | S_t2swhypo == 13 | S_t2swhypo == 52;
% %         In(idx_in, 16) = median(S_r2s(SM));
% %         SM = S_ntis == 11 | S_ntis == 50 | S_t2swhypo == 11 | S_t2swhypo == 50;
% %         In(idx_in, 17) = median(S_r2s(SM));
% %         SM = S_ntis == 12 | S_ntis == 51 | S_t2swhypo == 12 | S_t2swhypo == 51;
% %         In(idx_in, 18) = median(S_r2s(SM));
% %         SM = S_ntis == 14 | S_ntis == 55 | S_t2swhypo == 13 | S_t2swhypo == 55;
% %         In(idx_in, 19) = median(S_r2s(SM));
%         %
% %         In(idx_in, 20) = get_volume(logical(S_t2swhypo), F);
% %         fname = 'T2swHypo_bin_notreg_mni_mask';
% %         SM_t2swhypomni = logical(load_series(fullfile(Path_in, fname), []));
% %         NII = load_series(fullfile(Path_in, fname), 0);
% %         F = NII.hdr.dime.pixdim(2:4);        % Increment
% %         In(idx_in, 21) = get_volume(SM_t2swhypomni, F);
%         idx_in = idx_in + 1;
%     end
% % 	% Find corresponding annotation (chi number that matches path)
% %     f=@(in) (~isempty(strfind(Path_in, in)));
% %     M_pipe = cellfun(f, Chi_pipe);
% %     if sum(M_subj) == 1
% %         Age = str2double(SubjAnot{2}{M_subj});
% %         if sum(M_pipe) == 1 && Age < 65 && Age > 20
% %             ID = SubjPipe{1}{M_pipe};
% %             T1W = SubjPipe{2}{M_pipe};
% %             GRE = SubjPipe{3}{M_pipe};
% %             fprintf(fd, '%s, %s, %s,\n', ID, T1W, GRE);
% %         end
% %     end
% end
% %
% % fclose(fd);
% %
% In = In(1:idx_in-1, :);
% Paths = Paths(1:idx_in-1);
% 
% save(matfile, 'In', 'Paths');
%% Uncomment above for the first time

load(matfile);

idx_sel = 5; % what should be on the y-axis? 1: GP vol, 2: nGP vol, 3: total BG vol
Idx = [idx_sel 7 10];
M  = sum(In(:, Idx)>0, 2) == length(Idx); sum(M) % Exclude BGID volumes == 0
Paths = Paths(M);  %Subject
X = In(M, 7); % age
Y = In(M, idx_sel)./In(M, 10)*1e6; % scale per ICV

figure;
[idx_oli, p_med] = plot_quantreg(X, Y);
Out = [Paths(idx_oli), num2cell(X(idx_oli)), num2cell(Y(idx_oli))];
save_xls('volvsage_OLI', {'%s', '%d', '%0.2f'}, Out);
xlabel('\bf Age in years');
ylabel('\bf BGID volume per ICV in ppm');
%set(gca, 'YScale', 'log');
[rho, p] = corr(X, Y, 'type', 'spearman');
title(sprintf('N=%d, k=%0.3g, d=%0.3g; rho=%0.3g, p=%0.3g', sum(M), p_med(1), p_med(2), rho, p));

% M = sum(isnan(In),2) == 0 & In(:, 5) >= 65;
% Mat = In(M, [7, 8, 9, 10]);
% ha = mcdcov(Mat, 'plots', 1, 'alpha', 0.75);
% [~, Idx] = sort(ha.rd);
% save_xls('paths_sorted', {'%s', '%0.3f', '%d', '%0.2f', '%0.2f', '%0.2f'}, [Paths(Idx), num2cell([ha.rd(Idx)' Mat(Idx, :)])]);
% 
% for idx_type = [0 4]
%     figure;
%     for idx_plot = 1:3
%         subplot(2,2,idx_plot);
%         %scatter(X, In(M, 11+idx_type+idx_plot)); hold on;
%         [idx_oli, p_med] = plot_quantreg(X, In(M, 11+idx_type+idx_plot)); hold on;
%         switch idx_plot
%             case 1, R2=31.2*(1-exp(-0.069*X))+10.0; ti = 'GP';
%             case 2, R2=18.3*(1-exp(-0.023*X))+12.0; ti = 'C';
%             case 3, R2=36.4*(1-exp(-0.013*X))+12.8; ti = 'P';
%             case 4, R2=36.4*(1-exp(-0.013*X))+12.8; ti = 'GP';
%         end
%         [X_s, I_s] = sort(X);
%         plot(X_s, R2(I_s), 'g');
% %         ha = mcdregres(X, In(M, 11+idx_type+idx_plot), 'plots', 0);
% %         P = [ha.slope ha.int];
% %         Tmp = [min(X) max(X)];
% %         plot(Tmp, polyval(P, Tmp), 'r', 'linewidth', 1);
%         title(sprintf('%s: y = %0.2g x + %0.2g', ti, p_med(1), p_med(2)));
%     end
% end
% figure;
% blandAltmanPlot(In(M, 21), In(M, 20));

