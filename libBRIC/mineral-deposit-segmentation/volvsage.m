addpath ./NIFTI
addpath ./libBRIC
addpath ./LIBRA

close all; clear all

pathfile = 'subjects_asps_356.csv';
anotfile = '../../../mineral/annot.csv';
matfile = 'subjects_asps_356_volvsage_t1wr2s_q0_45_both.mat';


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
% N_sub = length(Subj{1});
% In = NaN(N_sub, 11);
% Paths = cell(N_sub, 1);
% Chi = SubjAnot{1};
% idx_in = 1;
% for idx_subj = 1:N_sub
%     Path_in = Subj{1}{idx_subj};
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
%         fname = 'T2swHypo_mask';
%         S_t2swhypo = load_series(fullfile(Path_in, fname), []);
% %         fname = 'T2swHypoT1wHypo_mask';
% %         S_tmp = load_series(fullfile(Path_in, fname), []);
% %         S_t2swhypo(logical(S_tmp)) = 0;
%         NII = load_series(fullfile(Path_in, fname), 0);
%         F = NII.hdr.dime.pixdim(2:4);
%         In(idx_in, 1) = get_volume(S_t2swhypo == 13 | S_t2swhypo == 52, F);
%         In(idx_in, 2) = get_volume(S_t2swhypo == 11 | S_t2swhypo == 50, F);
%         In(idx_in, 3) = get_volume(S_t2swhypo == 12 | S_t2swhypo == 51, F);
%         In(idx_in, 4) = get_volume(S_t2swhypo == 14 | S_t2swhypo == 55, F);
%         In(idx_in, 5) = get_volume(logical(S_t2swhypo), F);
%         In(idx_in, 6) = 0;
%         % Increment
%         idx_in = idx_in + 1;
%     end
% end
% 
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

% load('../../../mineral/subjects_asps_107_v2_0-95.mat');
% Hu = NaN(length(Ret), 5);
% for idx_sub = 1:length(Ret)
%     I_ntis_means = Ret(idx_sub).I_ntis_means{2};
%     Hu(idx_sub, 1) = In(idx_sub, 7);
%     for idx_plot = 1:4
%         Hu(idx_sub, idx_plot+1) = I_ntis_means(idx_plot, 1);
%     end
% end
% close all;
% figure;
% for idx_plot = 1:4
%     subplot(2,2,idx_plot);
%     scatter(Hu(:, 1), Hu(:, idx_plot+1)); hold on;
%     M = Hu(:, 1) > 20;
%     Ha = mcdregres(Hu(M, 1), Hu(M, idx_plot+1), 'plots',0);
%     P = [Ha.slope Ha.int];
%     X = [min(Hu(M, 1)) max(Hu(M, 1))];
%     plot(X, polyval(P, X), 'r', 'linewidth', 1);
%     title(sprintf('y = %0.2g x + %0.2g', P(1), P(2)));
% end

