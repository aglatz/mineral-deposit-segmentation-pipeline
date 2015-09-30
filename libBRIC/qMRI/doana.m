function [] = doana()

bdir = '/home/s1063233/roiana';
% sdirs = {'23499', '23578', '23634', '23706', '23738', '23788', '23899', '24174'};
sdirs = {'23499', '23634', '23706', '23738', '23788', '23899', '24174'};
TSs = {'2015/03/09', '2015/04/14', '2015/05/05', '2015/05/13', '2015/05/22', '2015/06/23', '2015/09/08'};
TSformat = 'yyyy/mm/dd';
N_rois = 9;
c = (0.05:0.02:0.22)';
Type = 'R2s';

R1ref = zeros(N_rois, length(sdirs));
for idx = 1:length(sdirs)
    fname = fullfile(bdir, sdirs{idx}, 'T2s_ref', 'R2s_roiana');
    load(fname);
    R1ref(:,idx)=ROIs(:,1);
end


R1 = zeros(N_rois, length(sdirs));
for idx = 1:length(sdirs)
    fname = fullfile(bdir, sdirs{idx}, 'T2s', 'R2s_roiana');
    load(fname);
    R1(:,idx)=ROIs(:,1);
end

plot_R_vs_Rref(R1, R1ref, Type);
plot_R_vs_c(R1, c, Type);
plot_R_vs_c(R1ref, c, [Type 'ref']);
plot_BA(R1, R1ref, Type);
plot_R_vs_t(R1, datenum(TSs, TSformat), Type);
plot_R_vs_t(R1ref, datenum(TSs, TSformat), [Type 'ref']);


% ------- plot_R_vs_Rref -------
function [] = plot_R_vs_Rref(R, Rref, Type)
Rref_mean = median(Rref, 2);
% Rref_std = std(Rref, 0, 2);
R_mean = median(R, 2);
% R_std = std(R, 0, 2);
% L = 1:size(R, 1);

figure; hold on;
scatter(Rref(:), R(:), 20, 'fill');
% for idx=L
%     plot([Rref_mean(idx)-Rref_std(idx) Rref_mean(idx)+Rref_std(idx)], [R_mean(idx) R_mean(idx)]);
%     plot([Rref_mean(idx) Rref_mean(idx)], [R_mean(idx)-R_std(idx) R_mean(idx)+R_std(idx)]);
% end
I_min_max = [min([R_mean; Rref_mean]) max([R_mean; Rref_mean])];
% text(Rref_mean, R_mean, num2str(L'));
plot(I_min_max, I_min_max, '--k'); % line of equality
xlabel(['\bf ' Type '_{ref} in s^{-1}']);
ylabel(['\bf ' Type ' in s^{-1}']);
P = robustfit(Rref(:), R(:));
plot(I_min_max, polyval([P(2) P(1)], I_min_max), '--b');
if P(1) < 0
    Tmp = '-';
else
    Tmp = '+';
end 
text(I_min_max(1), I_min_max(2), ...
    sprintf(['\\bf ' Type '=%0.2f ' Type '_{ref} %s %0.2f'], P(2), Tmp, abs(P(1))));
set(gcf, 'color', 'w');
axis([I_min_max, I_min_max]);

P = robustfit(R_mean, Rref_mean);
if P(1) < 0
    Tmp = '-';
else
    Tmp = '+';
end 
text(I_min_max(1), I_min_max(2)*.9, ...
    sprintf(['\\bf ' Type '_{cal}=%0.2f ' Type ' %s %0.2f'], P(2), Tmp, abs(P(1))));
set(gcf, 'color', 'w');
axis([I_min_max, I_min_max]);


% ------- plot_R_vs_c --------
function [] = plot_R_vs_c(R, c, Type)
figure; hold on
xlabel('\bf MnCl_2 concentration c in mMol');
ylabel(['\bf Relaxation rate ' Type ' in s^{-1}']);
set(gcf, 'color', 'w');

C = repmat(c, 1, size(R, 2));
scatter(C(:), R(:), 20, 'fill');
P = robustfit(C(:), R(:));
plot(c, polyval([P(2) P(1)], c), '--k');
if P(1) < 0
    Tmp = '-';
else
    Tmp = '+';
end
R_mean = mean(R, 2);
text(min(c), max(R_mean), ...
    sprintf(['\\bf ' Type '=%0.2f s^{-1}/mmol/l c %s %0.2f s^{-1}; r=%0.3f'], ...
            P(2), Tmp, abs(P(1)), corr(R_mean, c)));
set(gcf, 'color', 'w');


% ---- plot_BA -------
function [] = plot_BA(R, Rref, Type)
D = (R-Rref).*2./(R+Rref).*100;
A = (R+Rref)/2;

figure; hold on;
scatter(A(:), D(:), 20, 'fill');
P = robustfit(A(:), D(:));
hold on;
Tmp = [min(A(:)) max(A(:))];
plot(Tmp, polyval([P(2) P(1)], Tmp), '--k');
xlabel('\bf Average MnCl_2 concentration c in mMol');
ylabel(['\bf Rel. difference in relaxation rate of ' Type ' in %']);
set(gcf, 'color', 'w');


% --- plot_R_vs_t ---
function [] = plot_R_vs_t(R, t, Type)
t = t-min(t);
Tminmax = [min(t) max(t)];
T=repmat(t', size(R, 1), 1);


figure; hold on;
scatter(T(:), R(:), 20, 'fill');

for idx=1:size(R, 1)
    tbl = table(t, R(idx, :)', 'VariableNames', {'Time', 'Rate'});
    lm = fitlm(tbl, 'Rate~Time');
    P = lm.Coefficients.Estimate;
    %P = robustfit(t', R(idx, :));
    Rest = polyval([P(2) P(1)], Tminmax);
    plot(Tminmax,Rest, '--k');
    text(Tminmax(1), Rest(1), sprintf('p(I)=%0.2f p(k)=%0.2f', lm.Coefficients.pValue(1), lm.Coefficients.pValue(2)));
end

xlabel('\bf Acquisition time in days');
ylabel(['\bf Relaxation rate ' Type ' in s^{-1}']);
set(gcf, 'color', 'w');

