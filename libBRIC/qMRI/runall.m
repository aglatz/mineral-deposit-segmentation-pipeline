clear all; close all

addpath ~/mineral-deposit-segmentation-pipeline/libBRIC/qMRI
addpath ~/mineral-deposit-segmentation-pipeline/libBRIC/qMRI/LMFnlsq/
addpath ~/mineral-deposit-segmentation-pipeline/libBRIC/misc-matlab/
addpath /ISIS/proc1/aglatz/mineral/NIFTI/

cdir = '/home/s1063233/23899/T2s_ref';
fname = fullfile(cdir, 'R2s');
Nz = 1;
% R2s = dcmr2sfit(cdir);
%R2s = 1000./double(rot90(load_series(fname, []), 1));
R2s = double(rot90(load_series(fname, []), 1));
ismrpha=true;
I_lims = [0 2];
figure; imagesc(R2s(:,:,Nz), I_lims);
axis image
axis off
colormap('gray');
colorbar;
set(gcf, 'color', 'w');
% export_fig(['R2s_mag.pdf'], '-a1',  '-q101');

% load('../T2s_2/mask.mat');
S_roi = rot90(load_series([fname '_roi'], []), 1);

% for erroding the mask
% SE = strel('disk', 2);
% SM = imerode(logical(S_roi), SE);
% S_roi(~SM) = 0;

% for inverting the lables
% for idx=1:9
%     SM = S_roi == idx;
%     S_roi(SM) = 9-idx+11;
% end
% SM = S_roi > 10;
% S_roi(SM) = S_roi(SM) - 10;

[ROIs] = roiana(R2s, S_roi);
save([fname '_roiana'], 'ROIs');

Tmp = R2s;
Tmp(logical(S_roi)) = 0;
figure; imagesc(Tmp(:,:,Nz), I_lims);
axis image
axis off
colormap('gray');
colorbar;
set(gcf, 'color', 'w');
% export_fig(['R2s_mask.pdf'], '-a1',  '-q101');

if ismrpha
    c = (0.05:0.02:0.22)';
    %c = c(end:-1:1)';
    figure; errorbar(c, ROIs(:,1), ROIs(:,2));
    Type = '2';
    xlabel('\bf MnCl_2 concentration c in mMol');
    ylabel(['\bf Mean ROI relaxivitiy R_{' Type '} in s^{-1}']);
    set(gcf, 'color', 'w');
    save(fullfile(cdir, 'R2s.mat'), 'R2s', 'c', 'ROIs');

    P = robustfit(c, ROIs(:,1));
    hold on;
    plot(c, polyval([P(2) P(1)], c), '--k');
    if P(1) < 0
        Tmp = '-';
    else
        Tmp = '+';
    end
    text(min(c), max(ROIs(:,1)), ...
        sprintf(['\\bf R_{' Type '}=%0.2f s^{-1}/mmol/l c %s %0.2f s^{-1}; r=%0.3f'], ...
                P(2), Tmp, abs(P(1)), corr(ROIs(:,1), c)));
    set(gcf, 'color', 'w');
    % export_fig(['R2s_vs_c.pdf'], '-a1',  '-q101');
end
