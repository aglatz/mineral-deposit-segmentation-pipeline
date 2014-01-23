clear all; clc; close all

addpath ../../misc-matlab
addpath ../../mineral-deposit-segmentation/NIFTI/
addpath ../../mineral-deposit-segmentation/LIBRA/

b = 0.05; % in mm
[X, Y, Z] = meshgrid(-7+b:b:7, -6+b:b:6, -6+b:b:6);
dB = zeros(size(X));

Ryz = sqrt(Y.^2 + Z.^2);
R = sqrt(X.^2 + Ryz.^2);

a = 1.5; % Diameter in mm
xi = .5e-6; % susceptibility w/r water inside
xe = 0; % susceptibility w/r water outside
dx = xi - xe;
B0 = 1.5; % in T, at R->inf

% Field change
M = R < a;
dB(M) = xi/3 * B0;
M = R >= a;
dB(M) = (dx/3*(a./R(M)).^3.*(3*(X(M).^2)./(X(M).^2+Ryz(M).^2)-1)+xe/3)*B0;
dB = permute(dB, [3 1 2]);

% figure
% imagesc(dB); colorbar; axis xy; axis equal; title('\DeltaB in T');

NII = make_nii(single(dB*1e6), [b b b], [0, 0, 0], 16); % in uT
save_series_raw(NII, 'dB');

% Lamor frequency change
g = 42.58e6; % Hz/T
dF = g*dB; % Hz

NII = make_nii(single(dF), [b b b], [0, 0, 0], 16);
save_series_raw(NII, 'dF');

TE = 20e-3:5e-3:60e-3; % s
N_TE = length(TE);
N = [10 10 20];

Vloss = zeros(1, N_TE);
%dphi = zeros([size(dF) N_TE]);
Mxysum = zeros([ceil(size(dF)./N) N_TE]);

for idx_TE = 1:N_TE
    fprintf('Calculating phase change at %0.3fms ...\n', TE(idx_TE)*1e3);
    
    % Phase change
    dphi = -2*pi*TE(idx_TE)*dF;

    % figure
    % imagesc(dphi); colorbar; axis xy; axis equal; title('\phi@TE in rad');

    % Calculate Mxy of GE signal
    Mxy = exp(sqrt(-1) * dphi);
    for idx1 = 1:size(Mxysum, 1)
        for idx2 = 1:size(Mxysum, 2)
            for idx3 = 1:size(Mxysum, 3)
                Tmp = Mxy((idx1-1)*N(1)+1:(idx1)*N(1), ...
                          (idx2-1)*N(2)+1:(idx2)*N(2), ...
                          (idx3-1)*N(3)+1:(idx3)*N(3));
                Mxysum(idx1, idx2, idx3, idx_TE) = abs(sum(Tmp(:)))/prod(N);
            end
        end
    end
    
    SM = Mxysum(:, :, :, idx_TE) < 0.99;
    Vloss(idx_TE) = get_volume(SM, N*b);
end

% NII = make_nii(single(dphi), [b b b TE(2)-TE(1)], [0, 0, 0], 16);
% save_series_raw(NII, 'dphi');

NII = make_nii(single(Mxysum), [N*b TE(2)-TE(1)], [0, 0, 0], 16);
save_series_raw(NII, 'amxy');

figure;
TE_ms = TE*1e3;
plot(TE_ms, Vloss);
xlabel('\bf TE in ms');
ylabel('\bf Vloss in mm^3');
hold on;
% Tmp = mcdregres(TE_ms(:), Vloss(:));
% P = [Tmp.slope Tmp.int];
P = polyfit(TE_ms, Vloss, 1);
plot(TE_ms, polyval(P, TE_ms), 'r');
title(sprintf('Vloss = %0.2f TE + %0.2f', P(1), P(2)));



