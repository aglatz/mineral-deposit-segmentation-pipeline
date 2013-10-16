close all; clear all

bpath='/media/LP3TBdisk/Andreas_PhD/mineral-deposit-segmentation-pipeline/BRICpipe/asps114';
mname='asps114_gre6magni';
pname_tmp='asps114_gre6phase';
pname='asps114_gre6phase_new';

S_mag = single(load_series(fullfile(bpath, mname), [])); % au
system(sprintf('cd %s; fslmaths %s -mul 3.14 -div 4096 %s', bpath, pname_tmp, pname));
S_pha = single(load_series(fullfile(bpath, pname), [])); % rad

fd = fopen(fullfile(bpath, 'GRE.nii.gz.txt'), 'r');
Idx_trunc = fscanf(fd, '%d %d');
fclose(fd);
S_mag = S_mag(:, :, Idx_trunc(1):Idx_trunc(2), :);
S_pha = S_pha(:, :, Idx_trunc(1):Idx_trunc(2), :);

fd = fopen(fullfile(bpath, 'tmp_GRE_brain.nii.gz.txt'), 'r');
Idx_trunc = fscanf(fd, '%d %d');
fclose(fd);

N_echo = 6;
S_mag = S_mag(:, :, Idx_trunc(1):Idx_trunc(2), 1:N_echo);
S_pha = S_pha(:, :, Idx_trunc(1):Idx_trunc(2), 1:N_echo);

save_series(fullfile(bpath, mname), fullfile(bpath, 'GRE_mag'), single(S_mag), []);
save_series(fullfile(bpath, pname), fullfile(bpath, 'GRE_pha'), single(S_pha), []);

se = strel('disk', 3);
SM_brain = imdilate(imerode(S_mag(:, :, :, 1) > 100, se), se);

te = 4.92e-3; %s
dte = 4.92e-3; %s
T = te:dte:te+dte*(N_echo-1); %s

%S_dif_1 = zeros(size(SM_brain, 1), size(SM_brain, 2), size(SM_brain, 3), N_echo-1);
for idx=1:1 %N_echo-1
    S_dif_2 =	(S_mag(:, :, :, idx+1).*exp(sqrt(-1)*S_pha(:, :, :, idx+1))) ./ ...
                (S_mag(:, :, :, idx).*exp(sqrt(-1)*S_pha(:, :, :, idx)));
    S_dif_2(~SM_brain) = 0;
%    S_dif_1(:, :, :, idx) = exp(sqrt(-1)*angle(S_dif_2));
end
S_dif = S_dif_2; %prod(S_dif_1, 4).^(1/(N_echo-1));
save_series(fullfile(bpath, mname), fullfile(bpath, 'F0'), single(angle(S_dif)/(2*pi)*1e6), []);
% clear S_dif_1 S_dif_2

S_dif_1 = zeros(size(SM_brain, 1), size(SM_brain, 2), size(SM_brain, 3), 5*5*5);
Var = zeros(3, 5*5*5);
idx = 1;
for idx_z=-2:2
    for idx_y=-2:2
        for idx_x=-2:2
            X = get_idx(idx_x, size(SM_brain, 1));
            Y = get_idx(idx_y, size(SM_brain, 2));
            Z = get_idx(idx_z, size(SM_brain, 3));
            S_dif_1(:, :, :, idx) = angle(S_dif(X, Y, Z) ./ S_dif) / dte;
            Var(:, idx) = [idx_x idx_y idx_z]';
            idx = idx + 1;
        end
    end
end

% wa = 1/sum(Var(1, :).^2);
% wb = 1/sum(Var(2, :).^2);
% wc = 1/sum(Var(3, :).^2);
% wd = 1/(3*3*3);
Wx = zeros(size(SM_brain));
Wy = zeros(size(SM_brain));
Wz = zeros(size(SM_brain));
W0 = zeros(size(SM_brain));
for idx_z=1:size(S_mag, 3)
    fprintf('%d ... \n', idx_z);
    for idx_y=1:size(S_mag, 2)
        for idx_x=1:size(S_mag, 1)
            D = reshape(S_dif_1(idx_x, idx_y, idx_z, :), 1, 5*5*5);
            if ~sum(isnan(D)) %&& sum(D == 0) == 1
%                 Wx(idx_x, idx_y, idx_z) = wa*sum(D.*Var(1, :));
%                 Wy(idx_x, idx_y, idx_z) = wb*sum(D.*Var(2, :));
%                 Wz(idx_x, idx_y, idx_z) = wc*sum(D.*Var(3, :));
%                 W0(idx_x, idx_y, idx_z) = wd*sum(D);
                b = regress(D', Var');
                Wx(idx_x, idx_y, idx_z) = b(1);
                Wy(idx_x, idx_y, idx_z) = b(2);
                Wz(idx_x, idx_y, idx_z) = b(3);
            end
        end
    end
end
save_series(fullfile(bpath, mname), fullfile(bpath, 'Wx'), single(Wx), []);
save_series(fullfile(bpath, mname), fullfile(bpath, 'Wy'), single(Wy), []);
save_series(fullfile(bpath, mname), fullfile(bpath, 'Wz'), single(Wz), []);
% save_series(fullfile(bpath, mname), fullfile(bpath, 'W0'), single(W0), []);

R2s_1 = zeros(size(S_mag, 1), size(S_mag, 2), size(S_mag, 3), 2);
% S_mag_1 = zeros(size(S_mag), 'single');
R2s_2 = zeros(size(S_mag, 1), size(S_mag, 2), size(S_mag, 3), 2);
% S_mag_2 = zeros(size(S_mag), 'single');
% R2s_3 = zeros(size(S_mag, 1), size(S_mag, 2), size(S_mag, 3), 2);
S_tmp = S_mag(:, :, :, 1);
S_tmp(~SM_brain) = [];
I_thr = quantile(double(S_tmp), .05);
N = [3.1 3.1 3.1 2.7 3.1 2.7];
for idx_z=1:size(S_mag, 3)
    fprintf('%d ... \n', idx_z);
    for idx_y=1:size(S_mag, 2)
        for idx_x=1:size(S_mag, 1)
            if SM_brain(idx_x, idx_y, idx_z) && ...
               double(S_mag(idx_x, idx_y, idx_z, 1)) > I_thr
                % close all;
                % figure;
                Mag = double(reshape(S_mag(idx_x, idx_y, idx_z, 1:N_echo), 1, N_echo));
                %plot(T, Mag);
                %hold on;
                P = polyfit(T, log(Mag), 1);
                SNR = Mag./N;
                res = @(c) (exp(c(2)).*exp(c(1).*T') - Mag'); %.*(SNR'./max(SNR));
                if isnan(sum(res(P).^2))
                    continue;
                end
                [C, ssq] = LMFnlsq(res, P); %, 'Display', 1);
                R2s_1(idx_x, idx_y, idx_z, :) =  [-C(1), ssq];
                % S_mag_1(idx_x, idx_y, idx_z, :) = single(exp(polyval(C, T)));
                %plot(T, exp(polyval(C, T)), 'g');

                wz = double(Wz(idx_x, idx_y, idx_z));
                if wz == 0
                    F = ones(size(T));
                else
                    F = (sin(wz/2*T)./(wz/2*T));
                end
                res = @(c) (exp(c(2)).*exp(c(1).*T').*F' - Mag'); %.*(SNR'./max(SNR));
                if isnan(sum(res(P).^2))
                    continue;
                end
                [C, ssq] = LMFnlsq(res, P); %, 'Display', 1);
                R2s_2(idx_x, idx_y, idx_z, :) =  [-C(1), ssq];
                % S_mag_2(idx_x, idx_y, idx_z, :) = single(exp(polyval(C, T)));
                %plot(T, exp(polyval(C, T)), 'k');
                
%                 dB0_v = double(reshape(dB0(idx_x, idx_y, :), 1, size(dB0, 3)))';
%                 dx = 4; %mm
%                 x = (1:size(dB0, 3))'*dx-dx/2; %mm
%                 x0 = idx_z*dx-dx/2; %mm
%                 Res = zeros(size(T));
%                 %Res_0 = int_grad(x, dx, x0, dB0_v, 0);
%                 for idx = 1:length(T)
%                     Res(idx) = int_grad(x, dx, x0, dB0_v, T(idx)) ./ dx;
%                 end
%                 res = @(c) real(exp(c(2)).*exp(c(1).*T').*abs(Res)' - Mag');
%                 if isnan(sum(res(P).^2))
%                     continue;
%                 end
%                 [C, ssq] = LMFnlsq(res, P); %, 'Display', 1);
%                 R2s_3(idx_x, idx_y, idx_z, :) =  [-C(1), ssq];
            end
        end
    end
end
save_series(fullfile(bpath, mname), fullfile(bpath, 'R2s_1'), single(R2s_1), []);
% save_series(fullfile(bpath, mname), fullfile(bpath, 'S_mag_1'), single(S_mag_1), []);
save_series(fullfile(bpath, mname), fullfile(bpath, 'R2s_2'), single(R2s_2), []);
% save_series(fullfile(bpath, mname), fullfile(bpath, 'S_mag_2'), single(S_mag_2), []);
% save_series(fullfile(bpath, mname), fullfile(bpath, 'R2s_3'), single(R2s_3), []);

% figure;
% for hu=1:28
% dx = 4; %mm
% x = (1:size(dB0, 3))'*dx-dx/2; %mm
% x0 = hu*dx-dx/2; %mm
% gz = reshape(dB0(92, 150, :), size(dB0, 3), 1)*2; %rad/s
% Res = zeros(size(T));
% for idx = 1:length(T)
%     Res(idx) = exp(polyval(P, T(idx))).*int_grad(x, dx, x0, gz, T(idx));
% end
% plot(T, exp(polyval(P, T)), 'k'); 
% hold on
% plot(T, abs(Res)/dx, 'r');
% hold off;
% pause(1);
% end

% figure;
% R2s = 0;
% for f = [0 0.8 1.6 2.4 3.2 4]
%     x = (-10:10)'*10; %mm
%     dx = 10; %mm
%     x0 = 0; %mm
%     T = (0:1:50)*1e-3; %s
%     gz = f*x*2*3.1416; %rad/s
%     Res = zeros(size(T));
%     for idx = 1:length(T)
%         Res(idx) = exp(-R2s*T(idx))*int_grad(x, dx, x0, gz, T(idx));
%     end
%     plot(T, abs(Res)/dx); hold on
%     % plot(T, exp(-R2s*T).*(sin(gz/2*T)./(gz/2*T))); hold on;
% end




