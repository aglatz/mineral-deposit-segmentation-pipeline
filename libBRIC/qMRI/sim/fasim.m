function [out] = fasim()
% R1 and R2* values of the MnCl2 solutions
R1 = [0.7645 0.9238 1.0593 1.2012 1.3369 1.4837 1.6556 1.8018 1.9342];
R2s = [4.8077 6.2893 7.8740 9.1743 10.8696 12.3457 13.8889 15.3846 16.6667];

figure;
da = [1:0.02:1.2]; %[0, 0.1:0.2:2];
N_da = length(da);
Col = hsv(N_da);
out = zeros(N_da, 4);
for idx = 1:N_da
    % Generate GRE signal with two flip-angles that are 
    % slighly off the desired values.
    [S1, TR] = signal(R1, R2s, 2*da(idx));
    [S2] = signal(R1, R2s, 12*da(idx));

    % Recon with desired values
    R1r = recon([S1; S2], [2 12], TR);

    % Generate line of equality plot
    scatter(R1, R1r, 20, Col(idx, :), 'filled');
    P = robustfit(R1, R1r);
    hold on;
    minX = min([R1]);
    maxX = max([R1]);
    plot([minX maxX], polyval([P(2) P(1)], [minX maxX]), 'color', Col(idx, :));
    text(R1(end), R1r(end), sprintf('%0.1f%%', da(idx)*100-100));
    R1h = P(1) + P(2)*R1;
    out(idx, 1) = da(idx);
    out(idx, 2) = 1-norm(R1r-R1h)/norm(R1r-mean(R1h));
    out(idx, [3 4]) = P;
end
plot([minX maxX], [minX maxX], '--k');
xlabel('\bf Reference R1_{ref} in s^{-1}');
ylabel('\bf Reconstructed R1 in s^{-1}');
%text(minX, maxX/2, sprintf('R1=%0.2fR1_{ref}+(%0.2f)', P(2), P(1)));
set(gcf, 'color', 'white')



function [S, TR] = signal(R1, R2s, fa)
M0 = 1000; %au
TR = 8e-3; %ms
TE = 3.5e-3; %ms

S = M0.*sind(fa).*(1-exp(-TR.*R1)).*exp(-TE.*R2s)./(1-cosd(fa).*exp(-TR.*R1));

function [R1] = recon(S, fa, TR)
R1 = 1/TR.*log( (S(1,:).*sind(fa(2)).*cosd(fa(1)) - S(2,:).*sind(fa(1)).*cosd(fa(2))) ./ ...
                (S(1,:).*sind(fa(2)) - S(2,:).*sind(fa(1))) );