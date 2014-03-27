function [out] = fasim()
% R1 and R2* values of the MnCl2 solutions
R1 = [0.7645 0.9238 1.0593 1.2012 1.3369 1.4837 1.6556 1.8018 1.9342];
R2s = [4.8077 6.2893 7.8740 9.1743 10.8696 12.3457 13.8889 15.3846 16.6667];
% R2s = repmat(4.8077, 1, 9);

figure;
da = [0, 0.1:0.2:2];
N_da = length(da);
Col = hsv(N_da);
out = zeros(N_da, 4);
for idx = 1:N_da
    % Generate GRE signal with two flip-angles that are 
    % slighly off the desired values.
    [S1, TR] = signal(R1, R2s, 2+da(idx));
    [S2] = signal(R1, R2s, 12+da(idx));
    S_frac = S1./S2;
    R1r_ref = recon(S_frac, [2 12], TR);
    d = da(idx)/2;
    a = -TR*R1;
    x = 2/180*pi;
%    R1r = 1/TR * (log(( (1-cos(x)*exp(a)) ./ (1-exp(a)) * (1+d) - cos(x) )) - log(( (1-cos(x)*exp(a)) ./ (1-exp(a)) * (1+d) - 1 )));
%     R1r = -a/(d+1)/TR;
%    R1r = R1 + 1/TR * ((d*(exp(a)*cosd(12)-1) - d*exp(-a).*(exp(a)*cosd(12)-1))/(cosd(12)-1));
%    R1r = 1/TR*(log((1-cosd(12)*exp(a))./(1-exp(a))*(1+d)-cosd(12)) - log((1-cosd(12)*exp(a))./(1-exp(a))*(1+d)-1));
%    R1r = R1 - log(-d*cosd(12)-cosd(12)+d+1)/TR + log(-2*cosd(12)+d+1)/TR + R1*cosd(12)/(-2*cosd(12)+d+1);
%    R1r = 1/TR * log((-d*exp(-a)*cosd(12)-cosd(12)+d+1)./(d+exp(-a)*(-d*cosd(12)-cosd(12)+1)));
    R1r = R1 + 1/TR * log((1+d*exp(a))/(1+d));


    % Generate line of equality plot
    hold on;
    scatter(R1, R1r_ref, 20, 'k', 'filled');
    scatter(R1, R1r, 20, Col(idx, :), 'filled');
    P = robustfit(R1, R1r);
    minX = min([R1]);
    maxX = max([R1]);
    plot([minX maxX], polyval([P(2) P(1)], [minX maxX]), 'color', Col(idx, :));
    text(R1(end), R1r(end), sprintf('%0.1f%%', da(idx)));
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

function [R1] = recon(S_frac, fa, TR)
R1 = 1/TR.*log( (S_frac.*sind(fa(2)).*cosd(fa(1)) - sind(fa(1)).*cosd(fa(2))) ./ ...
                (S_frac.*sind(fa(2)) - sind(fa(1))) );
