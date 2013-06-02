function Idx = blandAltmanPlot(A,B, varargin)
%reference: Y H Chan, Biostatistics 104:
%Correlational Analysis,
%Singapore Med J 2003 Vol 44(12) : 614-619
meanAB=(A+B)./2;
if ~isempty(varargin)
    difff = varargin{1};
else
    difff = (A-B)./meanAB;
end
M = isnan(difff);
meanAB(M)=[];
difff(M)=[];
%Q = quantile(difff, [.05 .5 .95]);
meanDiff=mloclogist(difff);
stdDiff=mscalelogist(difff);

f = 1.95;
meanp2D=meanDiff+f*stdDiff;
meanm2D=meanDiff-f*stdDiff;
n=length(difff);
minD=min(meanAB)-0.1;
maxD=max(meanAB)+0.1;

figure;
scatter(meanAB, difff, 20, 'k');
%V_low=quantile(B(~M), .25);
M_low=B(~M)<40;
hold on;
M = abs(difff-meanDiff) > f*stdDiff;
Idx = find(M);
for idx = 1:sum(M)
    text(meanAB(Idx(idx)), difff(Idx(idx)), num2str(Idx(idx)));
end
scatter(meanAB(M_low), difff(M_low), 20, 'r');
%text(meanAB,difff,num2str((B(~M))));
hold on;
X=[minD maxD];
plot(X,ones(1,2)*meanp2D,'--k');
text(maxD*.75,meanp2D+0.01,sprintf('mu+%0.2f*sd=%0.3f', f, meanp2D));
hold on;
plot(X,ones(1,2)*meanm2D,'--k');
text(maxD*.75,meanm2D+0.01,sprintf('mu-%0.2f*sd=%0.3f', f, meanm2D));
hold on;
plot(X,ones(1,2)*meanDiff,'k');
text(maxD*.75,meanDiff,sprintf('mu=%0.3f', meanDiff));
Tmp = mcdregres(meanAB(:), difff(:), 'plots', 0);
Est_param1 = [Tmp.slope Tmp.int];
plot(X, polyval(Est_param1, X), '--b');
text(maxD*.5,polyval(Est_param1, maxD*.5)+0.01, sprintf('y=%fx+%f', Est_param1(1), Est_param1(2)));
% if ~isempty(C)
%     scatter(meanAB(:), C(:), 20, 'g');
%     Tmp = mcdregres(meanAB(:), C(:), 'plots', 0);
%     Est_param2 = [Tmp.slope Tmp.int];
%     plot(X, polyval(Est_param2, X), '--g');
%     text(maxD*.5,polyval(Est_param2, maxD*.5)+0.01, sprintf('y=%fx+%f', Est_param2(1), Est_param2(2)));
%     k = Est_param2(1)/Est_param1(1)
%     d = Est_param2(2)-k*Est_param1(2)
% end
xlim(X);
xlabel('Average volume in mm^3');
ylabel('(Vol-Vol_{ref}) per average volume');
legend('Vol_{ref}>=40mm^3','Vol_{ref}<40mm^3')