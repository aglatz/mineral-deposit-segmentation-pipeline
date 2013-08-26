% Creates a stacked bar plot from the values in Y1 and Y2.
% It's a workaround since transparent figures cannot be
% saved in eps figures.
% INPUTS: X - the bin centres
%         Y1, Y2 - Occurences
%         ti1, ti2 - titles corresponding to Y1 and Y2 for the legend
function plot_bar2(X, Y1, Y2, ti1, ti2)
Diff1 = (Y1(:)-Y2(:));
Diff1(Diff1 <= 0) = 0;

Diff2 = (Y2(:)-Y1(:));
Diff2(Diff2 <= 0) = 0;

Min12 = min([Y1(:) Y2(:)], [], 2);

Y = [Min12 Diff1, Diff2];
bar(X, Y, 0.8, 'stack');
% legend([ti1 '+' ti2], ti1, ti2);
    