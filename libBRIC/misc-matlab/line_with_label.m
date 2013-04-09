% Internal function

function [h] = line_with_label(val, lab, col, data_max, f_newlim, f_newpos, orient)
Lim = axis;
switch orient
    case 'v'
        if ~isempty(f_newlim)
            NewLim = Lim;
            NewLim(4) = (data_max-Lim(3))*f_newlim + Lim(3);
            axis(NewLim);
        end
        yposition = (data_max-Lim(3))*f_newpos + Lim(3);
        h = vline_ypos(val, yposition, lab, col);
    case 'h'
        if ~isempty(f_newlim)
            NewLim = Lim;
            NewLim(2) = (data_max-Lim(1))*f_newlim + Lim(1);
            axis(NewLim);
        end
        xposition = (data_max-Lim(1))*f_newpos + Lim(1);
        h = hline_xpos(xposition, val, lab, col);
end


%-------------------------------------------------------------------------
function [h] = vline_ypos(xpos, ypos, txt, col)
y=get(gca,'ylim');
plot([xpos xpos],y,col,'LineWidth', 1);
h = findobj(gca, 'type', 'line'); h = h(1);
text(double(xpos)+50, double(ypos), sprintf([txt '%.1f'], double(xpos)), 'color', 'k');

%-------------------------------------------------------------------------
function [h] = hline_xpos(xpos, ypos, txt, col)
x=get(gca,'xlim');
plot(x, [ypos ypos],col,'LineWidth', 1);
h = findobj(gca, 'type', 'line'); h = h(1);
text(double(xpos), double(ypos)+50, sprintf([txt '%.1f'], double(ypos)), 'color', 'k');



