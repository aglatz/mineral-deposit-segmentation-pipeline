function [] = save_ps_figure(basename, handle)
% Function for saving a figure as a Postscript graphic file.
% No valid handle means we should delete the file (clean-up). This
% has the advantage that just this function has to know about the
% specific file name.
% INPUTS: basename - filepath and basename of the postscript file
%         handle - filehandle of figure
%
[path, name, ext] = fileparts(basename);
if ~strcmp(ext, '.ps')
	ext = '.ps';
end
basename = [path '/' name];
if handle
    %set(handle, 'PaperType', 'a4');
    set(handle, 'PaperUnits', 'centimeters');
    %set(handle, 'PaperPosition', [0.5 0.5 27 20]); % landscape
    set(handle, 'PaperPosition', [0.5 0.5 20 27]); % portrait
    print(handle, [basename ext], '-dpsc2', '-append', '-r300');
	close(handle);
else
    delete([basename ext]);
end
