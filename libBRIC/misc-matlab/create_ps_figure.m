function [handle] = create_ps_figure()
% Creates an 'invisible' figure.
% INPUTS: None
% RETURNS: the figure handle
%
handle = figure('Visible', 'off');
