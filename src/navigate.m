% Function that allows for figure navigation

function navigate(src,evt)
    % Callback to parse keypress event data to adjust margins
    % Create figure environment with 
    %   f = figure('KeyPressFcn', @(obj, evt) navigate(obj,evt) );
    %
    % While the figure is selected:
    % Up/Down arrow keys adjust axis width by 50%
    % Left/Righ arrow keys shift axis by 50% of it's width
    %
    axisHandle = get(src, 'Children');
    range = get(axisHandle, 'XLim');
    switch logical(true)
        case strcmp(evt.Key, 'uparrow') % Widen margins
            set(axisHandle, 'XLim', range + 0.25 * diff(range) .*[-1 1] )
        case strcmp(evt.Key, 'downarrow') % Narrow margins
            set(axisHandle, 'XLim', range + 0.25 * diff(range) .*[1 -1] )
        case strcmp(evt.Key, 'leftarrow') % Shift margins left
            set(axisHandle, 'XLim', range - 0.5 * diff(range) ) 
        case strcmp(evt.Key, 'rightarrow') % Shift margins right
            set(axisHandle, 'XLim', range +  0.5 * diff(range) ) 
    end
end
