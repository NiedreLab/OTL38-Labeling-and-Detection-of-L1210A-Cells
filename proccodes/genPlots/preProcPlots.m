function preProcPlots(time, data, bg, data_bs, std_sig, std_proc, params)
%PREPROCPLOTS Generates the folowing plots. Will not generate figures for
% empty y-value arguments (data,...std_proc).

% - Original Data
% - Data Background
% - Smoothed, background subtracted data
% - Signal moving STD
% - Processed moving STD

% Plot original data
if ~isempty(data)
    figure();
    hh = plot(time, data);
    set(hh, {'color'}, {params.color}')
    title('Original Data')
    xlabel('Time (s)')
    ylabel(['Signal (' params(1).units ')'])
    legend(params.name, 'Interpreter', 'none')
end

% Plot background
if ~isempty(bg)
    figure()
    hh = plot(time, bg);
    set(hh, {'color'}, {params.color}')
    title('Background')
    xlabel('Time (s)')
    ylabel(['Signal (' params(1).units ')'])
    legend(params.name, 'Interpreter', 'none')
end

% Plot smoothed, background subtracted data
if ~isempty(data_bs)
    figure()
    hh = plot(time, data_bs);
    set(hh, {'color'}, {params.color}')
    title('Smoothed, Bkgnd Subtracted Data')
    xlabel('Time (s)')
    ylabel(['Signal (' params(1).units ')'])
    legend(params.name, 'Interpreter', 'none')
end

% Plot raw moving std
if ~isempty(std_sig)
    figure()
    hh = plot(time, std_sig);
    set(hh, {'color'}, {params.color}')
    title('Moving Standard Deviation')
    xlabel('Time (s)')
    ylabel(['Signal (' params(1).units ')'])
    legend(params.name, 'Interpreter', 'none')
end

% Plot processed moving std
if ~isempty(std_proc)
    figure()
    hh = plot(time, std_proc);
    set(hh, {'color'}, {params.color}')
    title('Processed Standard Deviation')
    xlabel('Time (s)')
    ylabel(['Signal (' params(1).units ')'])
    legend(params.name, 'Interpreter', 'none')
end

end

