 function Annotation_Matrix = MEP_Visualizer(data, config)
 
% A simple tool to visualise and inspect biosignals, e.g. Motor-evoked
% potentials (MEP). press 'h' to display key bindings.
% 
% 2021 - Maximilian Hommelsen
%
% - Press 'h' to show all key bindings.
% - Use 'a' and 'd' to navigate through the trials
% - Use 'w' and 's' to annotate a trial as good and bad, respectively
% - Keys 1, 2, and 3 can also be used for annotation; 1 = good, 2 = bad, 3 = other.
% - Peak-to-peak amplitude can be extracted based on a time window (config.PeakWindow)
% - Use 'f' to modify a Peak Window for an individual trial.
% - Press 'p' to save annotation_file


%%% Basic Usage:

%    annotation_matrix = MEP_Visualizer(data)
% or
%    annotation_matrix = MEP_Visualizer(data, config)   
% 
% 'data' (matrix) should either have 2d (Samples x Trials) or 3d (Samples x Trials x
% Channels)
% 
% 'config' (struct) contains additional parameters that can be provided/modified by the user.
% 
%   config.SamplingFrequency = 10000; % Sampling Fequency in Hz
%   config.TriggerTime       = 0; % Occurence of Pulse/Trigger; Time in [S] relative to the first sample; If provided, this is used as reference point for the PeakWindow
%   config.PeakWindow        = [0.02 0.04] % Time Window for calculation of peak to peak amplitude; >> Time in [S] relative to Trigger Time<<
%   config.ChannelNames      = {'Chan 1', 'Chan 2'} % Cell array with Channel Names. Size needs to be == size(data,3)
%   config.demean            = 'no' % Demean each single trial, 'yes' or 'no'
%   config.unit              = 'V' % Only for labeling of axis/text
%   config.Xlimit            = Set specific X-Axis Limits;
%   config.Ylimit            = Set specific Y-Axis Limits;
%
% - Returns a simple matrix of shape(Trials, 4, Channels), where the four values describe
%   Annotation Type (=1, 2, 3), Onset of PeakWindow, Offset of PeakWindow, Peak-to-peak Amplitude


     % Close Windows
    close all;
    
    % Storage
    annotation_matrix = [];
    
    % IInit Structure with Default Settings
    viz = struct;
    viz.SamplingFrequency = 10000;
    viz.PeakWindow = [0.02 0.03]; % ms after Trigger Pulse
    viz.TriggerTime = 0; % Time in ms (from first data sample) at which TMS pulse is aplied
    viz.ChannelNames = {}; % Names of Channels == size(data,3)
    viz.unit = 'V'; % Display Unit
    viz.demean = '';
    
    % Assign Data File
    viz.Data = data;
    viz.NumChannels = 0;
    viz.Xlimit = nan; % Time in Seconds
    viz.Ylimit = nan; % Voltage in Milivolt
    viz.Routine_Active = [];
   
    % Overwrite default values with config input
    fnames = fieldnames(config);
    for f_idx = 1:size(fnames,1)
        viz.(fnames{f_idx}) = config.(fnames{f_idx}); 
    end
    
    % Demean Data
    if strcmp(viz.demean, 'yes')
        fprintf('Demeaning Data \n')
        viz.Data = hf_demean_data(viz.Data);
    end
    
    % Replace mu in strings for labels
    if contains(viz.unit, '?')
        viz.unit = strrep(viz.unit,'?','\mu');
    end
        
    % Find number of Channels
    if size(viz.Data, 3) == 1
        viz.NumChannels = 1;
    else
        viz.NumChannels = size(viz.Data, 3);
    end
    
    % Create Channel Labels
    if isempty(viz.ChannelNames)
        for ch = 1:viz.NumChannels
            viz.ChannelNames{end+1} = sprintf('Chan %d', ch);
        end
    end
    
    % Define CMAP
    viz.Cmap = colormap(hsv(3));
    viz.Cmap = viz.Cmap([2, 1, 3], :); % Green, Red, Blue
        
    % Stores Trial Information: [Categorization, Interval_On, Interval_Off, Peak to Peak Voltage]
    storage_dims = [size(viz.Data,2) 4 viz.NumChannels];
    viz.Annotation_Matrix = nan(storage_dims);
    
    % Set X-Axis Values
    viz.Current.Time = linspace(0,size(viz.Data,1)/viz.SamplingFrequency,size(viz.Data,1)) - viz.TriggerTime;
    
    % Set XLimit
    if any(isnan(viz.Xlimit))
        viz.Xlimit = [viz.Current.Time(1) viz.Current.Time(end)];
    end
    
    % Set YLimit (if not provided, take average of all to approximate)
    if any(isnan(viz.Ylimit))
        avg_data = mean(mean(viz.Data,3),2);
        viz.Ylimit = [-max(avg_data), max(avg_data)];
    end
    
    % Set Step Size
    viz.YlimitStepSize = mean(abs(viz.Ylimit))/10;

    % initiliase Main Figure Window
    FigH = figure(); hold on;
    FigH.Tag = 'MainWindow';
    set(FigH,'units','normalized','outerposition',[0 0 1 1])
    xlim(viz.Xlimit)
    ylim(viz.Ylimit)
    grid on; grid minor;
    xlabel('Time in [S]', 'FontSize', 14)
    ylabel(sprintf('Voltage [%s]', viz.unit), 'FontSize', 14)

    % Initialise Figure Handles and Tags
    MainDataH = plot(1,1, 'LineWidth', 1.5); MainDataH.Tag = 'MainData';
    MinMarkerH = plot(1,1, 'o-','MarkerFaceColor','red','MarkerFaceColor',[.49 1 .63] , 'MarkerSize', 9,  'MarkerEdgeColor', 'black'); MinMarkerH.Tag = 'MinMarker';
    MaxMarkerH = plot(1,1, 'o-','MarkerFaceColor','red','MarkerFaceColor',[.49 1 .63] , 'MarkerSize', 9,  'MarkerEdgeColor', 'black'); MaxMarkerH.Tag = 'MaxMarker';
    LowBorderH = line([1 1], [1 1], 'LineWidth', 1, 'Color', [0.3 0.3 0.3]); LowBorderH.Tag = 'LowBorder';
    UpperBorderH = line([1 1], [1 1], 'LineWidth', 1, 'Color', [0.3 0.3 0.3]); UpperBorderH.Tag = 'UpperBorder';
    MainTitleH = title('empty', 'FontSize', 18); MainTitleH.Tag = 'MainTitle';
    TextBoxPeakH = text(viz.Xlimit(2)*0.85,viz.Ylimit(2)*0.85,'empty', 'FontSize', 18); TextBoxPeakH.Tag = 'TextBoxPeak';

    % Start Routine
    viz = start_routine(viz);
    
    % For Return
    Annotation_Matrix = viz.Annotation_Matrix;

 end

function obj = start_routine(obj)

    % Current Plot
    obj.Routine_Active = 1;
    obj.Current.Channel = 1;
    obj.Current.Trial = 1;

    % Start Routine
    while obj.Routine_Active == 1

        % Prepare Plot
        obj = prepare_for_plot(obj, obj.Current.Channel, obj.Current.Trial);

        % Plot Both Curves;
        obj = update_figure(obj);

        % wait for keypress
        keypress = waitingForButtonPress();

        % evaluate keypress
        obj = evaluateKeyPress(obj,keypress);

    end
    % Return Current Matrix
    annotation_matrix = obj.Annotation_Matrix;
end

% Select Session
function obj = prepare_for_plot(obj, channel, trial)

    % Set Trial to Good by Default;
    if isnan(obj.Annotation_Matrix(trial, 1, channel))
        obj.Annotation_Matrix(trial, 1, channel) = 1;
    end

    % Set
    obj.Current.Channel = channel;
    obj.Current.Trial = trial;

    % Select Data;
    obj.Current.Data = obj.Data(:, trial, channel);
end
        
function obj = update_figure(obj)

    % Select Current EMG_Selection
    Select = obj.Annotation_Matrix(:, :, obj.Current.Channel);

    % Plot Main Data;
    objHan = findobj('Tag', 'MainData');
    objHan.Color = obj.Cmap(Select(obj.Current.Trial,1),:);

    warning('off');
    objHan.XData = obj.Current.Time;
    objHan.YData = obj.Current.Data;
    warning('on');

    % Update Borders of Target Interval
    updateIntervalBorder(obj)

    % Update Min Max Marker
    [minX,minY,maxX,maxY] = updateMinMaxMarker(obj);

    % Update Peak to Peak Voltage
    p2p_voltage = maxY - minY;
    objHan = findobj('Tag', 'TextBoxPeak');
    objHan.String = sprintf('%s %s', num2str(round(p2p_voltage,3)), obj.unit);
    %objHan.Position = [objHan.Position(1) 0.85*obj.Ylimit(1)];
    objHan.Position = [obj.Xlimit(2)*0.85 obj.Ylimit(2)*0.85];

    % Print Peak to Peak Voltage
    obj.Annotation_Matrix(obj.Current.Trial, 4, obj.Current.Channel) = p2p_voltage;
    fprintf('[Channel: %s, Trial: %d] - Peak-to-peak Voltage: %d %s\n',  obj.ChannelNames{obj.Current.Channel}, obj.Current.Trial, obj.Annotation_Matrix(obj.Current.Trial, 4, obj.Current.Channel), obj.unit)

    % Update Peak to peak Text in Plot
    objHan = findobj('Tag', 'TextBoxPeak');
    objHan.String = sprintf('%s %s', num2str(round(p2p_voltage,4)), obj.unit);

    % Set PeakWindow as Default Value for Currrent Interval;
    Interval = obj.Annotation_Matrix(obj.Current.Trial, [2, 3], obj.Current.Channel);

    % Replace NAN
    if isnan(Interval)
        obj.Annotation_Matrix(obj.Current.Trial, [2, 3], obj.Current.Channel) = obj.PeakWindow;
    end

    % Update Title;
    objHan = findobj('Type', 'Axes');
    objHan.Title.String = sprintf('Channel: %s || Trial %d/%d', obj.ChannelNames{obj.Current.Channel}, obj.Current.Trial, size(obj.Data,2));

end

%%%% HELPER FUNCTIONS;

function [minX,minY,maxX,maxY] = findMinMaxVoltage(obj)

% Select Interval previously defined;
Interval = obj.Annotation_Matrix(obj.Current.Trial, [2, 3], obj.Current.Channel);

if isnan(Interval)
    Interval = obj.PeakWindow;
end

% Find Index of Low & Upper Border
[~, peakIndexLB] = min(abs(obj.Current.Time-Interval(1)));
[~, peakIndexUB] = min(abs(obj.Current.Time-Interval(2)));

% Select Data
mep_data = obj.Current.Data;


% Select IV Part of Data
iv = mep_data(peakIndexLB:peakIndexUB);

% Find Indices
[minVal, minIdx] = min(iv);
[maxVal, maxIdx] = max(iv);

% Compute p2p Value;
p2p = peak2peak(mep_data(peakIndexLB:peakIndexUB));

% Return Values;
minX = obj.Current.Time(peakIndexLB+minIdx - 1);
minY = mep_data(peakIndexLB+minIdx - 1);

maxX = obj.Current.Time(peakIndexLB+maxIdx - 1);
maxY = mep_data(peakIndexLB+maxIdx - 1);

end


function [minX,minY,maxX,maxY] = updateMinMaxMarker(obj)

% update Peak 2 Peak Markers
[minX,minY,maxX,maxY] = findMinMaxVoltage(obj);

objHan = findobj('Tag', 'MinMarker');
objHan.XData = minX;
objHan.YData = minY;

objHan = findobj('Tag', 'MaxMarker');
objHan.XData = maxX;
objHan.YData = maxY;

end

function updateIntervalBorder(obj)

% Select Interval
Interval = obj.Annotation_Matrix(obj.Current.Trial, [2, 3], obj.Current.Channel);

% Replace Nan with Default Interval;
if isnan(Interval)
    Interval = obj.PeakWindow;
end

objHan = findobj('Tag', 'LowBorder');
objHan.XData = [round(Interval(1),3) round(Interval(1),3)];
objHan.YData = [obj.Ylimit(1) obj.Ylimit(2)];
objHan = findobj('Tag', 'UpperBorder');
objHan.XData = [round(Interval(2),3) round(Interval(2),3)];
objHan.YData = [obj.Ylimit(1) obj.Ylimit(2)];

end

function keypress = waitingForButtonPress()

% Waits for certain Button Presses

buttons = {'w', 's', 'a','d', 'q', 'e', 'f', 'p', 't', 'm', 'n', 'h', '1', '2', '3', '0', '9', 'escape'};

while 1
    w = waitforbuttonpress;
    if w == 1
        keypress = get(gcf,'CurrentKey');
        
        if any(strcmp(buttons,keypress))
            break
        end
    end
end
end


function obj = evaluateKeyPress(obj,keypress)

% Evaluates Button Presses

% Assignment for easier handling;
channel = obj.Current.Channel;
trial = obj.Current.Trial;
emg = obj.Annotation_Matrix(:, :, channel);

% Do Action on KeyPress
if keypress == 'a'
    % Previous Trial;
    
    if trial ~= 1
        trial = trial - 1;
    end
    
elseif keypress == 'd'
    %D = show next trial
    
    % Default: Trial is good.
    if isnan(emg(trial,1))
        emg(trial,1) = 1;
    end
    
    if trial +1 > size(emg,1)
        disp('LAST Trial')
    else
        trial = trial + 1;
    end
    
elseif keypress == 'w'
    %W = mark trial as good
    emg(trial,1) = 1;
    
elseif keypress == 's'
    
    %S = mark trial as bad
    emg(trial,1) = 2;
    
elseif keypress == 'q'
    %Q = jump 10 trials backwards
    if trial > 11
        trial = trial - 10;
    else
        trial = 1;
    end
    
elseif keypress == 'e'
    %E = jump 10 trials foward
    if trial < size(emg,1)-10
        trial = trial + 10;
    else
        trial = 1;
    end
    
elseif keypress == '1'
    
    %1 = Mark MEP as Good
    emg(trial,1) = 1;
    
elseif keypress == '2'
    
    %2 = Mark MEP as Bad
    emg(trial,1) = 2;
    
elseif keypress == '3'
    
    %3 = Mark Trial as Other;
    emg(trial,1) = 3;
    
elseif keypress == 'f'
    %F = mark interval for computing peak 2 peak amplitude
    try
        rect = getrect;
        xLocations = [rect(1) rect(1)+rect(3)];
        
        % Get Index of Window
        [~, xloc_on] = min(abs(obj.Current.Time-xLocations(1)));
        [~, xloc_off] = min(abs(obj.Current.Time-xLocations(2)));
        
        % Get Time
        time_on = obj.Current.Time(xloc_on);
        time_off = obj.Current.Time(xloc_off);
        
        emg(trial,2:3) = [time_on time_off];
        sprintf('Interval set to %d - %d ms.', time_on, time_off)
        
    catch
        %do nothing
    end
    
elseif keypress == 't'
    %T = Reset Target Interval to Standard;
    emg(trial,2:3) = obj.PeakWindow;
    sprintf('[Chan: %d, Trial: %d] Interval set to %d - %d ms.', obj.Current.Channel, obj.Current.Trial, obj.PeakWindow(1), obj.PeakWindow(2))
    
elseif keypress == '9'
    %9 = reduce size of y axis
    obj.Ylimit = [obj.Ylimit(1)*0.9 obj.Ylimit(2)*0.9];

    ylim(obj.Ylimit)
    obj.Ylimit
    
elseif keypress == '0'
    %0 = increase size of y axis
    obj.Ylimit = [obj.Ylimit(1)*1.1 obj.Ylimit(2)*1.1];
    ylim(obj.Ylimit)
    obj.Ylimit
    
    
elseif keypress == 'escape'
    %ESC  = QUIT
    
    % Stop Routine
    obj.Routine_Active = 0;
    fprintf('Quit without saving.  \n')
    close all;
    
elseif keypress == 'm'
    
    % Switch Channel
    if channel == obj.NumChannels
        channel = 1;
    else
        channel = channel + 1;
    end

    % Reassign Channel, VERY IMPORTANT
    emg = obj.Annotation_Matrix(:, :, channel);
    
elseif keypress == 'n'

    % Switch Channel
    if channel == 1
        % Go To Max
        channel = obj.NumChannels;
    else
        channel = channel - 1;
    end

    % Reassign Channel, VERY IMPORTANT
    emg = obj.Annotation_Matrix(:, :, channel);
    
elseif keypress == 'h'
    % Call Button Help Function
    display_button_functions()
    
elseif keypress == 'p'
    
    % P = SAVE & QUIT
    % Call GUI
    [file,path,indx] = uiputfile('*.mat','Save File & Quit');

    % Select Data
    Annotation_Matrix = obj.Annotation_Matrix;

    % Save
    save([path file], 'Annotation_Matrix')
    
    % QUIT
    obj.Routine_Active = 0;
    fprintf('Saved %s & quit. \n', [path file])
    close all;
end

% Reassign;
obj.Annotation_Matrix(:, :, channel) = emg;
obj.Current.Trial = trial;
obj.Current.Channel = channel;
end


function output_data = hf_demean_data(data)

% storage
output_data = nan(size(data));

% For all Channels;
for nChan = 1:size(data,3)
    
    % Select Channel
    chan = data(:,:,nChan);
    
    % Detrend & Demean
    for nTrial = 1:size(chan,2)

        % Demean
        chan(:,nTrial) = chan(:,nTrial) - mean(chan(:,nTrial));
    end
    
    % Reassign;
    output_data(:,:,nChan) = chan;
end
end


function display_button_functions()

fprintf('h: show keymappings \n')
fprintf('w: annotate trial as good \n')
fprintf('s: annotate trial as bad \n')
fprintf('a: show previous trial \n')
fprintf('d: show next trial \n')
fprintf('m: show next channel \n')
fprintf('n: show previous channel \n')
fprintf('q: show previous trial (skip 10) \n')
fprintf('e: show next trial (skip 10) \n')
fprintf('f: select a window to compute peak-to-peak amplitude (drag over area) \n')
fprintf('t: reset interval to default value \n')
fprintf('1: annotate trial as good \n')
fprintf('2: annotate trial as bad \n')
fprintf('3: annotate trial as other \n')
fprintf('9: zoom in (y-axis) \n')
fprintf('0: zoom out (y-axis) \n')
fprintf('p: save selection & quit \n')
fprintf('escape: quit without saving \n')

end


