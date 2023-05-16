# What is EP_Visualizer?
EP_Visualizer is a simple MATLAB tool that allows user to visualise and inspect biosignals, such as motor-evoked potentials (MEP) or other evoked potentials (EP).
It can be used to extract the peak to peak amplitude of the signals in a pre-defined or manually selected time window and allows for annotation of low quality signals for later rejection.

# Installation
The visualizer is a simple script that should work by simply copying it into your current working folder or any other folder that is part of the MATLAB search path.

# How to use EP_Visualizer
The EP_Visualizer can be called with a simple one-line command:

output_matrix = MEP_Visualizer(data)\
or \
output_matrix = MEP_Visualizer(data, config) if an additional config structure is provided.

The input 'data' should either have two dimensions (Time Samples x Trials) or three dimensions (Time Samples x Trials x  Channels). The input sfreq corresponds to the sampling frequency of the recorded dataset.

The config structure allows the user to modify some basic settings, such as:\
config.TriggerTime       = 0; - Occurence of Pulse/Trigger; Time in [S] relative to the first sample; If provided, this is used as reference point for the PeakWindow\
config.PeakWindow        = [0.02 0.04] - Time Window for calculation of peak to peak amplitude; >> Time in [S] relative to Trigger Time<<\
config.ChannelNames      = {'Chan 1', 'Chan 2'} - Cell array with Channel Names. Size needs to be == size(data,3)\
config.demean            = 'no' - demean each single trial, can be 'yes' or 'no'\
config.unit              = 'V' - Only for labeling of axis/text\
config.Xlimit            = [min_value max_value] - Set specific X-Axis Limits;\
config.Ylimit            = [min_value max_value] - Set specific Y-Axis Limits;

The output of the function corresponds to a matrix of Dimensions (4, Input_Shape). The four values of the novel dimension in the output matrix describe some properties of the trial:\
1 - Annotation of the trial (1 = good, 2 = bad, 3 = other)\
2 - Onset in [s] of peak-to-peak amplitude calculation window\
3 - Offset in [s] of peak-to-peak amplitude calculation window\
4 - Peak to Peak Amplitude [µV]

Once the function is called, the user can navigate through the dataset and annotate the trials using several keys:\

% - Press 'h' to show all key bindings.\
% - Use 'a' and 'd' to navigate through the trials\
% - Use 'n' and 'm' to cycle through the channels\
% - Use 'w' and 's' to annotate a trial as good and bad, respectively\
% - Keys 1, 2, and 3 can also be used for annotation; 1 = good, 2 = bad, 3 = other\
% - Use 'f' to draw a time window with your mouse to define an individual Peak Window Time\
% - Press 't' to reset the Peak Window Time to default value\
% - Press 'p' to save annotation_file

![alt text](https://github.com/XamHom/MEP_Visualizer/blob/main/image1.png?raw=true)

