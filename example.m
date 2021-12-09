% Close all windows & clear variables
close all; clear;

% Load Example Data File. Contains 10 Motor-evoked potentials (MEP) recorded in 3 different channels.
% Each MEP contains of 750 samples, acquired at a sample rate of 10kHz
load('data.mat')

% Scaling of data from V to mV
data = data * 1000;

% Config to set parameters
config.PeakWindow = [0.02 0.035]; % Target Interval for calculation of peak-to-peak amplitude
config.SamplingFrequency = 10000; % Sampling Frequency in Hz
config.ChannelNames = {'EMG1', 'EMG2', 'EMG3'}; % Some simple Channel Names
config.TriggerTime = 0.005; % Time of TMS pulse. Pulse was delievered 5ms after the first data point.
config.demean = 'yes'; % Demean each trial
config.unit = 'mV'; % Unit for axis-labels; 'V', 'mV', '?V'

% Call Function
Annotation_Matrix = MEP_Visualizer(data, config);

