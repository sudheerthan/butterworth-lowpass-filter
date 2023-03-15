% Define time and frequency parameters
fs = 100000; % Sampling rate
t = 0:1/fs:1-1/fs; % Time vector
f1 = 6000; % Frequency of 6kHz wave
f2 = 6000; % Frequency of sine/cosine waves
 t2 = 0:0.1:.01;
% Create signals
x = 3*sin((2*pi*f1*t)); % 6kHz wave

s1 = 10*sin(2*pi*f1*t); % Sine wave
c1 = 10*cos(2*pi*f1*t); % Cosine wave

% Multiplied sine signals
y1 = x .* s1; % Sine multiplication
subplot(5,1,1);
plot(t, y1);
title('multi sin Signal 1');


% Multiplied cosine signals
y2 = x .* c1; % Cosine multiplication
subplot(5,1,2);
plot(t, y2);
title('multi cosine Signal 1');


% Define filter parameters
f_cutoff = 100;                      % Cutoff frequency
[b,a] = butter(2, f_cutoff/(fs/2));  % 2nd order Butterworth LPF coefficients

% Filter the signal in phase
yfs = filter(b, a, y1);
subplot(5,1,3);
plot(t, yfs);
title('filter sine Signal 1');


% Filter the signal quadarature 
yfc = filter(b, a, y2);
subplot(5,1,4);
plot(t, yfc);
title('filter cosine Signal 1');

%calculate the input voltage
c = sqrt(yfs.^2 + yfc.^2);
%yfs = filter(b, a, c);
subplot(5,1,5);
plot(t, c);
title('amplitude 1');
