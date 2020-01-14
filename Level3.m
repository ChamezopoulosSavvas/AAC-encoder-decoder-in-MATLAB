%% init

close all;
clear variables;
clc;

%Encode
fprintf('Level 3 demonstration:\n\n');
tic;
AACSeq3 = AACoder3('LicorDeCalandraca.wav','AACSeq3.mat');
tcode = toc;
fprintf('Coding Elapsed Time: %f s\n',tcode);

%Decode
tic;
x = iAACoder3(AACSeq3, 'Level3Demo.wav');
tdecode = toc;
fprintf('Decoding Elapsed Time: %f s\n',tdecode);

%Evaluate quality
[SNR, bitrate, compression] = demoAAC3('LicorDeCalandraca.wav',...
                                       'Level3Demo.wav',...
                                       'AACSeq3.mat');
fprintf('Compresion ratio: %.2f %% \n', compression*100);
fprintf('Channel 1 SNR: %f dB\nChannel 2 SNR: %f dB\n',SNR(1),SNR(2));