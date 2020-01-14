function [SNR, bitrate, compression] = demoAAC3(fNameIn, fNameOut, fNameAACoded)

%function [SNR, bitrate, compression] = demoAAC3(fNameIn, fNameOut, frameAACoded)
%
% function that demonstrates the third level of decoding
% 
% arguments:
%     fNameIn : char array that contains the name of the original file
%     fNameOut: char array that contains the name of the coded-decoded file
%     fNameAACoded: .mat file that AACSeq3 stuct is stored
%
%return values:
%     SNR: the SNR value for the two tracks
%     bitrate: bits per sec of coded-decoded track
%     compression: bitrate before coding over bitrate after coding
    
    
    %% SNR
    [y, fsy] = audioread(fNameIn); %#ok<ASGLU>
    [x, fsx] = audioread(fNameOut); %#ok<ASGLU>

    %diff = length(x) - length(y);
    %channel 1
    sigma_x=var(x(:,1));
    sigma_e=var(x(:,1) - y(:,1));
    SNR(1) = 10 * log(sigma_x^2/sigma_e^2);
    %channel 2
    sigma_x=var(x(:,2));
    sigma_e=var(x(:,2) - y(:,2));
    SNR(2) = 10 * log(sigma_x^2/sigma_e^2);
    
    %% BitRate
    load(fNameAACoded , 'AACSeq3');
    
    lenR = zeros(length(AACSeq3),1);
    lenL = zeros(length(AACSeq3),1);
    for i=1:length(AACSeq3)
        lenR(i) = length(AACSeq3(i).chr.stream);
        lenL(i) = length(AACSeq3(i).chl.stream);
    end
    
    bitrate = mean([mean(lenR) ; mean(lenL)]);
    
    %% Compression
    infoIn = audioinfo(fNameIn);
    infoOut = audioinfo(fNameOut);
    
    bitsIn = infoIn.BitsPerSample*infoIn.TotalSamples*infoIn.NumChannels;
    KbytesIn = bitsIn/(8*1024);
    
    bitsOut = infoOut.BitsPerSample*infoOut.TotalSamples*infoOut.NumChannels;
    KbytesOut = bitsOut/(8*1024);
    
    KbytesOut = KbytesOut/8;
    
    fprintf('Unompressed audio: %f MB (%d bits)\n', KbytesIn/1024, bitsIn); 
    fprintf('Compressed struct: %f KB (%d bits)\n', KbytesOut, bitsOut); 
   
    compression = KbytesOut/KbytesIn;
    

end