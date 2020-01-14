function AACSeq3 = AACoder3(fNameIn, fnameAACoded)

% function AACSeq3 = AACoder3(fNameIn, fnameAACoded)
%
% This function implements level 3 coding. Takes a 2-channel uncompressed
% .wav file sampled at 48kHz and returns AACSeq3 struct. It also saves the
% AACSeq3 into an external file.
%
% arguments:
%   fNameIn = (string) the filename of the track
%   fnameAACoded = (string) the fileneme the struct will be saved
% return value:
%    AACSeq3: struct of size Kx1 where K is the number of coded frames
%         AACSeq3(i).frameType     : The type of the frame:
%                                     "OLS","LSS","ESH","LPS"
%         AACSeq3(i).winType       : The method selected for coding the frame:
%                                     "KBD" or "SIN"
%         AACSeq2(i).chl.frameF    : The MDCT of the left-channel frame
%                                     AFTER TNS
%         AACSeq2(i).chr.frameF    : The MDCT of the right-channel frame
%                                     AFTER TNS
%         AACSeq3(i).chl.TNScoeff : The quantized TNS coefficents of the
%                                     left channel
%         AACSeq3(i).chr.TNScoeff : The quantized TNS coefficents of the
%                                     right channel
%         AACSeq3(i).chl.T: The thresholds of the left channel frames
%         AACSeq3(i).chr.T: The thresholds of the right channel frames
%         AACSeq3(i).chl.G: The quantized global gains for the left
%                               channel
%         AACSeq3(i).chr.G: The quantized global gains for the right
%                               channel
%         AACSeq3(i).chl.sfc: The Huffman coded sfc  sequence for the left
%                               channel
%         AACSeq3(i).chr.sfc: The Huffman coded sfc  sequence for the right
%                               channel
%         AACSeq3(i).chl.stream: The Huffmann coded sequence of the
%                                   quantized MDCT coefficients for the
%                                   left channel
%         AACSeq3(i).chr.stream: The Huffmann coded sequence of the
%                                   quantized MDCT coefficients for the
%                                   right channel
%         AACSeq3(i).chl.codebook: ?? ?uffman codebook used for the left
%                                   channel
%         AACSeq3(i).chr.codebook: ?? ?uffman codebook used for the right
%                                   channel

    info = audioinfo(fNameIn);
    if(strcmp(info.CompressionMethod, 'Uncompressed') == 0 ||...
            info.SampleRate ~= 48000)
        disp("ERROR: Wrong compression method and/or sample rate, exiting...");
        return;
    end

    % load boards B.2.1.9.a/b
    vars = {'B219a','B219b'};
    load('Tableb219.mat', vars{:});
    
    %incr first 3 columns by 1 so they match matlab array indexes
    B219a(:,1:3) = B219a(:,1:3) + 1; %#ok<NODEF>
    B219b(:,1:3) = B219b(:,1:3) + 1; %#ok<NODEF>

    [y, fs] = audioread(fNameIn); %#ok<ASGLU>

    frameLength = 2048;

    %separate the 2 channels
    R = y(:,1);
    L = y(:,2);

    %get the frames
    framesR = buffer(R,frameLength,frameLength/2);
    framesL = buffer(L,frameLength,frameLength/2);
    
    numFrames =  size(framesR,2);
    SMR_R = ones(42*8,numFrames);
    SMR_L = ones(42*8,numFrames);

    %Construct struct
    ch = struct('frameF', zeros(frameLength,1),...
        'TNScoeff', zeros(4*8,1),...
        'T', zeros(length(B219b)*8,1), 'G', zeros(1,8), ...
        'sfc', zeros(length(B219b)*8,1), ...
        'stream', ' ', 'codebook', 0.0);

    for i=1:numFrames
        AACSeq3(i,1) = struct('frameType', "", 'winType', "",...
            'chl', ch, 'chr', ch); %#ok<AGROW>
    end

    % load Look-Up Table
    huffLUT = loadLUT();


    %% Sequence Segmentation Control
    for i=1:numFrames

        %select previous frame type
        if(i>1)
            prevFrameType = AACSeq3(i-1).frameType;
        else
            prevFrameType = "OLS";
        end

        %select next frame type
        if(i<numFrames)
            nextFrameT = [framesR(:,i+1), framesL(:,i+1)];
        else
            nextFrameT = zeros(frameLength,2);
        end

        frameT = [framesR(:,i), framesL(:,i)];
        %characterize frame
        AACSeq3(i).frameType = SSC(frameT, nextFrameT, prevFrameType);
        %choose window type
        %AACSeq3(i).winType = "SIN";
        AACSeq3(i).winType = "KBD";
        
        %choose prev1 and prev2 frames
        if(i == 1)
            frameTprev1R = zeros(frameLength,1);
            frameTprev2R = frameTprev1R;
            frameTprev1L = zeros(frameLength,1);
            frameTprev2L = frameTprev1L;
        elseif (i == 2)
            frameTprev1R = framesR(:,i-1);
            frameTprev2R = zeros(frameLength,1);
            frameTprev1L = framesL(:,i-1);
            frameTprev2L = zeros(frameLength,1);
        else
            frameTprev1R = framesR(:,i-1);
            frameTprev2R = framesR(:,i-2);
            frameTprev1L = framesL(:,i-1);
            frameTprev2L = framesL(:,i-2);
        end

        %% psycho
        SMR_R(:,i) = psycho(framesR(:,i), AACSeq3(i).frameType, frameTprev1R, frameTprev2R);
        SMR_L(:,i) = psycho(framesL(:,i), AACSeq3(i).frameType, frameTprev1L, frameTprev2L);
        
        %% Filterbank
        frameF = filterbank( frameT, AACSeq3(i).frameType, AACSeq3(i).winType);
        
        %store it appropiately
        frameF_R = frameF(:,1);
        frameF_L = frameF(:,2);

        %% TNS
        [frameF_R, AACSeq3(i).chr.TNScoeff] = ...
            TNS(frameF_R, AACSeq3(i).frameType);

        [frameF_L, AACSeq3(i).chl.TNScoeff] = ...
            TNS(frameF_L, AACSeq3(i).frameType);
        
        
        
         AACSeq3(i).chr.frameF = frameF_R;
         AACSeq3(i).chl.frameF = frameF_L;
        %% Threshold

        if (AACSeq3(i).frameType == "ESH")
            AACSeq3(i).chr.T = threshold(frameF_R , ...
                AACSeq3(i).frameType, ...
                SMR_R(:,i), B219b);
            AACSeq3(i).chl.T = threshold(frameF_L , ...
                AACSeq3(i).frameType, ...
                SMR_L(:,i), B219b);
        else
            AACSeq3(i).chr.T = threshold(frameF_R , ...
                AACSeq3(i).frameType, ...
                SMR_R(:,i), B219a);
            AACSeq3(i).chl.T = threshold(frameF_L , ...
                AACSeq3(i).frameType, ...
                SMR_L(:,i), B219a);
        end
        
        
        %% quantizer
        [coeffSec_R, AACSeq3(i).chr.sfc , AACSeq3(i).chr.G] = ...
                                AACquantizer(frameF_R,...
                                             AACSeq3(i).frameType,...
                                             SMR_R(:,i));
        [coeffSec_L, AACSeq3(i).chl.sfc, AACSeq3(i).chl.G] = ...
                                AACquantizer(frameF_L,...
                                             AACSeq3(i).frameType,...
                                             SMR_L(:,i));
                                         
        AACSeq3(i).chr.frameF = coeffSec_R;
        AACSeq3(i).chl.frameF = coeffSec_L;
        
        
        %% Huffman
        [AACSeq3(i).chr.stream, AACSeq3(i).chr.codebook] = ...
            encodeHuff(coeffSec_R, huffLUT); %, forcedCodebook);

        [AACSeq3(i).chl.stream, AACSeq3(i).chl.codebook] = ...
            encodeHuff(coeffSec_L, huffLUT); %, forcedCodebook);
    end
    
    save(fnameAACoded, 'AACSeq3');
end