function x = iAACoder3(AACSeq3, fNameOut)
% function x = iAACoder3(AACSeq3, fNameOut)
% 
% This function decodes the track originally coded by AACoder1
% arguements:
%   AACSeq3: struct of size Kx1 where K is the number of coded frames
%         AACSeq3(i).frameType     : The type of the frame:
%                                     "OLS","LSS","ESH","LPS"
%         AACSeq3(i).winType       : The method selected for coding the frame:
%                                     "KBD" or "SIN"
%         AACSeq3(i).chl.TNScoeff : The quantized TNS coefficents of the
%                                     left channel
%         AACSeq3(i).chr.TNScoeff : The quantized TNS coefficents of the
%                                     right channel
%         AACSeq3(i).chl.T: The thresholds of the left channel frames
%         AACSeq3(i).chr.T: The thresholds of the right channel frames
%         AACSeq3(i).chl.G: The quantized global gains for the left
%                               channel
%         AACSeq3(i).chl.G: The quantized global gains for the right
%                               channel
%         AACSeq3(i).chl.sfc: The Huffman coded sfc  sequence for the left
%                               channel
%         AACSeq3(i).chl.sfc: The Huffman coded sfc  sequence for the right
%                               channel
%         AACSeq3(i).chl.stream: The Huffmann coded sequence of the
%                                   quantized MDCT coefficients for the
%                                   left channel
%         AACSeq3(i).chr.stream: The Huffmann coded sequence of the
%                                   quantized MDCT coefficients for the
%                                   right channel
%         AACSeq3(i).chl.codebook: The huffman codebook used for the left
%                                   channel
%         AACSeq3(i).chr.codebook: The huffman codebook used for the right
%                                   channel
%   fNameOut: char array containing the name of the file name for the
%               decoded signal
%
% Assumptions: track is 2-channel signal sampled at 48 kHz

    %find number of frames to decode
    numFrames = length(AACSeq3);
    
    %arrays for right and left frames, since its about 2-channel signals
    frameTR = zeros(2048,numFrames);
    frameTL = frameTR;
    
    % load Look-Up Table
    huffLUT = loadLUT();
    
    
    for i=1:numFrames
        
        %% Huffmann Decoding
        decCoeffs_R = decodeHuff( AACSeq3(i).chr.stream, ...
                                AACSeq3(i).chr.codebook,...
                                huffLUT);
        decCoeffs_L = decodeHuff( AACSeq3(i).chl.stream, ...
                                AACSeq3(i).chl.codebook,...
                                huffLUT);

        %% De-quantizer
         frameF_R = iAACquantizer(decCoeffs_R, ...
                                  AACSeq3(i).chr.sfc,...
                                  AACSeq3(i).chr.G,...
                                  AACSeq3(i).frameType);
         frameF_L = iAACquantizer(decCoeffs_L, ...
                                  AACSeq3(i).chl.sfc,...
                                  AACSeq3(i).chl.G,...
                                  AACSeq3(i).frameType);
        
        %% iTNS

        %apply iTNS to left and right channel of the frame
        frameF_R = iTNS(frameF_R, ...
                        AACSeq3(i).frameType,...
                        AACSeq3(i).chr.TNScoeff);
        frameF_L = iTNS(frameF_L, ...
                        AACSeq3(i).frameType,...
                        AACSeq3(i).chl.TNScoeff);
        
        %get frameF
        frameF = [frameF_R , frameF_L];
        %% iFilterbank
        
        frameT = ifilterbank( frameF, AACSeq3(i).frameType,...
            AACSeq3(i).winType);
        %split the output into right and left sections
        frameTR(:,i) = frameT(:,1);
        frameTL(:,i) = frameT(:,2);
    end
    
    frameTR(:,1) = [zeros(1024,1) ; frameTR(1025:end,1)];
    frameTL(:,1) = [zeros(1024,1) ; frameTL(1025:end,1)];
    %so far we have the time domain frames
    
    %unbuffer them
    frameLength = 2048;
    R = zeros(frameLength + numFrames*(frameLength/2),1);
    L = zeros(frameLength + numFrames*(frameLength/2),1);
    
    for i=1:numFrames
        pos = (i-1)*(frameLength/2);
        R((pos+1):(pos+frameLength)) = R((pos+1):(pos+frameLength)) + frameTR(:,i);
        L((pos+1):(pos+frameLength)) = L((pos+1):(pos+frameLength)) + frameTL(:,i);
    end
    
    %trim extra points
    R = R(1025:(end - 1024 - 670));
    L = L(1025:(end - 1024 - 670));
    
    x = [R , L];
    
    audiowrite(fNameOut, x, 48000); 

end