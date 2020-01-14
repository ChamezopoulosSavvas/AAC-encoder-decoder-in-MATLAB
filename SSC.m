% frameType = SSC(frameT, nextFrameT, prevFrameType)
% 
% function that determines the type of a frame
% 
% arguments:
%     frameT:         frame(i) in time domain. array of size 2048x2
%     nextFrameT:     frame(i+1) in time domain. array of size 2048x2
%     prevFrameType:  type of frame (i-1). string. one of 
%                     "OLS", "LSS", "ESH", "LPS"
%  return value:
%     frameType: the type of the frameT. one of:
%                "OLS": for ONLY_LONG_SEQUENCE
%                "LSS": for LONG_START_SEQUENCE
%                "ESH": for EIGHT_SHORT_SEQUENCE
%                "LPS": for LONG_STOP_SEQUENCE

function frameType = SSC(frameT, nextFrameT, prevFrameType)
    
    %set the constants
    frameLength = length(frameT);
    subframeLength = 128;
    initGap = 448;
    endGap = 448;
    %high pass filter
    H = tf([0.7548 -0.7548],[1 -0.5095],-1);
    
    %split the frame
    R = frameT(:,1);  %#ok<NASGU>
    L = frameT(:,2); %#ok<NASGU>
    
    %split the next frame and check if its "ESH"
    nextR = nextFrameT(:,1);
    nextL = nextFrameT(:,2);
    eshNextR = checkESH(nextR, H, initGap, endGap, frameLength, subframeLength);
    eshNextL = checkESH(nextL, H, initGap, endGap, frameLength, subframeLength);
     
    nextFrameTypeR = "N/A";
    nextFrameTypeL = "N/A";
    if(eshNextR == 1)
        nextFrameTypeR = "ESH";
    end
    if(eshNextL == 1)
        nextFrameTypeL = "ESH";
    end
    %find type for the two channels separately
    frameTypeR = labelingFR(prevFrameType, nextFrameTypeR);
    frameTypeL = labelingFR(prevFrameType, nextFrameTypeL);
    %combine the two types to find the final type of the frame
    frameType = labelingRL(frameTypeR, frameTypeL);
end