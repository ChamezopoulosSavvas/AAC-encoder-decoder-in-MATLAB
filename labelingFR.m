% function frameType = labelingFR(prevFrameType, nextFrameType)
% 
% function that determines the type of a frame based on the type
% of the previous and the next frame
% 
% arguments:
%     prevFrameType: (string) the type of the previous frame.
%     one of the following:
%         "OLS": for ONLY_LONG_SEQUENC
%         "LSS": for LONG_START_SEQUENCE
%         "ESH": for EIGHT_SHORT_SEQUENCE
%         "LPS": for LONG_STOP_SEQUENCE
%     nextFrameType: same as prevFrameType, only for the following frame.  
%           
%     
% return value:
%    frameType: (string) same as prevFrameType, but for the current frame.

function frameType = labelingFR(prevFrameType, nextFrameType)
   
    %fail return value
    frameType = "N/A";
    
    %case 1
    if(prevFrameType == "OLS")
        if(nextFrameType == "ESH")
            frameType = "LSS";
        else
            frameType = "OLS";
        end
    
    
    %case 2
    elseif(prevFrameType == "ESH")
        if(nextFrameType == "ESH")
            frameType = "ESH";
        else
            frameType = "LPS";
        end
    %case 3
    elseif(prevFrameType == "LSS")
        frameType = "ESH";
    elseif(prevFrameType == "LPS")
        frameType = "OLS";
    end
end