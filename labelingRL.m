% function frameType = labelingRL(frameTypeR, frameTypeL)
% 
% function that determines the type of a stereo signal frame 
% based on the type of the right and left frame
% 
% arguments:
%     frameTypeR: (string) the type of the right frame.
%     should be one of the following:
%         "OLS": for ONLY_LONG_SEQUENC
%         "LSS": for LONG_START_SEQUENCE
%         "ESH": for EIGHT_SHORT_SEQUENCE
%         "LPS": for LONG_STOP_SEQUENCE
%     frameTypeL: same as FrameTypeR, only for the left frame.  
%           
%     
% return value:
%    frameType: (string) same as FrameTypeR, but for the total frame.

function frameType = labelingRL(frameTypeR, frameTypeL)
    %types = ["OLS"; "LSS"; "ESH"; "LPS"];
    
    %type if fail
    frameType = "N/A";
    
    if( frameTypeR == "OLS" ) 
        frameType = frameTypeL;
    elseif ( frameTypeR == "LSS" ) 
        if(frameTypeL == "ESH" || frameTypeL == "LPS")
            frameType = "ESH";
        else
            frameType = "LSS";
        end
    elseif ( frameTypeR == "ESH" )
            frameType = "ESH";
    elseif ( frameTypeR == "LPS" )
        if ( frameTypeL == "OLS" || frameTypeL == "LPS" )
            frameType = "LPS";
        else
            frameType = "ESH";
        end
    end
end