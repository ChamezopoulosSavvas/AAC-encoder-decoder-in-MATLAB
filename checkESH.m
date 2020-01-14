% function esh = checkESH( frame, H, initGap, endGap, frameLength, subframeLength)
% 
% This function checks a given frame (2048x1) if it's of type "ESH"
% arguments:
%     frame:          the frame under investigation (2048x1)
%     H:              the high pass filter needed
%     initGap:        the number of values to skip before the start of the
%                     segmentation of the frame into subframes
%     endGap:         the number of values to stop before the end of the frame
%                     segmenting into subframes
%     frameLength:    the length of the frame
%     subframeLength: the length of the subframe
% return value:
%     esh: a flag set to 1 if the frame is of type "ESH", 0 else.
%

function esh = checkESH( frame, H, initGap, endGap, frameLength, subframeLength)

    esh = 0;
    
    %time steps
    t_frame = (0:(frameLength - 1))';
    %filter frame
    filtered = lsim(H, frame, t_frame);
    
    %segment frame into subframes
    subframes = buffer(filtered(...
        (initGap + (subframeLength/2)):(end - endGap - 1)),...
        subframeLength);
    %skip last column as it is half filled with zeroes
    subframes = subframes(:,1:(end-1));
    
    %calculate (s_l)^2
    ssq_l = (sum(subframes.^2))';
    
    %calculate attack values
    for j=2:length(ssq_l)
        dssq_l = ssq_l(j)/((1/(j))*sum(ssq_l(1:(j-1))));
        if((dssq_l > 10) && (ssq_l(j) > 1e-03))
            %if the criteria are fulfilled set the flag and stop
            esh = 1;
            break;
        end
    end
end