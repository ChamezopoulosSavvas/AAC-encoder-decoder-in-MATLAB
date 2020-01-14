function s_w = hannWin(frameT, frameType)

%function s_w = hannWin(frameT, frameType)
%
% This function takes a frame and passes it through a symmetric Hanning
% Window.
% arguments:
%   frameT:    the frame in time-domain (2048x1)
%   frameType: the type of the frame
% return value:
%   s_w: the frame after it passes through the window

    s_w = zeros(size(frameT));
    
    if frameType == "ESH"
        
        s = buffer(frameT,256);
        N = length(s);    
        s_tmp = zeros(size(s));
        
        for c=1:size(s,2)
            for i=1:N
                s_tmp(i,c) = s(i,c)*(0.5 - 0.5*cos((pi*(i + 0.5))/(N/2)));
            end
            pos = (c-1)*256;
            s_w(pos+1:pos+256) = s_tmp(:,c);
        end
    else
        
        N = length(frameT);
        
        for i=1:N
            s_w(i) = frameT(i)*(0.5 - 0.5*cos((pi*(i + 0.5))/(N/2)));
        end
    end
end