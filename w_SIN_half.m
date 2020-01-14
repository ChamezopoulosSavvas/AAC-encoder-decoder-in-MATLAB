% function that creates ether the left or the right half of SIN Window
% N : number of samples
% half : Side of window to be calculated , gets values "left" or "right"

function w_sin = w_SIN_half(N, half)

    w_sin = zeros(N/2,1);
    
    if half == "left"
        %calc left win
        for n=0:(N/2-1)
            w_sin(n+1,1) = sin((pi/N) * (n + 1/2));
        end
    elseif half == "right"
        %calc right win
        for n=N/2:(N-1)
            w_sin(n-N/2+1,1) = sin((pi/N) * (n + 1/2)); %this in order the return vector to be 0ews128
        end
    end
end