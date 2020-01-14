% function that creates ether the left or the right half of WDB Window
% N : number of samples
% a : factor for kaiser function of matlab, gets values 6 or 4
% half : Side of window to be calculated , gets values left or right

function w_kbd = w_KBD_half(N, a, half)

    w_kbd = ones(N/2,1);
    w_kaiser = kaiser(N/2+1,a*pi) ;
    
    %denumerator is common for left and right windows
    sum2 = sum(w_kaiser(1:((N/2)+1)));
    if half == "left"
        %calc left Win
        %sum1=0; 
        %sum2 = sum(w_kaiser(1:((N/2)+1)));
        for n=0:((N/2)-1)
           sum1= sum(w_kaiser(1:(n+1)));
%            for i=0:n
%                sum1 = sum1 + w_kaiser(i+1);
%            end
%            for i=0:N/2
%                sum2 = sum2 + w_kaiser(i+1);
%            end
           w_kbd(n+1,1) = sqrt(sum1/sum2);
        end
        
    elseif half == "right"
        %calc right win
        %sum1=0; 
        %sum2 = sum(w_kaiser(1:((N/2)+1)));
        
        for n=N/2:(N-1)
            
            sum1 = sum(w_kaiser(1:(N-n)));
%            for i=0:N-n
%                sum1 = sum1 + w_kaiser(i+1);
%            end
%            for i=0:N/2
%                sum2 = sum2 + w_kaiser(i+1);
%            end
            w_kbd(n-N/2+1,1) = sqrt(sum1/sum2);
        end
    end
    %length(w_kbd);
end
    
        
    
