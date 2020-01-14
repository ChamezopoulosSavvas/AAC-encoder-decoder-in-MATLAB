function T = threshold(frameF , frameType, SMR, B219)

% function T = threshhold(frameF , frameType, SMR, B219)
% 
% Calculates the thresholds for each band of a frame
% arguments:
%   frameF: the frame in frequency domain AFTER TNS (1024x1)
%   frameType: the type of the frame as decided from SSC
%   SMR: Signal to Mask Ratio, (42x8) converted to (336x1) for "ESH"
%                              (69x1) converted to (336x1) for "OLS",
%                              "LSS", "LPS" frameTypes
%   B219: the B219 table
%
% return value:
%   T : the thresholds, (42x8) converted to (336x1) for "ESH"
%                       (69x1) converted to (336x1) for "OLS",
%                        "LSS", "LPS" frameTypes
    
    if frameType == "ESH"
        B219b = B219;
        %split frameF into subframes
        X = buffer(frameF, 128);

        %split SMR into subframes
        SMR = buffer(SMR, length(B219b));
        
        %Calculate T for every subframe
        T_t = zeros(length(B219b),size(X,2));
        P = zeros(length(B219b),size(X,2));
        for c = 1:size(X,2) %1:8
            for b=1:length(B219)     %1:128
                P(b,c) = sum((X(B219b(b,2):B219b(b,3),c)).^2);
                T_t(b,c) = P(b,c) / SMR(b,c);
            end
        end
        
        T = zeros(42*8,1);
        for c = 1:size(T_t,2)
            pos = (c-1)*length(B219b);
            T(pos+1:pos+length(B219b)) = T_t(:,c);
        end
        
    else
        B219a = B219;
        %discard empty part of SMR
        SMR = SMR(1:length(B219a));
        
        %calculate T
        P = zeros(length(B219),1);
        T = zeros(length(B219),1);
        for b=1:length(B219)     
            P(b) = sum(frameF(B219(b,2):B219(b,3)).^2);
            T(b) = P(b)/SMR(b);
        end
        
        T = [T ; zeros(42*8 - 69,1)];
    end
end