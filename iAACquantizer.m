function frameF = iAACquantizer(S, sfc, G, frameType)

%function frameF = iAACquantizer(S, sfc, G, frameType)
% 
% This function reverses the quantization step
% arguments:
%   S: 1024x1 array with the quantization symbols
%   sfc: ScaleFactor coefficients: (42x8) converted to (336x1) for "ESH"
%                                  (69x1) converted to (336x1) for "OLS",
%                                   "LSS", "LPS" frameTypes
%   G: Global gain of current frame: 1x8 for "ESH" ,scalar converted to 1x8
%                                    for "OLS","LSS", "LPS" frameTypes
%   frameType: the type of the frame as decided from SSC
% return value:
%   frameF: 1024x1 array with the frame in MDCT coefficients BEFORE iTNS

    vars = {'B219a','B219b'};
    load('Tableb219.mat', vars{:}); 
    %incr first 3 columns by 1 so they match matlab array indexes
    B219a(:,1:3) = B219a(:,1:3) + 1; %#ok<NODEF>
    B219b(:,1:3) = B219b(:,1:3) + 1; %#ok<NODEF>
    frameF = zeros(1024,1);
        
    if frameType == "ESH"
    
        % buffer S, sfc
        S = buffer(S,128);
        sfc = buffer(sfc, length(B219b));
        a = zeros(size(sfc));
        X_est = zeros(size(S));
        
        for c=1:size(sfc,2)
            % acquire a_j
            a(1,c) = G(1,c);
            for i=2:size(a,1)
                a(i,c) = sfc(i,c) + a(i-1,c);
            end
            % calculate X_est
            for i=1:size(a,1)
                 X_est((B219b(i,2)):(B219b(i,3)),c) = ...
                     (sign(S((B219b(i,2)):(B219b(i,3)),c))).* ...
                     (abs(S((B219b(i,2)):(B219b(i,3)),c)).^(4/3)).*...
                     (2.^(a(i,c)/4));   
            end
        end
        
        % un-buffer
        for c=1:size(S,2)
            pos = (c-1)*size(S,1);
            frameF(pos+1:pos+size(S,1)) = X_est(:,c);
        end
    else
        sfc = sfc(1:length(B219a));
        G = G(1);
        a = zeros(size(sfc));
        frameF = zeros(size(S));
        
        %acquire a_j
        a(1) = G;
        for i=2:length(a)
            a(i) = sfc(i) + a(i-1);
        end
        
        % acquire frameF
        for i=1:size(a,1)
                 frameF((B219a(i,2)):(B219a(i,3))) = ...
                     (sign(S((B219a(i,2)):(B219a(i,3))))).* ...
                     (abs(S((B219a(i,2)):(B219a(i,3)))).^(4/3)).*...
                     (2.^(a(i)/4));   
        end
        
        frameF = frameF';
    end
end