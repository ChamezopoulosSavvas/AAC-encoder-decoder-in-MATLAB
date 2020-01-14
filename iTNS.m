function frameFout = iTNS(frameFin, frameType, TNScoeffs)

% function frameFout = iTNS(frameFin, frameType, TNScoeffs)
% 
% reverses the TNS process
% arguments:
%   frameFin: the MDCT coefficients after TNS for a single channel 
%             array 1024x1 if the frame is of type "OLS", "LPS", "LSS"
%             array 128x8 converted to 1024x1 if the frame is of type "ESH"
%   frameType: the type of the frame as chosen at the SSC stage
%   TNScoeffs: the quantized TNS coefficients
%              array 32x1 if frame type is "OLS", "LSS", "LPS" 
%              only the first 4 cells are non-zero.
%              array 4x8 converted to 32x1 if frame type is "ESH"
%
% return value:
%   frameFout:the MDCT coefficients before TNS for a single channel
%             array 1024x1 if the frame is of type "OLS", "LPS", "LSS"
%             array 128x8 converted to 1024x1 if the frame is of type "ESH".

    if (frameType == "OLS" || frameType == "LSS" || frameType == "LPS")
        
        %take only the first 4 elements scince these are the coeffs we need
        TNScoeffs = TNScoeffs(1:4);
        %create and apply inverse filter
        frameFout = filter(1,[1, - TNScoeffs'],frameFin);
        %frameFout = real(filter(1,[1, TNScoeffs'],frameFin));

    elseif frameType == "ESH"
        
        %buffer frameFin and the TNS coefficients
        X_TNS = buffer(frameFin, 128);
        X = zeros(size(X_TNS));
        TNScoeffs = buffer(TNScoeffs, 4);
        frameFout = zeros(size(frameFin));
        
        for c=1:size(X,2) %1:8
            %create and apply filter
            %X(:,c) = X_TNS(:,c);
            X(:,c) = filter(1,[1, - TNScoeffs(:,c)'],X_TNS(:,c));
            %X(:,c) = real(filter(1,[1, TNScoeffs(:,c)'],X_TNS(:,c)));
            %store it in return value
            pos = (c-1)*size(X,1);
            frameFout(pos+1:pos+size(X,1)) = X(:,c);
            
        end
    end
end