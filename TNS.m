function [frameFout, TNScoeffs] = TNS(frameFin, frameType)

% function [frameFout, TNScoeffs] = TNS(frameFin, frameType)
% This function implements the TNS stage for a single channel
% arguments:
%   frameFin: the MDCT coefficients BEFORE TNS for a single channel 
%             array 1024x1 if the frame is of type "OLS", "LPS", "LSS"
%             array 128x8 if the frame is of type "ESH"
%   frameType: the type of the frame as chosen at the SSC stage
% 
% return values:
%   frameFout: the MDCT coefficients AFTER TNS for a single channel
%              array 1024x1 if the frame is of type "OLS", "LPS", "LSS"
%              array 128x8 converted to 1024x1if the frame is of type "ESH"
%   TNScoeffs: the quantized TNS coefficients
%              array 32x1 if frame type is "OLS", "LSS", "LPS" 
%              only the first 4 cells are non-zero.
%              array 4x8 converted to 32x1 if frame type is "ESH"
    vars = {'B219a','B219b'};
    load('Tableb219.mat', vars{:}); 
    
    %incr first 3 columns by 1 so they match matlab array indexes
    B219a(:,1:3) = B219a(:,1:3) + 1; %#ok<NODEF>
    B219b(:,1:3) = B219b(:,1:3) + 1; %#ok<NODEF>
    
    X = frameFin;
    
    if frameType == "ESH" 
        
        %split frameFin into subframes
        X = buffer(X, 128);
        order = 4;
        a = zeros(order, size(X,2));            %4x8
        TNScoeffs1 = a;                         %4x8
        TNScoeffs = zeros(order*size(X,2),1);   %32x1
        
        %step 1
        P = zeros(length(B219b),size(X,2));
        S_w = zeros(size(X));
        
        for c = 1:size(X,2) %1:8
            for j=1:length(B219b)     %1:128
                
                P(j,c) = sum(X(B219b(j,2):B219b(j,3),c).^2);
                for k=B219b(j,2):B219b(j,3)
                    S_w(k,c) = sqrt(P(j,c));
                end
            end
            
            %S_w smoothing
            for k=(size(X,1)-1):-1:1
                S_w(k,c) = (S_w(k,c) + S_w(k+1,c))/2;
            end
            for k=2:size(X,1)
                S_w(k,c) = (S_w(k,c) + S_w(k-1,c))/2;
            end
        end
        
        %normalize X
        X_n = X./S_w;
       
        %step 2
        for c=1:size(X_n,2)
            r = autocorr(X_n(:,c),order);
            R = toeplitz(r(1:4), r(1:4));
            r = r(2:(order+1));
            a(:,c) = R\r;
            
            %step 3
            b = 4; %number of bits for quantizer
            %apply quantizer
            TNScoeffs1(:,c) = quantizer(a(:,c),b);
            
            %step 4 - fir filter on non-Normalized frameFin
            X(:,c) = filter([1, - TNScoeffs1(:,c)'],1,X(:,c)); 
            %X(:,c) = real(filter([1, TNScoeffs1(:,c)'],1,X(:,c))); 
           
        end
        
        frameFout = zeros(size(frameFin));
        for c=1:size(X,2) %1:8
            pos = (c-1)*size(X,1);
            frameFout(pos+1:pos+size(X,1)) = X(:,c);
            pos = (c-1)*size(TNScoeffs1,1);
            TNScoeffs(pos+1:pos+size(TNScoeffs1,1)) = TNScoeffs1(:,c);
        end
        
        
    else
        
       %step 1
        P = zeros(length(B219a),1);
        S_w = zeros(length(X),1);
        for j=1:length(B219a)     
            P(j) = sum(X(B219a(j,2):B219a(j,3),1).^2);
            for k=B219a(j,2):B219a(j,3)
                S_w(k) = sqrt(P(j));
            end
        end
        
        %S_w smoothing
        for k=1023:-1:1
            S_w(k) = (S_w(k) + S_w(k+1))/2;
        end
        for k=2:1024
            S_w(k) = (S_w(k) + S_w(k-1))/2;
        end
        %normalize X
        X_n = X./S_w;
        
        %step 2
        order = 4;
        r = autocorr(X_n,order);
        R = toeplitz(r(1:4), r(1:4));
        r = r(2:(order+1));
        a = R\r;
        
        %step 3
        b = 4; %number of bits for quantizer        
        
        %apply quantizer
        TNScoeffs = quantizer(a,b);
        
        %fir filter
        frameFout = filter([1, - TNScoeffs'],1,X);
        
        %frameFout = real(filter([1, TNScoeffs'],1,X));
        
        TNScoeffs = [TNScoeffs ; zeros(4*7,1)];
    
    end
end