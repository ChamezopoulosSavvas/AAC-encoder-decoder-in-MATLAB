function [S, sfc, G] = AACquantizer(frameF, frameType, SMR)

% function [S, sfc, G] = AACquantizer(frameF, frameType, SMR)
% This function implements the quantization process
% arguments:
%   frameF: 1024x1 array with the frame in MDCT coefficients AFTER TNS
%   frameType: the type of the frame as decided from SSC
%   SMR: the Signal to Mask ratio as produced from psycho
% return values:
%   S: 1024x1 array with the quantization symbols
%   sfc: ScaleFactor coefficients: (42x8) converted to (336x1) for "ESH"
%                                  (69x1) converted to (336x1) for "OLS",
%                                   "LSS", "LPS" frameTypes
%   G: Global gain of current frame: 1x8 for "ESH" ,scalar converted to 1x8
%                                    for "OLS","LSS", "LPS" frameTypes

    %% load tables
    vars = {'B219a','B219b'};
    load('Tableb219.mat', vars{:});
    %incr first 3 columns by 1 so they match matlab array indexes
    B219a(:,1:3) = B219a(:,1:3) + 1; %#ok<NODEF>
    B219b(:,1:3) = B219b(:,1:3) + 1; %#ok<NODEF>

    X = frameF;
    MQ = 8191;
    MagicNumber = 0.4054;

    if frameType == "ESH"
        
        X = buffer(X,128);
        
        %Calculate T for every subframe
        T = threshold(frameF , frameType, SMR, B219b);
        
        %split T into subframes
        T = buffer(T, length(B219b));
        
        %% step 1.1
        a = ones(length(B219b),8);
        Ss = zeros(length(X),8);
        X_est = Ss;
        P_e = zeros(length(B219b),8);
        P_e_prev = P_e;
        
        %% step 1.2
            
        
        for c=1:size(X,2)
            
            a(:,c) = ((16/3)*log2((max(X(:,c)).^(3/4))/MQ));
            flag = zeros(length(B219b),1);
            anythingMoved = 1;
            while(anythingMoved == 1)
                
                anythingMoved = 0;
                for b=1:length(B219b)
                    %calculate S, X_est, P_e
                    P_e_prev(b,c) = P_e(b,c);
                    P_e(b,c) = 0;
                    for k=B219b(b,2):B219b(b,3)
                        %calculate S
                        Ss(k,c) = sign(X(k,c))*floor((abs(X(k,c))*(2^(-a(b,c)/4)))^(3/4) + MagicNumber);
                        %calculate X_est
                        X_est(k,c) = (sign(Ss(k,c)))*(abs(Ss(k,c))^(4/3))*(2^(a(b,c)/4));
                        %calculate P_e
                        P_e(b,c) = P_e(b,c) + (X(k,c) - X_est(k,c))^2;
                    end
                    
                    %check for maxdiff=60
                    flag(b) = 0;
                    if (b < length(B219b) && b >1 )
                        if (abs(a(b+1,c) - a(b,c)) > 60) || (abs(a(b-1,c) - a(b,c)) > 60)
                            flag(b) = 1;
                        end
                    elseif (b == length(B219b) )
                        if (abs(a(b-1,c) - a(b,c)) > 60) 
                            flag(b) = 1;
                        end
                    elseif (b == 1 )
                        if (abs(a(b+1,c) - a(b,c)) > 60) 
                            flag(b) = 1;
                        end
                    end
                    
                    %&& (P_e_prev(b,c) - P_e(b,c) ~= 0)
                    if (P_e(b,c) < T(b,c) && flag(b) == 0 && (P_e_prev(b,c) - P_e(b,c) ~= 0))
                        a(b,c) = a(b,c) + 1;
                        anythingMoved = 1;
                    end
                end
            end
        end
        
        %% step 2.1
        G(1:8) = a(1,:);
        
        %% step 2.2
        sfcc = zeros(length(B219b),8);
        sfcc(1,:) = a(1,:);
        for i=2:length(B219b)
            sfcc(i,:) = a(i,:) - a(i-1,:);
        end
        
        %% unbuffer
        sfc = zeros(42*8,1);
        S = zeros(1024,1);
        for c=1:size(X,2)
            pos1 = (c-1)*size(X,1);
            pos2 = (c-1)*size(sfcc,1);
            S(pos1+1:pos1+size(X,1)) = Ss(:,c);
            sfc(pos2+1:pos2+length(B219b)) = sfcc(:,c); 
            
        end
    else

        T = threshold(frameF , frameType, SMR, B219a);
        T = T(1:69);
       
        %% step 1.1
        a = ones(length(B219a),1)*((16/3)*log2((max(X)^(3/4))/MQ));

        S = zeros(size(X));
        X_est = S;
        P_e = zeros(length(B219a),1);
        P_e_prev = P_e;
        %% step 1.2
        
        flag = zeros(length(B219a),1);
        anythingMoved = 1;
        counter = 1;
        while(anythingMoved == 1)
            
            anythingMoved = 0;
            for b=1:length(B219a)
                %calculate S, X_est, P_e
                P_e_prev(b) = P_e(b);
                P_e(b) = 0;
                for k=B219a(b,2):B219a(b,3)
                    %calculate S
                    S(k) = sign(X(k))*floor((abs(X(k))*(2^(-a(b)/4)))^(3/4) + MagicNumber);
                    %calculate X_est
                    X_est(k) = (sign(S(k)))*(abs(S(k))^(4/3))*(2^(a(b)/4));
                    %calculate P_e
                    P_e(b) = P_e(b) + (X(k) - X_est(k))^2;
                end
                
                flag(b) = 0;
                if (b < length(B219a) && b >1 )
                    if (abs(a(b+1) - a(b)) > 60) || (abs(a(b-1) - a(b)) > 60)
                        flag(b) = 1;
                    end
                elseif (b == length(B219a) )
                    if (abs(a(b-1) - a(b)) > 60) 
                        flag(b) = 1;
                    end
                elseif (b == 1 )
                    if (abs(a(b+1) - a(b)) > 60) 
                        flag(b) = 1;
                    end
                end
                
                %if (P_e(b) < T(b) && (P_e_prev(b) - P_e(b) ~= 0) && flag(b) == 0)
                if (P_e(b) < T(b) && flag(b) == 0 && (P_e_prev(b) - P_e(b) ~= 0))
                    a(b) = a(b) + 1;
                    anythingMoved = 1;
                end                
            end
            
            counter = counter + 1;
        end  
        
        
        %% step 2.1
        G(1:8) = a(1);
        
        %% step 2.2
        sfc = zeros(length(B219a),1);
        sfc(1) = a(1);
        for i=2:length(B219a)
            sfc(i) = a(i) - a(i-1);
        end
        sfc = [sfc ; zeros(length(B219b)*8 - length(B219a),1)];
        
    end


end