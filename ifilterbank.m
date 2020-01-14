%
% frameT : frame under 
% frameType : type of frame as it is defined from sequence segmentation control
% winType : type of window that is selected

%OLS, LPS, LSS, ESH
function frameT = ifilterbank( frameF, frameType, winType)

% function frameT = ifilterbank( frameF, frameType, winType)
% 
% function that reverses the filterbank process
% arguments:
%   frameF: the transformed frame (1024x2)
%   frameType: the type of the frame (OLS, LSS, ESH, LPS)
%   winType: the type of the window (KBD, SIN)

    %% Initialize
    N = 2 * length(frameF(:,1));
    %total window for OLS,LSS,LPS
    w = zeros(N,1);
    %window for ESH, legth=256
    w_s = zeros(256,1);
       
    %% Create windows
    if frameType == "OLS"
        
        if winType == "KBD"
            a=6;
            %first: 1024 left w_wdb
            w(1:N/2,1) = w_KBD_half(N, a, "left");
            %at the end: 1024 right w_wdb
            w(N/2+1:N,1) = w_KBD_half(N, a, "right");
            
        elseif winType == "SIN"
            %first: 1024 left SIN
            w(1:N/2,1) = w_SIN_half(N, "left");
            %at the end: 1024 right SIN
            w(N/2+1:N,1) = w_SIN_half(N, "right");
        end
        
    elseif frameType == "LSS"
        
        if winType == "KBD"
            %first: 1024 left w_wdb
            a=6;
            w(1:N/2,1) = w_KBD_half(N, a, "left");
            %next: 448 ones
            w(N/2+1 : N/2 + 448,1) = 1;
            %next: 128 right w_wdb
            a=4;
            w(N/2+448+1 : N/2+448 + 128,1) = w_KBD_half(256, a, "right");
            %at the end: 448 zeros
            w(N/2+448+1+128 : N,1) = zeros(448,1);
            
        elseif winType == "SIN"
            %first: 1024 left SIN
            w(1:N/2,1) = w_SIN_half(N, "left");
            %next: 448 ones
            w(N/2+1 : N/2 + 448,1) = ones(448,1);%1;
            %next: 128 right SIN
            w(N/2+448+1 : N/2+448 + 128,1) = w_SIN_half(256, "right");
            %at the end: 448 zeros
            w(N/2+448+1+128 : N,1) = zeros(448,1);
        end     
        
    elseif frameType == "LPS"
        
        if winType == "KBD"
            %first: 448 zeros
            w(1:448,1) = zeros(448,1);
            %next: 128 left w_wdb
            a=4;
            w(448+1 : 448 + 128,1) = w_KBD_half(256, a, "left");
            %next: 448 ones
            w(448+128+1 : 448+128 + 448,1) = ones(448,1);%1;
            %at the end: 1024 right w_wdb
            a=6;
            w(448+128+448+1 : N,1) = w_KBD_half(N, a, "right");
            
        elseif winType == "SIN"
            %first: 448 zeros
            w(1:448,1) = zeros(448,1);
            %next: 128 left SIN
            w(448+1 : 448 + 128,1) = w_SIN_half(256, "left");
            %next: 448 ones
            w(448+128+1 : 448+128 + 448,1) = ones(448,1);%1;
            %at the end: 1024 right SIN
            w(448+128+448+1 : N,1) = w_SIN_half(N, "right");
        end    
        
    elseif frameType == "ESH"
        
        if winType == "KBD"
            a=4;      
            w_s(1:128,1) = w_KBD_half(256, a, "left");
            w_s(129:256,1) = w_KBD_half(256, a, "right");
            
        elseif winType == "SIN"
            w_s(1:128,1) = w_SIN_half(256, "left");
            w_s(129:256,1) = w_SIN_half(256, "right");           
        end
    end
    

    %% Calculate IMDCT coefficients
    
    if (frameType == "OLS" || frameType == "LSS" || frameType == "LPS")
        z = imdct4(frameF);
        %z = imdctl(frameF);
        
    elseif frameType == "ESH"
        z_sub_R = zeros(256,8);
        z_sub_L = zeros(256,8);
        %for both columns of frame frameF
        for i=1:8
            pos = (i-1)*128;
            z_sub_R(:,i) = imdct4(frameF(pos+1:pos+128,1));
            z_sub_L(:,i) = imdct4(frameF(pos+1:pos+128,2));
            %z_sub_R(:,i) = imdctl(frameF(pos+1:pos+128,1));
            %z_sub_L(:,i) = imdctl(frameF(pos+1:pos+128,2));
        end
    end
    
    %% Create frameT
    frameT = zeros(N,2);
    
    %After IMDCT we multiply our signal with the proper window depending...
    %   on the frametype to produce the frameT
    if (frameType == "OLS" || frameType == "LSS" || ...
            frameType == "LPS")
        %multiply with the window
        frameT = w .* z;
        
    elseif frameType == "ESH"
        %multiply each subframe with the window
        subframes_R = zeros(256,8);
        subframes_L = zeros(256,8);
        for i=1:8
            subframes_R(:,i) = w_s .* z_sub_R(:,i);
            subframes_L(:,i) = w_s .* z_sub_L(:,i);
            %subframes(:,i+8+1) = w_s .* z_sub(:,i+8+1);
        end
        
        %create frameT
        for i=1:8
            pos = 448+(i-1)*128;
            frameT(pos+1 : pos+256,1) = frameT(pos+1 : pos+256,1) + subframes_R(:,i) ;
            frameT(pos+1 : pos+256,2) = frameT(pos+1 : pos+256,2) + subframes_L(:,i) ;
        end
    end    
end

    
    
    
    
    