function frameF = filterbank( frameT, frameType, winType)
% function frameF = filterbank( frameT, frameType, winType)
% 
% function that implements the filterbank stage of encoding
% arguments:
%   frameT : the frame of the two channel signal in time domain.
%            array 2048x2
%   frameType : type of frame as it is defined from sequence
%               segmentation control. One of "OLS", "LPS", "LSS", "ESH"
%   winType : type of window that is selected

    %% initialize
    
    % number of samples
    N = length(frameT(:,1));
    %total window for OLS,LSS,LPS
    w = zeros(N,1);
    %window for ESH, legth=256
    w_s = zeros(256,1);
   
    %% create windows
    if frameType == "OLS"
        
        if winType == "KBD"
            a=6;
            %first: 1024 left w_wdb
            w(1:(N/2),1) = w_KBD_half(N, a, "left");
            %at the end: 1024 right w_wdb
            w((N/2+1):N,1) = w_KBD_half(N, a, "right");
            
        elseif winType == "SIN"
            %first: 1024 left SIN
            w(1:(N/2),1) = w_SIN_half(N, "left");
            %at the end: 1024 right SIN
            w((N/2+1):N,1) = w_SIN_half(N, "right");
        end
        
    elseif frameType == "LSS"
        
        if winType == "KBD"
            %first: 1024 left w_wdb
            a=6;
            w(1:(N/2),1) = w_KBD_half(N, a, "left");
            %next: 448 ones
            w((N/2+1) : (N/2 + 448),1) = ones(448,1);%1;
            %next: 128 right w_wdb
            a=4;
            w((N/2+448+1) : (N/2+448 + 128),1) = w_KBD_half(256, a, "right");
            %at the end: 448 zeros
            w((N/2+448+1+128) : N,1) = zeros(448,1);
            
        elseif winType == "SIN"
            %first: 1024 left SIN
            w(1:N/2,1) = w_SIN_half(N, "left");
            %next: 448 ones
            w((N/2+1) : (N/2 + 448),1) = ones(448,1);%1;
            %next: 128 right SIN
            w((N/2+448+1) : (N/2+448 + 128),1) = w_SIN_half(256, "right");
            %at the end: 448 zeros
            w((N/2+448+1+128) : N,1) = zeros(448,1);
        end     
        
    elseif frameType == "LPS"
        
        if winType == "KBD"
            %first: 448 zeros
            w(1:448,1) = zeros(448,1);
            %next: 128 left w_wdb
            a=4;
            w((448+1) : (448 + 128),1) = w_KBD_half(256, a, "left");
            %next: 448 ones
            w((448+128+1) : (448+128 + 448),1) = ones(448,1);%1;
            %at the end: 1024 right w_wdb
            a=6;
            w((448+128+448+1) : N,1) = w_KBD_half(N, a, "right");
            
        elseif winType == "SIN"
            %first: 448 zeros
            w(1:448,1) = zeros(448,1);
            %next: 128 left SIN
            w(448+1 : 448 + 128,1) = w_SIN_half(256, "left");
            %next: 448 ones
            w((448+128+1) : (448+128 + 448),1) = ones(448,1);%1;
            %at the end: 1024 right SIN
            w((448+128+448+1) : N,1) = w_SIN_half(N, "right");
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
    
    %% Calculate z
    if (frameType == "OLS" || frameType == "LSS" || frameType == "LPS")
        z = w .* frameT;
        
    elseif frameType == "ESH"
        
        %Create array subframes to store subrames for both chanels and put them on columns
        subframes_R = zeros(256,8);
        subframes_L = zeros(256,8);
        %get 8 subframes for each channel of frameT
        for i=1:8
            pos = 448+(i-1)*128;
            subframes_R(:,i) = frameT(pos+1 : pos+256,1);
            subframes_L(:,i) = frameT(pos+1 : pos+256,2);
        end
        
        %Create z_sub. z_sub takes the values of subframes multiplied
        %element-wise with the filter w_s
        z_sub_R = zeros(256,8);
        z_sub_L = zeros(256,8);
        for i=1:8
            z_sub_R(:,i) = w_s .* subframes_R(:,i);
            z_sub_L(:,i) = w_s .* subframes_L(:,i);
        end
    end
    
    %% Calculate MDCT coefficients
    
    frameF = zeros(N/2,2);
    
    if (frameType == "OLS" || frameType == "LSS" || frameType == "LPS")
        
        frameF = mdct4(z);
        %frameF = mdctl(z);   
        
    elseif frameType == "ESH"

        %For each window of both columns (channels) of frameT, calculate simultaneously mdct...
        %   coefficients and store them in frameF as it follows:
        %In first column (right channel) of frameF put one after another,...
        %   the coefficients for every suframe of the right channel.
        %   Do the same for the second column
        for i=1:8
            %call mdctl that returns the needed part of frameF
            pos=(i-1)*128;
            %right side
            frameF(pos+1:pos+128, 1) = mdct4(z_sub_R(:,i));
            %frameF(pos+1:pos+128, 1) = mdctl(z_sub_R(:,i));
            %left side
            frameF(pos+1:pos+128, 2) = mdct4(z_sub_L(:,i));
            %frameF(pos+1:pos+128, 2) = mdctl(z_sub_L(:,i));
        end
    end
end
    