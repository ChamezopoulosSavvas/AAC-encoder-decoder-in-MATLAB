function SMR = psycho(frameT, frameType, frameTprev1, frameTprev2)

% This function implements the Psychoacoustic model for a single channel
% arguments:
%   frameT: the frame in time domain (2048x1)
%   frameType: the type of the frame as decided from SSC
%   frameTprev1: the previous frame in time domain (2048x1)
%   frameTprev2: the previous frame of frameTprev1 in time domain (2048x1)
%
% return value:
%   SMR: Signal to Mask Ratio, (42x8) converted to (336x1) for "ESH"
%                              (69x1) converted to (336x1) for "OLS",
%                              "LSS", "LPS" frameTypes
    %% step 0
    vars = {'B219a','B219b'};
    load('Tableb219.mat', vars{:}); 
    %incr first 3 columns by 1 so they match matlab array indexes
    B219a(:,1:3) = B219a(:,1:3) + 1; %#ok<NODEF>
    B219b(:,1:3) = B219b(:,1:3) + 1; %#ok<NODEF>
    
    %% step 1
    %load appropriate spread sheet
    vars = {'sprda','sprdb'};
    load('spreads.mat',vars{:});
        
    if frameType == "ESH"
        %% step 2
        %we only need the last 2 subframes from the previous frames
        % process the previous frame
        frameTprev1_w = hannWin(frameTprev1, frameType); 
        % process the frame before previous frame
        %frameTprev2_w = hannWin(frameTprev2, frameType);
            
        %buffer the previous frames
        subframesTprev1_w = buffer(frameTprev1_w(end - 512 + 1:end),256);
        %subframesTprev2_w = buffer(frameTprev2_w,256);
        
        %subframesFprev1_w = zeros(length(subframesTprev1_w), 1);
        %subframesFprev2_w = zeros(length(subframesTprev1_w), 1);
        r_1 = zeros(length(subframesTprev1_w), 8);
        f_1 = r_1; r_2 = r_1; f_2 = r_1;
        
        %prepare frameT
        frameT_w = hannWin(frameT, frameType);
        subframesT_w = buffer(frameT_w, 256);
        subframesF_w = zeros(size(subframesT_w,1), size(subframesT_w,2));
        r_pred = zeros(size(subframesT_w,1), size(subframesT_w,2));
        f_pred = r_pred; r = r_pred; f = r_pred; c = r_pred;
        
        %init for step 5
        e = zeros(length(B219b),size(subframesT_w,2));
        c_b = zeros(length(B219b),size(subframesT_w,2));
        
        %init for step 6
        ecb = zeros(length(B219b),size(subframesT_w,2));
        ct = zeros(length(B219b),size(subframesT_w,2));
        cb = ct;
        en = ct;
        
        %init for step 7
        tb = cb;
        
        % init for step 8
        SNR = tb;
        
        % init for step 9
        bc = SNR;
        
        % init for step 10
        nb = bc;
        
        % init for step 11
        qthr = nb;
        npart = qthr;
        
        % init for step 12
        SMRc = npart;
        
        for cl=1:size(subframesT_w,2)
            
            %% step 2
            subframesF_w(:,cl) = fft(subframesT_w(:,cl));
            r(:,cl) = abs(subframesF_w(:,cl));
            f(:,cl) = angle(subframesF_w(:,cl));
            
            if cl == 1
                
                subframesFprev1_w = fft(subframesTprev1_w(:,2));
                subframesFprev2_w = fft(subframesTprev1_w(:,1));
        
                r_1(:,cl) = abs(subframesFprev1_w);
                r_2(:,cl) = abs(subframesFprev2_w);
            
                f_1(:,cl) = angle(subframesFprev1_w);
                f_2(:,cl) = angle(subframesFprev2_w);
                
                %% step 3
                r_pred(:,cl) = 2*r_1(:,cl) - r_2(:,cl);
                f_pred(:,cl) = 2*f_1(:,cl) - f_2(:,cl);
            
            elseif cl == 2
                
                %subframesFprev1_w = fft(subframesTprev1_w(:,2));
                subframesFprev2_w = fft(subframesTprev1_w(:,2));
        
                r_1(:,cl) = r(:,cl-1);
                r_2(:,cl) = abs(subframesFprev2_w);
            
                f_1(:,cl) = f(:,cl-1);
                f_2(:,cl) = angle(subframesFprev2_w);
                
                %% step 3
                r_pred(:,cl) = 2*r_1(:,cl) - r_2(:,cl);
                f_pred(:,cl) = 2*f_1(:,cl) - f_2(:,cl);
            
            else
            
                %pred the rest of the subframes
                %% step 3
                r_pred(:,cl) = 2*r_pred(:,cl-1) - r_pred(:,cl-2);
                f_pred(:,cl) = 2*f_pred(:,cl-1) - f_pred(:,cl-2);
            
            end
            
            %% step 4
            c(:,cl) = (sqrt(((r(:,cl).*(cos(f(:,cl)))) - (r_pred(:,cl).*(cos(f_pred(:,cl))))).^2 ...
                + ((r(:,cl).*(sin(f(:,cl)))) - (r_pred(:,cl).*(sin(f_pred(:,cl))))).^2))./...
                (r(:,cl) + abs(r_pred(:,cl)));
            
            
            %% step 5
            for b=1:length(B219b)
                e(b,cl) = sum(r(B219b(b,2):B219b(b,3),cl).^2);
                c_b(b,cl) = sum(c(B219b(b,2):B219b(b,3),cl).*(r(B219b(b,2):B219b(b,3),cl).^2));            
            end
            
            %% step 6
            for b=1:length(B219b)
                den = 0;
                for bb=1:length(B219b)
                    ecb(b,cl) = ecb(b,cl) + e(bb,cl)*sprdb(bb,b);
                    ct(b,cl) = ct(b,cl) + c_b(bb,cl)*sprdb(bb,b);
                    den = den + sprdb(bb,b);
                end
                cb(b,cl) = ct(b,cl)/ecb(b,cl);
                en(b,cl) = ecb(b,cl)/den;
            end
            
            %% step 7
            for b=1:length(B219b)
                tb(b,cl) = -0.299 - 0.43*log(cb(b,cl));
            end
            
            %% step 8
            NMT = 6;
            TMN = 18;
            for b=1:length(B219b)
                SNR(b,cl) = TMN*tb(b,cl) + (1 - tb(b,cl))*NMT;
            end
            
            %% step 9
            for b=1:length(B219b)
                bc(b,cl) = 10^(-SNR(b,cl)/10);
            end
            
            %% step 10
            for b=1:length(B219b)
                nb(b,cl) = en(b,cl)*bc(b,cl);
            end
            
            %% step 11
            for b=1:length(B219b)
                qthr(b,cl) = eps*(2048/2)*10^(B219b(b,6)/10);
                npart(b,cl) = max(nb(b,cl), qthr(b,cl));
            end
           
            %% step 12
            for b=1:length(B219b)
                SMRc(b,cl) = e(b,cl)/npart(b,cl);
            end
        end 
        
        %% final step: format SMR
        SMR = zeros(42*8,1);
        for c = 1:size(SMRc,2)
            pos = (c-1)*length(B219b);
            SMR(pos+1:pos+length(B219b)) = SMRc(:,c);
        end
        
    else
        %% step 2
        % process the previous frame
        frameTprev1_w = hannWin(frameTprev1, frameType);
        frameFprev1_w = fft(frameTprev1_w);
        %find angle and phase
        r_1 = abs(frameFprev1_w);
        f_1 = angle(frameFprev1_w);
        
        %process the frame before the previous frame
        frameTprev2_w = hannWin(frameTprev2, frameType);
        frameFprev2_w = fft(frameTprev2_w);
        %find angle and phase
        r_2 = abs(frameFprev2_w);
        f_2 = angle(frameFprev2_w);
        
        %process current frame
        frameT_w = hannWin(frameT, frameType);
        frameF_w = fft(frameT_w);
        %find angle and phase
        r = abs(frameF_w);
        f = angle(frameF_w);
        
        %% step 3
        r_pred = 2*r_1 - r_2;
        f_pred = 2*f_1 - f_2;
        
        %% step 4
        c = (sqrt(((r.*(cos(f))) - (r_pred.*(cos(f_pred)))).^2 ...
            + ((r.*(sin(f))) - (r_pred.*(sin(f_pred)))).^2))./...
            (r + abs(r_pred));
        
        %% step 5
        e = zeros(length(B219a),1);
        c_b = zeros(length(B219a),1);
        for b=1:length(B219a)
            e(b) = sum(r(B219a(b,2):B219a(b,3)).^2);
            c_b(b) = sum(c(B219a(b,2):B219a(b,3)).*(r(B219a(b,2):B219a(b,3)).^2));            
        end
        
        %% step 6
        ecb = zeros(length(B219a),1);
        ct = zeros(length(B219a),1);
        en = ct;
        for b=1:length(B219a)
            den = 0;
            for bb=1:length(B219a)
                ecb(b) = ecb(b) + e(bb)*sprda(bb,b);
                ct(b) = ct(b) + c_b(bb)*sprda(bb,b);
                den = den + sprda(bb,b);
            end
            en(b) = ecb(b)/den;
        end
        cb = ct./ecb;
        
        %% step 7
        tb = -0.299 - 0.43*log(cb);
        
        %% step 8
        NMT = 6;
        TMN = 18;
        
        SNR = TMN*tb + (1 - tb)*NMT;
        
        %% step 9
        bc = 10.^(-SNR/10);
        
        %% step 10
        nb = en.*bc;
        
        %% step 11
        qthr = eps*(2048/2)*10.^(B219a(:,6)/10);
        npart = max(nb, qthr);
        
        %% step 12
        SMRc = e./npart;
        
        %% final step: format SMR
        SMR = [SMRc ; zeros(42*8 - 69,1)];
    end
end