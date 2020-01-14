function a_q = quantizer(a,b)
% function a_q = quantizer(a,b,maxx)
% quantizes a vector a with
% a uniform, summetrical quantizer of step 0.1
% also, the function calculates whether or not the quantized coefficients
% will lead to an unstable inverse FIR filter, and adjust the coefficients
% accordingly, before they are quantized again.
%
% arguments:
%   a: the linear prediction coefficients
%   b: the number of bits of the quantizer
% return value:
%   a_q: the quantized coefficients

    %number of levels of the quantizer
    N = 2^b;
    delta = 0.1;
    
    d = (-N*delta/2:delta:N*delta/2)';
    
    %limits of zones
    d(1) = -inf;
    d(end) = inf;
    %quantization levels
    x = (-N*delta/2:delta:N*delta/2)';
    
    a_q = a;
    m_r = 2;
    while(m_r >=1)
        for i=1:length(a_q)
            for j = 1:length(d)
                if a_q(i) <= d(j)
                    a_q(i) = x(j);
                    break;
                end
            end
        end

        %check if coefficients are gonna lead to stable inverse fir filter
        
        % denominator
        den = [1, - a_q'];
        %abs of roots of denominator
        r_abs = abs(roots(den));
        %find root with max abs
        m_r = max(r_abs);
        %stablize polynomina, if necessary
        if m_r >= 1 
            den = polystab(den);
            a_q = den(2:end)';
        end
        %and re-quantize the coeffs
     end
end