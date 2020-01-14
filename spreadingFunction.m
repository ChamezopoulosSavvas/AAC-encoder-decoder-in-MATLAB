function x = spreadingFunction(i,j,bval_i, bval_j)
% This func calculates the spread between two different bands
% arguments:
%   i : the band that spreads
%   j : the band that is spread in
%   bval_i: the central frequency of the i band
%   bval_j: the central frequency of the j band
% return value:
%   x: the spread between the bands

    if i>=j
        tmpx = 3*(bval_j - bval_i);
    else
        tmpx = 1.5*(bval_j - bval_i);
    end

    tmpz = 8*min((tmpx - 0.5)^2 - 2*(tmpx - 0.5),0);
    tmpy = 15.811389 + 7.5*(tmpx + 0.474) - ...
           17.5*sqrt(1 + (tmpx + 0.474)^2);
    if tmpy < -100
        x = 0;
    else
        x = 10^((tmpz+tmpy)/10);
    end
end
