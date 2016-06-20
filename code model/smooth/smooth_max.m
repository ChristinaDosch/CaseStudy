function [m, gradX, gradY, hessXX, hessYY, hessXY] = smooth_max(x,y,epsilon)
% SMOOTH_MAX(x,y,epsilon) berechnet das Maximum von x und y, falls
% |x-y|>=epsilon, oder benutzt smooth_ppart, um die max-Funktion zu
% gl?tten, falls |x-y|<epsilon.
% 
% Input:
%       x - n by m matrix      
%       y - n by m matrix
%       epsilon - pos. scalar << 1            
%
% Output:
%       m - n by m matrix           value of the smooth max(x,y)

if all(size(x) ~= size(y)), error('size(x) unequal size(y)'); end

[m, gradX, hessXX, ~] = smooth_ppart(x-y, epsilon);
m = m + y;
gradY = -gradX + 1;
hessYY = hessXX;
hessXY = -hessXX;

