function m = smooth_min(x,y,epsilon)
% SMOOTH_MIN(x,y,epsilon) berechnet das Minimum von x und y, falls
% |x-y|>=epsilon, oder benutzt smooth_ppart, um die min-Funktion zu
% gl?tten, falls |x-y|<epsilon.
% 
% Input:
%       x - n by m matrix      
%       y - n by m matrix
%       epsilon - pos. scalar << 1            
%
% Output:
%       m - n by m matrix           value of the smooth min(x,y)

if all(size(x) ~= size(y)), error('size(x) unequal size(y)'); end
[m, ~, ~] = smooth_ppart(x-y, epsilon);
m = -m - x;