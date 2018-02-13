function [status, l1, l2] = intriangle2(P, T)

% INTRIANGLE.M
% 
% PURPOSE: Point (p) in triangle (T) test.
%
% INPUTS:
%   p : (2 x 1) (double) array of (x,y) Catesian coordinates of a 
%       vector in R^3.
%
%   T : (3 x 2) (double) array of 2-simplex. Each row T(i,:) = 1,2,3 is an 
%       (x,y)-coordinate of triangle vertex vi.
%
% OUTPUT:
%   status : 0 (point p is outside triangle T)
%            1 (point p is in/on edge of triangle T)
%   l1     : Barycentric coordinate (x-basis) for p(1)
%   l2     : Barycentric coordinate (y-basis) for p(2)
%
% EXAMPLE USAGE:  
%       p = [ 0.5, 0.5] 
%       T = [ 1 0; 0 1; 0 0 ]; 
%   Then:
%       [status, l2, l2] = intriangle(p, T)
%   should return:
%       >>      status = 1 (i.e. point p is in triangle T)
%               l1 = 0.5
%               l2 = 0.5
%
% MORE INFORMATION:
%
% Cartesian coordinates of point p = (p1, p2) where p in unit 2-simplex; 
% Note: in 3D we have (p1, p2, p3) where p3 = 1-p1-p2.
%
% Barycentric coordinates lambda = (l1, l2) of point p; and lambda in
% unit 2-simplex as well.
%
% 2-simplex triangle vertices in (x,y) Cartesian coordinates: 
%   T = [    v1;
%            v2;
%            v3   ]
%
% where (x,y) Cartesian coordinates are:
% 
%   v1 = (x1, y1), z1 = 1-x1-y1
%   v2 = (x2, y2), z2 = 1-x2-y2
%   v3 = (x3, y3), z3 = 1-x3-y3
%
% Linear transform between Cartesian coordinates p and Barycentric
% Coordinates lambda: p = lambda * V; or more verbosely
%
%   p1 = l1*x1 + l2*x2 + l3*x3
%   p2 = l1*y1 + l2*y2 + l3*y3
%
% Since l3 = 1-l1-l2, we can write
%
%   p1 = l1*x1 + l2*x2 + (1-l1-l2)*x3
%   p2 = l1*y1 + l2*y2 + (1-l1-l2)*y3
%
% Or compactly, the linear transform is T such that
%
%   T * (l1, l2)' = (p1, p2)' - (x3, y3)', 
%
% so then solving for (l1, l2):
%
% (l1, l2)' = inv(T) * [ (p1, p2)' - (x3, y3)' ]
%
% THEOREM:
%   (i)  Point p strictly in triangle iff (l1 + l2 < 1), l1 > 0 and  l2 > 0
%   (ii) Point p on an edge of triangle if 0 <= li <= 1 and
%               l1 = 0, or l2 = 0 or l1 + l2 = 1
%   (ii) (Otherwise) Point p strictly outside triangle if:
%                   (l1 > 1) or (l2 > 1) or (l1 + l2) > 1    
%
% In this function (i) or (ii) yield OUTPUT = 1 (i.e. "in the triangle")
%
% (c) 2013, T. Kam (tcy.kam@gmail.com)
%
% See also INPOLYGON (MATLAB in-built)

% Cartesian coordinates of triangle vertices

x1 = T(1,1); y1 = T(1,2); % v1 = (x1, y1)
x2 = T(2,1); y2 = T(2,2); % v2 = (x2, y2)
x3 = T(3,1); y3 = T(3,2); % v3 = (x3, y3)

% Set up hard-coded matrix inv(T):
iT_11 = (y2 - y3)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2);
iT_12 = -(x2 - x3)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2);
iT_21 = -(y1 - y3)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2);
iT_22 = (x1 - x3)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2);

iT = [  iT_11, iT_12;
        iT_21, iT_22 ];
    
% Barycentric coordinates: l1 and l2 (Note: l3 = 1 - l1 - l2 )
% Note: l2 and l2 are basis vectors corresponding to point p.

l12 = iT * (p - T(3,:))';

l1 = l12(1);
l2 = l12(2);

if (l1 >= 0) && (l2 >= 0) && (l1 + l2 <= 1)
    status = 1;
else
    status = 0;
end