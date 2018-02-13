function output = randpolyfill(xPolyVerts, Nsim, varargin)

% RandPolyFill.M
%
% Function computes realizations of uniform random variables
% constrained to be inside or on bondary of a given polygon.
% For any polygon summarized by the clockwise-ordered vertices
% in array xPolyVerts. Based on a simple pentagon example by Roger
% Stafford.
%
%
% INPUT:
% 
%   xPolyVerts:     Clockwise orderding; each row is a vertex.
%
%   Nsim:           Number of random points to simulate.
%
%   varargin:       Optional inputs:
%                   {1}. plotresults ( = 1: 'on'); (= 0: 'off')
%                   {2}. fixrandseed (= 1: fixed seed); (= 0 random)
% OUTPUT:
% 
%   output:         Uniformly distributed points in polygon.
%
% (c) 2013, Timothy Kam. Email: tcy.kam__at__gmail.com 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Use subject to GNU LGPL licensing terms. Cite this header and 
% author completely in subsequent re-use and modifications.
% =======================================================================
% $Revision: 5.0.0 $  $Date: 2013/09/11 12:45:20 $ 
%
% See also RAND, ABS, DET (Matlab in-built)

[ N_verts, N_dim ] = size(xPolyVerts);

if nargin == 4
    plotresults = varargin{1};
    fixrandseed = varargin{2};
elseif nargin == 3
    plotresults = varargin{1};
    fixrandseed = 0;
elseif nargin == 2
    plotresults = 0;
    fixrandseed = 0;
    disp('You can set third input as PLOTRESULTS = 1 to display output');
else
    error('RANDPOLYFILL:nargin','You need minimum TWO inputs.')  
end

if N_dim ~= 2
    % Throw error and quit if polygon not 2D:
    error('RANDPOLYFILL:N_dim','This only works for 2D polygons!')
end

% Number of triangular partitions:
NtriParts = N_verts - 2;

% Storage for triangular partition elements of polygon:
DoubleAreaTriangle = zeros(NtriParts);

    % Fix vertex #1 as origin:
    p1 = xPolyVerts(1,:);

    for index = 1 : NtriParts
        % Vectors (u and v) from p1 to vertices of current triangle
        CurrentTriVec = [ xPolyVerts(index+1,:) - p1;       % Vector u
                          xPolyVerts(index+2,:) - p1    ];  % Vector v

        % 2 x area of a triangle formula:
        % Area of triangle given by vectors u and v from an origin is
        % given by: (1/2)*| det (u v) | = (1/2)*| u x v |.
        % See e.g.: http://mathworld.wolfram.com/TriangleArea.html
        DoubleAreaTriangle(index) = abs( det( CurrentTriVec ) );
        
    end

    % Total area of polygon = total sum of area of triangles:
    DoubleAreaPolygon = sum(DoubleAreaTriangle);

    % Normalized cumulative area of triangle partition elements:
    NormCumArea = [0; cumsum(DoubleAreaTriangle)/DoubleAreaPolygon];
    
    % Simulate Nsim number of uniform r.v.'s on [0,1] support:
    if fixrandseed
        s = RandStream('mt19937ar','Seed',0);
        RandStream.setGlobalStream(s);
    end

    u = rand(Nsim,1);    % To select which triangle to fall into
    s = sqrt(rand(Nsim,1)); % Use s2 to select points within ...
        s2 = [s,s];

    t = rand(Nsim,1);    % ... given triangle ...
        t2 = [t, t];

    % Generalized vertices
    pa = zeros(Nsim, N_dim);
    pb = zeros(Nsim, N_dim);
    
    for index = 1 : NtriParts
        % pa = (r<=a)*p2 + ((a<r)&(r<=b))*p3 + (b<r)*p4; 
        % pb = (r<=a)*p3 + ((a<r)&(r<=b))*p4 + (b<r)*p5;
        u1 = NormCumArea(index);
        u2 = NormCumArea(index+1);
        pa = pa + (( u1 < u) & (u <= u2)) * xPolyVerts(index+1,:);
        pb = pb + (( u1 < u) & (u <= u2)) * xPolyVerts(index+2,:);
    end
   
    
    % Resulting random points in polygon:
    output = (1-s)*p1 + s2.*((1-t2).*pa + t2.*pb); % The random points
    
    % Plot results
    if plotresults
        c = xPolyVerts([1:N_verts,1],:);

        figure
            plot(c(:,1),c(:,2),'ro',c(:,1),c(:,2),'r-',...
                    output(:,1),output(:,2),'b*')
            axis equal
    end
end
