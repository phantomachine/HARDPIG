function [ dt, DT, K_DT, ic ] = simplex_delaunay(D,n)

% SIMPLEX_DELAUNAY.M: Generate Delaunay Triangulations
% (c) T.Kam, 2013
%
% INPUT:
%       D{n}    : Structure matrix of 3D vertex coordinates
%       n       : index of n-th partition scheme, partitions K = 4^n
% 
% OUTPUT:
%       dt      : Delaunay structure variable
%       DT      : 3D coordinates of Delaunay trinagular-complex vertices
%       K_DT    : numbers of unique points in Delaunay triangulation
%       ic      : in-centers of each Delaunay triangle partition element


% Generate Delaunay triangulation based on unique points in D:
    x = D{n}(:,1);
    y = D{n}(:,2);

    dt = DelaunayTri(x,y);

    % 3D coordinates of simplex:
    DT = [dt.X, 1 - sum(dt.X,2)];
    K_DT = size(dt.Triangulation,1);

    ic = incenters(dt);