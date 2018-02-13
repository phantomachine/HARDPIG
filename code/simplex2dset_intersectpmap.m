function [xQP, xTriIndex, xPolyVerts, xPolyLcons, xPolyRand, xPolyQjRand, ...
            xPolyIndexQjRand ] = ...
               simplex2dset_intersectpmap(NA, P_a, dt, DtriVerts, Nsim, ...
                                            show_waitbar, varargin)
                    
%% SIMPLEX2DSET_INTERSECTPMAP
%
% PURPOSE.  Notation: 
%   a:  is (1 x N_Z) vector of action profile of agents; row element of
%   action profile set.
%   P(a):   is Markov map under each $a \in A$; a (N_Z x N_Z) nonsingular
%   array.
%   Q_{k}:  
% Function computes intersection of each triangular image $P(a)(Q_{k})$ 
% with every partition element $Q_{j}, j=1,...,K$ of unit simplex 
% $D := \cup_{k = 1,...,K} Q_{k}$. Intersections vary with each $(a, k)$.
%
% INPUT:
% 
% NA:         Number of action profiles $a \in A$
% P_a:        (N_Z x N_Z x NA) array of stacked $P(a)$
% dt:         Object from DelaunayTri in SIMPLEX2DSET_TRIPART.M ( set D )
% DtriVerts:  3D vectors of vertices of partition elements of D
% Nsim     :  Number of random uniform points in each PolyVerts
% show_waitbar: 'on' or 'off'
%
% OUTPUT:
% 
% xQP:        All intersections in Structure Array: vertices, indexes
% xTriIndex:  All intersections stored as indexes $j = 1,...,K$ of $Q(j)$
% xPolyVerts: All intersections stored as vertices of polygons, Poly(j)
% xPolyLcons: All intersections stored as linear-weak inequalities: Poly(j)
%
% USAGE:
%
% xPolyVerts{a}{k}{i} : action profile a => P(a) where a = 1,...,NA (index)
%                      Current state partition k <=> Q(k)
%                      i => j(i) <=> Q(j) intersects Q(k)
%
%                       Note: 
%                       * Polygon xPolyVerts{a}{k}{i} in triangle Q(j), 
%                         with number j := index{a}{k}(i); and j = 1,...,K.
% 
%                       * Index i = 1,...,N(a,k); where ...
%                         N(a,k) := numel(xPolyVerts{a}{k}) <= K.
%
% xPolyLcons{a}{k}{i} : is linear (weak) inequality constraints 
%                      polygon vertices from xPolyVerts{a}{k}{i}
%                      Used: function VERT2LCON by Matt Jacobsen (2011).
%                      Note:
%                      * At each {a}{k}{i}, AB = xPolyLcons{a}{k}{i} is an 
%                        array such that:
%                           A = AB(:,1:3);
%                           b = AB(:,4);
%                        and polygon stored as vertices xPolyVerts{a}{k}{i}
%                        is equivalently described by the set of all point
%                        $\lambda \in D \subset R^n$ such that:
%
%                           A * \lambda <= b.
% 
% See also DELAUNAYTRI, SIMPLEX2DSET_TRANSIT, SIMPLEX2DSET_INTERSECT
% % =======================================================================
% %     (c) 2013-- T.Kam
% %
% % DEPENDENCIES:
% %
% %     * KSMOD3 class for model and SSE solution.
% %         
% % Email: tcy.kam__at__gmail.com
% % =======================================================================
% % $Revision: 4.0.3 $  $Date: 2011/03/08 11:38:20 $
    
    if nargin == 6
        hit_and_run = [];
    elseif nargin == 7
        hit_and_run = varargin{1};
    end
    
    % if NA <= MAXIT display GUI WaitBar:
    MAXIT = 500;
    
    % Initialize:
    xQP = struct([]);       % All intersections, each P(a)
    xTriIndex = cell(NA,1); % All intersections' indexes, each P(a)
    xPolyVerts = cell(NA,1);% All intersections' vertices,each P(a)
    xPolyLcons = cell(NA,1);% All vertices as ineq. constraints
    xPolyRand = cell(NA,1); % Array of uniform r.v.'s in each cell
    xPolyQjRand = cell(NA,1); % NSIM Random vectors in partition Qj
    xPolyIndexQjRand = cell(NA,1); % Random vectors in partition Qj 
                                   % that intersect with Qk for each
                                   % P(a)(Qj)
    
    %% Main Loop: Pick each linear Markov operator P(a), then:
    %   1. simplex2dset_transit: computes all images P(a)(Q(k)), for every 
    %      (a,k) where a = 1,...,NA and k = 1,...,K.
    %   2. simplex2dset_intersect: computes are intersections of all 
    %      P(a)(Q(k)), at each (a,k), with all partition elements Q(j), 
    %      where j = 1,...,K.
 
  
if NA <= MAXIT && strcmp(show_waitbar, 'on')
    h = waitbar(0,'SIMPLEX2D-INTERSECTPMAP.M: Please wait...');
end

for n = 1:NA

    % Current P(a) matrix under action profile a := a(n):
    P = P_a(:,:,n);

    % Map partitions of D under current P(a) operator:
    PoD = simplex2dset_transit(P, DtriVerts);  
                                           % Row := each vertex' 3D 
                                           %        coordinates

    % Intersection(P(a)(D), D) and record them:                                       
    [ xTriIndex{n},xPolyVerts{n},xPolyLcons{n}, ...
        xPolyRand{n}, xPolyIndexQjRand{n}, xQP{n}, xPolyQjRand{n} ] = ...
              simplex2dset_intersect(dt.Triangulation, ...
                                     DtriVerts, PoD, P, Nsim, hit_and_run);

      if NA <= MAXIT && strcmp(show_waitbar, 'on')
            waitbar(n/NA,h,['SIMPLEX2D-INTERSECTPMAP.M: ',...
                                    num2str(n*100/NA), '% done'])
      end

end

if NA <= MAXIT && strcmp(show_waitbar, 'on')
    close(h)
end