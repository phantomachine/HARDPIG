function [ index, PolyVerts, PolyLcon, PolyRand, ...
                PolyIndexQjRand, IntersectQP, PolyQjRand  ] = ...
                 simplex2dset_intersect(trivert, Q, PoD, P, Nsim, varargin)

%% SIMPLEX2DSET_INTERSECT.M
%
% Let $D = \cup_{k} Q_{k}$ be a 2D unit-simplex. Let $j,k = 1,...,K$.
% This function finds intersections of all j triangles $P \circ Q_j$,
% with all triangular partition elements $Q_{k} \in D$.
% 
% Usage will involve each $P(Q_{k})$, where $P$ is defined in the function
% SIMPLEX2DSET_TRANSIT.M
%
% INPUT: Warning ... Run SIMPLEX2DSET_TRIPART to get object dt beforehand!
% 
% trivert:  Indices of vertices of Delaunay triangular partition elements.
%           Run SIMPLEX2DSET_TRIPART beforehand: trivert = dt.Triangulation 
%
% Q      :  K x 3 array of vertex coordinates of triangles $Q_k$. 
%           Run SIMPLEX2DSET_TRIPART beforehand: Q = [ dt.X, 1-sum(dt.X)]
%
% PoD    :  Run SIMPLEX2DSET_TRANSIT to calculate image of $P(a) circ D$,
%           where $P(a)$ is Markov matrix that depends on action parameters
%           $a$.
% 
% P      :  Short for current Markov matrix P(a) given profile a
% 
% Nsim   :  Number of random uniform points in each PolyVerts
%
% OUTPUT:
%
% IntersectQP:   Cell array with K cell elements. Usage:
%                IntersectQP{j} : j-th cell; contains intersections of
%                                  $P \circ Q_j$, with all $Q_k \in D$.
%                                  e.g. m-th intersection of $Q_j$ with 
%                                  m = 1,...,M (and M is nonconstant)
%                                  is a polygon $QP(m; j,k)$ with:
%
%                                  (1) IntersectQP{j}(m).IndexQk : 
%                                  Index k, of partition element $Q_k$
%                                  containing polygon $QP(m; j,k)$; 
% 
%                                  (2) IntersectQP{j}(m).VertsQk : 
%                                  2D coordinates of vertices of polygon 
%                                  $QP(m; j,k)$.
%
% index     :   Indexes of partition elements Q_{k} that intersect with 
%               every $P \circ Q_j$
%
% USAGE:
% PolyVerts{a}{k}{i} : action profile a => P(a)
%                      Current state partition k <=> Q(k)
%                      i => j(i) <=> Q(j) intersects Q(k)
%
%                       Note: 
%                       * Polygon PolyVerts{a}{k}{i} in triangle Q(j), 
%                         with number j := index{a}{k}(i); and j = 1,...,K.
% 
%                       * Index i = 1,...,N(a,k); where ...
%                         N(a,k) := numel(PolyVerts{a}{k}) <= K.
%
% PolyLcons{a}{k}{i} : is linear (weak) inequality constraints 
%                      polygon vertices from IntersectPolyVerts{a}{k}{i}
%                      Uses function VERT2LCON
%                      Note:
%                      * At each {a}{k}{i}, AB = PolyLcons{a}{k}{i} is an 
%                        array such that:
%                           A = AB(:,1:3);
%                           b = AB(:,4);
%                        and polygon stored as vertices PolyVerts{a}{k}{i}
%                        is equivalently described by the set of all point
%                        $\lambda \in D \subset R^n$ such that:
%
%                           A * \lambda <= b.
%                     
% (c) 2012, 2013, Timothy Kam. Email: tcy.kam__at__gmail.com 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Use subject to GNU LGPL licensing terms. Cite this header and 
% author completely in subsequent re-use and modifications.
% =======================================================================
% $Revision: 5.0.0 $  $Date: 2013/09/11 12:45:20 $ 
%
% See also SIMPLEX2DSET_TRANSIT, SIMPLEX2DSET_TRIPART, POLYBOOL, VERT2LCON
    
    if nargin < 6
        hit_and_run = 0;
    elseif nargin == 6
        hit_and_run = varargin{1};
    end
    
    %% Define triangular-polygon: simple-cycle over vertices
    cycle = [ 1, 2, 3, 1 ];
    
    %% Number of partition elements in dt
    K = size(trivert,1);

    %% Preallocate cell memory address:
    
    % j loop             
    index = cell(K,1);       % Indexes Qk intersection with current Qj
    IntersectQP = index;     % Both index and polyx together
    PolyVerts = index;       % Vertices Qk intersection with current Qj
    PolyLcon = index;        % Linear inequality representation of Vertices
    PolyRand = index;        % Uniform r.v.s in each intersecting polygon
    PolyQjRand = index;      % Store 
    PolyIndexQjRand = index;
    
    % k loop
    % Output:        #1         #2      #3          #4      #5
    field_names = {'IndexQk','VertsQk','LconsQk','RandQk','IndexQjRand'};
    empty_cells = repmat(cell(1),1,numel(field_names));
    entries = {field_names{:} ; empty_cells{:}};
    s = struct(entries{:});
    IntersectQkPart_j = repmat(s,K,1);                
    
    %% Loop over each partition element $Q_j$, j = 1,...,K:
    for j = 1:K % Current states

        % Pull out current k-th triangular partition $Q_j$:
        vert_index = trivert(j,:);         % vertices index
        Qj = Q(vert_index,:);              % Current-state partition $Q_j$
        Triangle = Qj(cycle,:);            % Simple walk over $Qj$ vertices
        
        % Also simulate random vectors (uniform) on 2D projection of $Q_j$:
        
        Triangle = Triangle(:,1:end-1); % 2D projection of 3D coordinates
        
        if hit_and_run
            [A,b,Aeq,beq] = vert2lcon(Triangle);
            At = [ A; Aeq ];
            bt = [ b; beq ];
            RandUpoints = cprnd(Nsim,At,bt);
        else
            % Own function: RANDPOLYFILL
            RandUpoints = randpolyfill(Triangle,Nsim,0,1);
        end
        
        % Index to each row vector in RandUpoints
        IndexRand = 1:size(RandUpoints,1);
        
        % Pull back up to 3D
        X = [RandUpoints, 1 - sum(RandUpoints,2)];
        PolyQjRand{j} = X;
        
        % Apply Markov map $P(a)$ to partition $Qj$ ...
        % Map current Qj under mapping P, TPoQj = Qj*P, pre-stored as PoD.
        % Map image is also a triangle:
        TPoQj = PoD(vert_index,:);                  % index-to-3D-coords
                                         
        xj = TPoQj(cycle,1);                        % 2D coordinates
        yj = TPoQj(cycle,2);
        
        % Now apply map P to every point in X as well ...
        Xplus = X*P;        
        
        % Initialize k-counter: for non-empty intersections(Qj,Qk)
        intersect_index = 0;    
        
        % Loop over each $Q_k$:
        for k = 1:K % Next period states

            % Pull out current k-th triangular partition $Q_k$:
            vert_index_k = trivert(k,:); % vertex indices of $Q_k$
            Qk = Q(vert_index_k,:); % Vertex indices to Cartesian coords
            
            % Current Qk coordinates
            xk = Qk(cycle,1);
            yk = Qk(cycle,2);
            %[xk,yk] = poly2cw(xk,yk); % Clockwise reordering: vertices

            % Intersect each $Q_k$ with currently fixed $P \circ Q_j$:
            [xi, yi] = polybool('intersection', xj, yj, xk, yk);
            
            % Store non-empty Intersections(Qj,Qk) ...
            % NOTE: variable max(#Intersections) stored in intersect_index:
            if ~isempty(xi) && ~isempty(yi)
                
                % 1. Index of partition element containing this polygon:
                intersect_index = intersect_index + 1;
                IntersectQkPart_j(intersect_index).IndexQk = k;
                
                % 2. Polygon vertices:
                %verts = [ xi, yi, 1-xi-yi ];   % 3D original
                verts = [ xi, yi ];             % 2D projection
                IntersectQkPart_j(intersect_index).VertsQk = verts;
                
                % 3. Convert verts to linear (weak) inequalities:
                [A,b,Aeq,beq] = vert2lcon(verts);
                Atemp = [ A; Aeq ];
                btemp = [ b; beq ];
                Ab = [ Atemp, btemp ];
                IntersectQkPart_j(intersect_index).LconsQk = Ab;
                
                % 4. Check if uniform random points in each polytope/gon:
                [in, on] = inpolygon(Xplus(:,1),Xplus(:,2),xi,yi);
                IntersectQkPart_j(intersect_index).RandQk = ...
                                  [ Xplus(in,1), Xplus(in,2) ];
                 
                % 5. Indexes to row numbers of X(j) s.t. P(a)(X(j)) \cap Q(k)
                IntersectQkPart_j(intersect_index).IndexQjRand = ...
                              sort([IndexRand(on == 1), IndexRand(in==1)]);
                                  
            end 

        end

        %% OUTPUT #1 ...
        % Current $P(a) \circ Q_j$ intersections with ALL $Q_k$ (i.e. $D$):
        IntersectQP{j} = IntersectQkPart_j(1:intersect_index);
        
        
        % Collapse down to a cell array:
        indtmp = struct2cell(IntersectQP{j});
        
        %% OUTPUT #1 ...
        % First field 'IndexQk':
        index{j} = cell2mat(indtmp(1,:));
        
        %% OUTPUT #2 ...
        % Second field 'VertsQk':
        verts_tmp = indtmp(2,:);
        PolyVerts{j} = verts_tmp;
        
        %% OUTPUT #3 ...
        % Third field 'LconsQk':
        PolyLcon{j} = indtmp(3,:);
        
        %% OUTPUT #4 ...
        % Fourth field 'PolyQjRand':
        PolyRand{j} = indtmp(4,:); 
        
        %% OUTPUT #5 ...
        % Fifth field 'PolyIndexQjRand':
        PolyIndexQjRand{j} = indtmp(5,:);
       
    end

    
%% MORE NOTES:
% -------------------------------------------------------------------------
% dt     :      Object from DelaunayTri's triangular partitions;
%               where:
%               * dt.X   :  (NX x 2) array of 2D triangle vertices.
%                           NX < K*3, is number of unique/common partition
%                           element vertices (in 2D coordinates).
%               * dt.Triangulation  :  (K x 3) array of indices to dt.X
%
% Q is equivalent to [ dt.X, 1 - sum(dt.X,2) ];
