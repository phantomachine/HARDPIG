function pival = punish_outer_lpset( self, PolyTriIndex, Fset, lb, ub,...
                                                      vartype,H, C, P_a)

%% PUNISH_OUTER.M
%
% Let $D = \cup_{k} Q_{k}$ be a 2D unit-simplex. Let $j,k = 1,...,K$.
% Calculate Max-Min punishment value for government player for ALL state
% space partition elements Q(k), given current estimate of correspondence
% stored in (H, C).
%
% INPUT: 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% self   :  Model class from KSMOD03           
%
% xTriIndex{a}{k}(i)  : action profile a => P(a)
%                       Current state partition k <=> Q(k)
%                           i => j(i) <=> Q(j) intersects Q(k)
%
%                       Note: 
%                       * Polygon PolyVerts{a}{k}{i} in triangle Q(j), 
%                         with number j := PolyTriIndex{a}{k}(i); and 
%                         j = 1,...,K.
% 
%                       * Index i = 1,...,N(a,k); where ...
%                         N(a,k) := numel(PolyVerts{a}{k}) <= K.
%
% xPolyLcons         : is linear (weak) inequality constraints 
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
% V_MIN             : Natural upper bound on average discounted payoffs
%
% V_MAX             : Natural lower bound on average discounted payoffs
%
% IndexProfileSet   : [ ia, ib ] index to action profiles (a,b) \in A x B
%
% H                 : Set of subgradients (L x NP) array
%
% C                 : Levels at direction h(l) and at partition k: 
%                           H*w <= C(:,k).            
%                     An (L x K) array
%
% P_a               : Set of P(a) Markov maps; (N_Z x N_Z x NA) array.
%
% OUTPUT:
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% pival   :          (1 x K) double array of max-min payoffs to large
%                    player at each current state-space partition Q(k). 
%                    This is $\{ \pi(k) | k = 1,...,K \}$ in the paper.
%                     
% (c) 2012, 2013, Timothy Kam. Email: tcy.kam__at__gmail.com 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Use subject to GNU LGPL licensing terms. Cite this header and 
% author completely in subsequent re-use and modifications.
% =========================================================================
% $Revision: 5.0.0 $  $Date: 2013/09/11 12:45:20 $
%
% See also PUNISHK
    
    %% Preliminaries:
    
    K = size(C,2);
    
    %% Solve bilinear programs at each state-space partition element Q(k):
    
    pival = zeros(1,K);
        
    % Let k be index of partition elements Q(k) containing next period's
    % state: $\lambda_{+} \in Q(k)$.
    %
    % Loop over all partition elements Q(k):
    
    fprintf('\tPUNISH_OUTER.M: ');
        
        % Restart MATLABPOOL
        matlabpool

        for k = 1:K
           %pival_k_temp = punishk( self, k, xTriIndex, xPolyLcons, ...
           %                   V_MIN, V_MAX, IndexProfileSet,H, C, P_a );
           
           pival_k_temp  = punishk_lpset(self, k, PolyTriIndex, Fset,...
                                            lb, ub, vartype, H, C, P_a);
           if isnan(pival_k_temp)
               pival(k) = V_MIN;
           else
               pival(k) = pival_k_temp;
           end
        end
        
        % Garbage cleanup for parallel workers:
        delete(gcp)

    fprintf('\tPUNISH_OUTER.M: Done ...\n');

end
