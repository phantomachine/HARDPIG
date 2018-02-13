function [ pival_k ] = punishk_bilinear(self, k, PolyTriIndex, PolyLcons, ...
                                 V_MIN, V_MAX, IndexProfileSet,H, C, P_a)

%% PUNISHK.M
%
% Let $D = \cup_{k} Q_{k}$ be a 2D unit-simplex. Let $j,k = 1,...,K$.
% Calculate Max-Min punishment value for government player at state-space
% partition element Q(k).
%
% INPUT: 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% self   :  Model class from KSMOD03           
%
% k      :  index for partition element Q(k) of D.
%
% PolyTriIndex{a}{k}(i) : action profile a => P(a)
%                         Current state partition k <=> Q(k)
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
% V_MIN             : Natural upper bound on average discounted payoffs
%
% V_MAX             : Natural lower bound on average discounted payoffs
%
% IndexProfileSet : [ ia, ib ] index to action profiles (a,b) \in A x B
%
% H       :  Set of subgradients (L x NP) array
%
% C       :  Levels at direction h(l) and at partition k: H*w <= C(:,k).            
%            An (L x K) array
%
% P_a     :  Set of P(a) Markov maps; (N_Z x N_Z x NA) array
%
% OUTPUT:
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% pival_k   : scalar double giving the max-min payoff to large
%                    player at current state-space partition Q(k). This is
%                    $\pi(k)$ in the paper.
%                     
% (c) 2012, 2013, Timothy Kam. Email: tcy.kam__at__gmail.com 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Use subject to GNU LGPL licensing terms. Cite this header and 
% author completely in subsequent re-use and modifications.
% =======================================================================
% $Revision: 5.0.0 $  $Date: 2013/09/11 12:45:20 $ 
%
% See also PUNISH_OUTER, SOLVESDP (YALMIP)

%% Extract Inputs
    %self = model;
    %PolyTriIndex = xTriIndex; 
    %PolyLcons = xPolyLcons;
     
    % Number of agent classes
    NZ = size(H,2);
    
    % Number of -j states
    maxN = self.N;
    
    % Bounds on total payoffs:
    %Vmin = repmat(V_MIN, NZ, 1);
    %Vmax = repmat(V_MAX, NZ, 1);
    
    % Number of action profiles (a,b) \in A x B:
    Num_ab_Profiles = size(IndexProfileSet,1); 

    % Extract set of private agent action profiles A:
    A = self.ProfileSetA; 
    NA = self.NA;

    % Extract set of government action profiles B:
    B = self.ProfileSetB;
    NB = self.NB;
    
    % Wage parameter in utility function for agents type j > 0:
    wage = self.WAGE;
    
    % Discount factor:
    delta = self.DELTA;
    delstar = self.DELSTAR; % 1 - delta

%% Construct Punishment Value at each k \mapsto Q_{k}

    pival = zeros(Num_ab_Profiles,1);
    
    %% Main loop over all action profiles (a,b):
    %tic
     
    fprintf('\t\tPUNISHK.M: Optimizing on partition k = %i\n', k);

    parfor index_ab = 1 : Num_ab_Profiles

        % NOTE: Local data for each processor slice ("worker")
        
        %% 1. Parameters
        
        N = maxN;       % Max number of -j states
        NZ_local = NZ;   % Number private agents states, card(Z)
        NA_local = NA;  % Number of private action profiles a
        A_local = A;    % Set of private action profiles a
        B_local = B;    % Set of government action profiles b
        
        delta_local = delta;        % discount factor
        delstar_local = delstar;    % 1 - delta
        wage_local = wage;          % wage level
        
        P_a_local = P_a;    % N x N x card(A) set of Markov matrices
        
        C_local = C;
        PolyTriIndex_local = PolyTriIndex;
        
        Vmin = V_MIN;
        Vmax = V_MAX;
        
        %% 2. Action Profiles (a,b)
        
        % Index to action profiles (ia, ib) \mapsto (a,b) \in A x B:
        IndexProfileSet_local = IndexProfileSet; % Product space: (ia, ib)
                                                 
        ia = IndexProfileSet_local(index_ab,1);  
        ib = IndexProfileSet_local(index_ab,2);
           
        % Current P(a) matrix under a:
        P = P_a_local(:,:,ia);
        P_ia = repmat(P,NA_local-1,1);
        
        % Current action profile (a,b) - "parameters" here
        a = A_local(ia,:);             % row vector
        b = B_local(ib,:);             % row vector
        
        % All other a# \neq a, and P(a#) matrices:
        ia_Not = setdiff(IndexProfileSet_local(:,1), ia); % All a# \neq a
        P_ia_Not = P_a_local(:,:,ia_Not);   % All P(a#) \neq P(a), 
                                            % N x N x card(A)-1
        P_ia_not = reshape(P_ia_Not, NZ_local*(NA_local-1), NZ_local);
        
        %% 3. Current payoffs from profile (a,b)
        
        % Current payoff from b:
        cons = b;
        cons(N+1:end) = cons(N+1:end) + wage_local;
        Util_b = ucons(self, cons);     
        
        % Current payoff from a:
        Util_a = -phia(self,a); 
        Util_ia = repmat(Util_a, NA_local-1, 1);
        
        % Current payoffs from all deviations from a:
        Util_ia_Not = -phia(self, A_local(ia_Not,:) );
        
        % Total current payoff:        
        Util_now = delstar_local*(Util_a + Util_b); % row vector
        
        % Agents' incentive constraints:
        DevP = delta_local*(P_ia_not - P_ia);
        DevA = delstar_local*(Util_ia - Util_ia_Not);
        [nr, nc ] = size(DevA);
        DevA = reshape(DevA', nr*nc, 1); 
        
        %% 4. Set up as YALMIP Separable Bilinear Programs
        
        % Check # partitions intersect with current P(a)(Q_k), (a,k):
        PolyLcons_local = PolyLcons;
        J = numel(PolyLcons_local{ia}{k});
        
        % Storage for payoff values:
        pival_temp = zeros(J,1);
        
        for j = 1:J
            
            % Polygon feasible set for l*P(a): i.e. Aj * P(a)'l' <= dj
            Aj = PolyLcons_local{ia}{k}{j}(:,1:end-1);
            dj = PolyLcons_local{ia}{k}{j}(:,end);
            
            % Index to Q(j) that intersects with P(a)(Q(k)):
            % Use later in H*w <= C(:,kp)
            
            kp = PolyTriIndex_local{ia}{k}(j);
            C_kp =  C_local(:,kp);
            
            % Initialize YALMIP variables
            w = sdpvar(NZ_local,1);      % (N x 1) promised utilities
            l = sdpvar(1,NZ_local-1);    % (1 x N) probability distribution
            
            % Auxiliary variable
            %lp = P' * l';               % next-period distribution l+
            l = [ l, 1-sum(l)];
            
            % Define YALMIP objective function -- bilinear in (l,w)
        
            objective = l*Util_now' + delta_local*l*P*w; % bilinear form

            % Define YALMIP constraint set:
            %
            % (1) Consistency w.r.t. W: $w \in W(Q_j) :=: W(:,k')$
            % (2) Feasibility of $\lambda' \in P(a)(Q_k) \cap Q_j$
            % (3) Government busget constraint $\lambda b \geq 0$
            % (4) Agents' Competitive best response:
            %       $\delta [P(a#)-P(a)] w \leq [ phi(a#)-P(a) ]$
            %     for every $a# \neq a$.
            % (5) l(:) must be probabilities -- bounded optimizers
            % (6) w(:) have natural bounds   -- bounded optimizers
            
            Pl = P'*l';
            Pl = Pl(1:end-1);
                    
            constraint = [ H*w <= C_kp,...     % (1) L constraints
                           Aj*Pl <= dj,...  % (2) Variable # constraints
                           -l * b' == 0, ...   % (3) 1 constraint
                           DevP*w <= DevA,...  % (4) NZ*(card(A)-1) constrs
                           0 <= l(1:end-1) <= 1, ... % (5) l(:) prob distribution
                           Vmin <= w(:) <= Vmax ]; % (6) Bounded payoffs

            sdp_options = sdpsettings(  'verbose',0,...
                                        'cachesolvers',1,...
                                        'solver','bmibnb',...
                                        'bmibnb.roottight',0|1,...
                                        'bmibnb.numglobal', 20,...
                                        'bmibnb.uppersolver','snopt',...
                                        'bmibnb.lowersolver','glpk');

            % Minimum government payoff at (k,a,b) w.r.t. (l,w)
            sol = solvesdp(constraint,objective,sdp_options);
            
            % Throw solution status (0 = successful; 3 = max-iter):
            %yalmiperror(sol.problem)


            %if sol.problem == 0 || sol.problem == 3
            val = double(objective);
            
            %if ~isnan(val) && ~isinf(val)
                %disp('(w,l) should have an optimal value')
                %wstar = double(w);
                %lstar = double(l);
                %optval = lstar*Util_now' + delta_local*lstar*P*wstar;
                %punishval_j = optval;
                %punishval_j = double(objective);
						  
			    pival_temp(j) = val;
			
            %else
            %elseif sol.problem ~= 0 && sol.problem ~= 3
                %disp([k, index_ab, j])
                %disp(yalmiperror(sol.problem))
                %disp('Solver not found: (w,l) is not optimized')
                %double(w)
                %double(l)
                %punishval_j = +Inf;  % +Inf for MIN problem below
			%	pival_temp(j) = Vmax; % +Inf; TEST: 25/12/2013
		    %end
                       
            %pival_temp(j) = punishval_j;
        end
        
        % At each (a,b),minimize over j :=: (l,w)
        pival(index_ab) = min(pival_temp);        
        
    end
    %toc
    
    %% Construct Max-Min payoff w.r.t. (b,a) at current partition Q_k:
    % Reshape PIVAL:
    % Note: PIVAL := { pi(ia,ib) } 
    
    pival_k = reshape(pival, NA, NB);
    
    % Minimize over a:
    pival_k = min(pival_k,[],1);
    
    % Maximize over b:
    %pival_k(pival_k == Inf) = -Inf; % -Inf for MAX problem below
    pival_k = max(pival_k,[],2);
    
    % Final trap for errors before exit:
    %if ~isreal(pival_k) && ( isnan(pival_k) || isinf(pival_k) )
    %    disp(['Partition k =', k ])
    %    error('Optimization Problem Not Well Defined')
    %end
    
    % Clear garbage: avoid memory leaks
    % yalmip('clear')
    
