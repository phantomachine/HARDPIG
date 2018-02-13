function Cnew = admit_outer(self, xTriIndex, xPolyLcons, V_MIN, V_MAX, ...
                                         IndexProfileSet, H, C, P_a,...
                                         pival, epsi, iter)

%%
% ADMIT_OUTER.M:
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
% pival             : (1 x K) array of punishment values.
%
% OUTPUT:
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Cnew   :          (L x K) double array containing updates of C. Recall:
%                   (H,C) represents extreme points of W correspondence.  
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

    [L, K] = size(C);   % No. of subgradients (L) and partition elements (K)
    NZ = size(H,2);     % Number of agent classes
    maxN = self.N;      % Max. number of -j agent states
    
    % Extract set of private agent action profiles A:
    A = self.ProfileSetA; 
    NA = self.NA;
    
    % Extract set of government action profiles B:
    B = self.ProfileSetB;
    
    % Number of action profiles (a,b) \in A x B:
    Num_ab_Profiles = size(IndexProfileSet,1); 
    
    % Wage parameter in utility function for agents type j > 0:
    wage = self.WAGE;
    
    % Discount factor:
    delta = self.DELTA;
    delstar = self.DELSTAR; % 1 - delta

    % Preallocate memory:
    Cnew = C;
    cplus_ab_temp = zeros(Num_ab_Profiles,1);
    %cplus_l_temp = zeros(L,1);
    
%% Main Loop: Admissibility w.r.t. (C)
    
    fprintf('\tADMIT_OUTER.M: ');
    
    for k = 1 : K               % Loop: partition elements Q(k)
       
         matlabpool 
         fprintf('\t iter = %i | k = %i | epsi(iter-1) = %6.3g\n',...
                                                            iter, k, epsi);
          %fprintf('\t\tIn direction l = \n');


        parfor l = 1 : L        % loop: subgradients h(l)
            
            %fprintf('\t\t\t< %i >\n',l);
            
            % Local worker copies of data:
            
            %k_local = k;
            
            N = maxN;       % Max. number of -j states
            NZ_local = NZ;   % Number private agents states, card(Z)
            NA_local = NA;  % Number of private action profiles a
        
            delta_local = delta;        % discount factor
            delstar_local = delstar;    % 1 - delta
            wage_local = wage;          % wage level
            
            C_local = C;
            H_local = H;
            h = H_local(l,:);   % Current subgradient
            A_local = A;        % Set of private action profiles a
            B_local = B;        % Set of government action profiles b
            
            IndexProfileSet_local = IndexProfileSet;
                                % Indices (ia, ib) to all (a,b)
                                
            P_a_local = P_a;    % N x N x card(A) set of Markov matrices
            PolyLcons_local = xPolyLcons;
            PolyTriIndex_local = xTriIndex;
            
            % Natural bounds on total average discounted payoffs:
            Vmin = V_MIN;
            Vmax = V_MAX;
            
            % (K x 1) vector of punishment values from PUNISH_OUTER.M:
            pival_local = pival;
            
            % Preallocate:
            cplus_ab_temp_local = cplus_ab_temp;
            %cplus_l_temp_local = cplus_l_temp;
            
            for index_ab = 1 : Num_ab_Profiles % Loop: profiles (a,b)
                
                %fprintf('\t\tAction profile No. %i of %i\n', index_ab, ...
                %                                        Num_ab_Profiles);
                
                %% 2. Action Profiles (a,b)

                ia = IndexProfileSet_local(index_ab,1);  
                ib = IndexProfileSet_local(index_ab,2);

                % Current P(a) matrix under a:
                P = P_a_local(:,:,ia);
                P_ia = repmat(P,NA_local-1,1);

                % Current action profile (a,b) - "parameters" here
                a = A_local(ia,:);             % row vector
                b = B_local(ib,:);             % row vector

                % All other a# \neq a, and induced P(a#) matrices:
                ia_Not = setdiff(IndexProfileSet_local(:,1), ia); 
                                                    % All a# \neq a
                P_ia_Not = P_a_local(:,:,ia_Not);   % All P(a#) \neq P(a), 
                                                    % N x N x card(A)-1
                P_ia_not = reshape(P_ia_Not,NZ_local*(NA_local-1),NZ_local);

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

                % Check # partitions intersect with current P(a)(Q_k):
                J = numel(PolyLcons_local{ia}{k});

                % Storage for payoff values:
                cplus_j_temp = zeros(J,1);

                for j = 1:J

                    % Polygon feasible set for l*P(a): Aj * P(a)'l' <= dj
                    Aj = PolyLcons_local{ia}{k}{j}(:,1:end-1);
                    dj = PolyLcons_local{ia}{k}{j}(:,end);

                    % Index to Q(j) that intersects with P(a)(Q(k)):
                    % Use later in H*w <= C(:,kp)
                    kp = PolyTriIndex_local{ia}{k}(j);
                    C_kp =  C_local(:,kp);
                    pival_k = pival_local(k);
                    
                    % TEST!!!
                    rep_pival_k = repmat(pival_k,NZ_local,1);
                    
                    % Initialize YALMIP variables
                    w = sdpvar(NZ_local,1);    % (N x 1) promised utilities
                    lam = sdpvar(1,NZ_local-1); % (1 x N) probability distn
                    lam = [lam, 1-sum(lam)];
                    
                    % Define YALMIP objective function -- bilinear in (l,w)
                    Evec = Util_now' + delta_local*P*w; 
                                          % (N x 1) vector of agent payoffs
                    objective = h*( Evec );
                                          % ... weighted by directions h
                    
                    
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
                    % (7) Government incentive compatibility

                    % TEST! 27/12/2013: removed agent optimality?
                    Pl = P'*lam';
                    Pl = Pl(1:end-1);
                    
                    constraint = [ H*w <= C_kp,...           % (1)
                                   Aj*Pl <= dj,...           % (2)
                                   -lam * b' <= 0, ...       % (3)
                                   0 <= lam(1:end-1) <= 1, ...     % (5)
                                   DevP*w <= DevA,...        % (4)
                                   Vmin <= w(:) <= Vmax, ... % (6)
                                  -Evec <= -rep_pival_k];    % (7)

                    %constraint = [ H*w <= C_kp,...           % (1)
                    %               Aj*P'*lam' <= dj,...      % (2)
                    %               -lam * b' <= 0, ...       % (3)
                    %               DevP*w <= DevA,...        % (4)
                    %               0 <= lam(:) <= 1, ...     % (5)
                    %               Vmin <= w(:) <= Vmax, ... % (6)
                    %               pival_k <= lam*Evec ];    % (7)

                    sdp_options = sdpsettings('verbose',0,...
                                              'cachesolvers',1,...
                                              'solver','bmibnb',...
                                              'bmibnb.roottight',0|1,...
                                              'bmibnb.numglobal', 20,...
                                              'bmibnb.uppersolver',...
                                                               'snopt',...
                                              'bmibnb.lowersolver','glpk');

                    % Maximum w.r.t (lam, w):
                    sol = solvesdp(constraint,-objective,sdp_options);

                    % Throw solution status (0 = successful; 3 = max-iter):
                    %yalmiperror(sol.problem)

                    %if sol.problem == 0 || sol.problem == 3
                     %if sol.problem == 0 || sol.problem == 3
                    val = double(objective);
            
                    %if ~isnan(val) && ~isinf(val) 
                    %if isnan(val)
                    %    cplus_j_temp(j) = -Inf;

                    %elseif sol.problem ~= 0 && sol.problem ~= 3
                    %else
                        cplus_j_temp(j) = val;
                    %end

                    %cplus_j_temp(j) = cval_j;   % c+[k,l,(a,b),j]
                    
                    
                end % EndFor j
                
                % max over all j:=: j(a,b):
                cplus_ab_temp_local(index_ab) = max(cplus_j_temp); 
                                                        % c+[k,l,(a,b)]

            end % EndFor index_ab
            
            cp_temp =  max(cplus_ab_temp_local);

            % Ad-hoc Error trap:
            %if isinf(cp_temp)
           %     Cnew(l,k) = Vmin; %c+(l,k)
           % else
                Cnew(l,k) = cp_temp;
           % end
                    
        end % EndFor l
        
        % Clear garbage: avoid memory leaks
        %yalmip('clear')
        matlabpool close
        

    end % EndFor k

    fprintf('\t\tADMIT_OUTER.M: Done!\n');
    
end
