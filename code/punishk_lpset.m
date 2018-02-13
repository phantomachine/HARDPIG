function [ pival_k ] = punishk_lpset(self, k, PolyTriIndex, Fset, lb, ub,...
                                                        vartype, H, C, P_a)

%% PUNISHK_LPSET.M
%
% Let $D = \cup_{k} Q_{k}$ be a 2D unit-simplex. Let $j,k = 1,...,K$.
% Calculate Max-Min punishment value for government player at state-space
% partition element Q(k). This version uses:
%
%   (1) Monte Carlo* simulation to evaluate maxima/minima at finite members
%   of the partition elements Q(k);
%   (2) Thus, (1) renders the approximate problem as a set of standard LP 
%   optimization problems.
%   (3) Since (1) convergence in probability to uniform distribution on
%   each convex polytope partition, the optima in simulation-based LPs: 
%   (1) and (2) should converge in probability to true global optimum--i.e.
%   outer approximation of punishment set.
%
%   * REFERENCES: 
%
%   (1) Lovasz, Laszlo (1999): "Hit and Run Mixes Fast", Math. Prog.
%       Ser. A. 86: 443-461.
%   (2) Smith, R.L. (1984): "Efficient Monte-Carlo Procedures for
%       Generating Points Uniformly Distributed Over Bounded Regions",
%       Oper. Res., 32: 1296-1308.

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
% PolyRand{a}{k}{i} : is linear (weak) inequality constraints 
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
% (c) 2013-, Timothy Kam. Email: tcy.kam__at__gmail.com 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Use subject to GNU LGPL licensing terms. Cite this header and 
% author completely in subsequent re-use and modifications.
% =======================================================================
% $Revision: 5.0.0 $  $Date: 2013/09/11 12:45:20 $ 
%
% See also PUNISH_OUTER

%% Preliminaries:

    % GLPK options
    param.msglev=1; % Output GLPK messages on workspace      
    param.save=1;   % Set save options
    
    % Settings:
    NZ = size(H,2);     % Number of agent classes
    maxN = self.N;      % Max. number of -j agent states
    
    % Extract set of private agent action profiles A:
    A = self.ProfileSetA; 
    NA = self.NA;
    
    % Extract set of government action profiles B:
    B = self.ProfileSetB;
    
    % Wage parameter in utility function for agents type j > 0:
    wage = self.WAGE;
    
    % Discount factor:
    delta = self.DELTA;
    delstar = self.DELSTAR; % 1 - delta

    
%% Main Loop: Admissibility w.r.t. (C)
fprintf('\t\tPUNISHK_LPSET.M: Optimizing on partition k = %i\n', k);    

    parfor index_a = 1 : NA        % loop: subgradients h(l)

%       fprintf('\t\t\t< %i >\n',l);

        % Local worker copies of data:
        
        N = maxN;       % Max. number of -j states
        NZ_local = NZ;   % Number private agents states, card(Z)
        NA_local = NA;  % Number of private action profiles a
        
        param_local = param;
        %ctype_local = ctype;
        vartype_local = vartype;
        lb_local = lb;
        ub_local = ub;

        delta_local = delta;        % discount factor
        delstar_local = delstar;    % 1 - delta
        wage_local = wage;          % wage level

        C_local = C;
        H_local = H;
        A_local = A;        % Set of private action profiles a
        B_local = B;        % Set of government action profiles b

        % Feasible sets F{a}{k}{i} where i := i(a,k)
        F = Fset;           % Feasible (lamda,b) sets: GBC >= 0

        % N x N x card(A) set of Markov matrices                               
        P_a_local = P_a;    
        PolyTriIndex_local = PolyTriIndex;

        % Natural bounds on total average discounted payoffs:
        % Vmin = V_MIN; % See lb
        % Vmax = V_MAX; % See ub


            % Current P(a) matrix under a:
            P = P_a_local(:,:,index_a);
            P_ia = repmat(P,NA_local-1,1); % Tile it NA-1 times!

            % Current action profile a:
            a = A_local(index_a,:);         % row vector
            %b = B_local(ib,:);             % row vector

            % All other a# \neq a, and induced P(a#) matrices:
            Index_A = 1:NA;                     % indexes to { a }
            ia_Not = setdiff(Index_A, index_a); 
                                                % All a# \neq a
            P_ia_Not = P_a_local(:,:,ia_Not);   % All P(a#) \neq P(a), 
                                                % N x N x card(A)-1
            P_ia_not = reshape(P_ia_Not,...
                                NZ_local*(NA_local-1),NZ_local);
                                            % Tile this too!
            
            %% Agent's best response:
            %
            % Current payoff from a: effort disutility
            Util_a = -phia(self,a); 
            Util_ia = repmat(Util_a, NA_local-1, 1);
            
            % Current payoffs from all deviations from a:
            Util_ia_Not = -phia(self, A_local(ia_Not,:) );
            
            % Agents' incentive constraints at current {k,l,a}
            DevP = delta_local*(P_ia_not - P_ia);
            DevA = delstar_local*(Util_ia - Util_ia_Not);
            [nr, nc ] = size(DevA);
            DevA = reshape(DevA', nr*nc, 1);
           

            %% 4. Solve feasible LP on each intersecting sub-state-space

            % Check # partitions intersect with current P(a)(Q_k):
            J = numel(F{index_a}{k}); % Each given P(a)[Q_k] intersects with 
                                 % J other partition elements Q_j

            % Storage for payoff values:
            pival_j_temp = zeros(J,1);
            
            j = 1;
            while j <= J

                % Recall F encodes GBC feasibility correspondence
                set_lambda = F{index_a}{k}{j}.Xi; 
                                          % Feasible realizations of states
                                          
                % Feasbile pairs (lambda, b)
                index_lambda = F{index_a}{k}{j}.ind_Xi;
                                          % Index to particular lambda 
                index_B = F{index_a}{k}{j}.ind_B;
                                          % Index to particular b
                balance_nonnegative = F{index_a}{k}{j}.balance;
                                          % Inner products <lambda,b> >= 0
                
                if isempty(balance_nonnegative) % Infeasible GBC set at j
                    
                    % Skip solving LP in intersection "j" since it has
                    % no feasible vectors (lambda,b) or no realizations
                    % drawn on that partition element
                    pival_j_temp(j) = -inf;
                    j = j + 1;
                    
                else % Feasible GBC set at j

                    % Current number of feasible pairs (lambda,b) at
                    % {k,l,a}
                    N_feasible = numel(index_lambda);
                    
                    % Catch GLPK errors in current j
                    status_vec = zeros(N_feasible,1);
                   
                    pival_n_temp = status_vec;

                    for n = 1:N_feasible
                        
                        %% Set up current payoff U(a,b)
                        % Current state lambda:
                        lambda = set_lambda( index_lambda(n), :); 

                        % Current payoff from b:
                        b_current_feasible = B_local(index_B(n),:);
                        cons = b_current_feasible;
                        cons(N+1:end) = cons(N+1:end) + wage_local;
                        Util_b = ucons(self, cons);  

                        % Total current payoff:        
                        Util_now = delstar_local*( Util_a + Util_b )'; 
                                                        % column vector

                        % Index to Q(j) that intersects with P(a)(Q(k)):
                        kp = PolyTriIndex_local{index_a}{k}(j);
                        
                        % Map kp to relevant step correspondence level
                        C_kp =  C_local(:,kp); %Use later in H*w <= C(:,kp)
                             
                        % Pick relevant punishment value at current state k
                        %pival_k = pival_local(k);
                        
                        %% Now set up LP for GLPK
                        %
                        % Define objective function for maximization:
                        %
                        %       objective' * w
                        %
                        % where actual (NZ x 1) array OBJECTIVE is:
                        %
                        %   objective = h* (Util_now' + delta_local*P) ...

                        
                       % Evec = Util_now + delta_local*P;
                        
                       
                        % .... but Util_now is a constant for GLPK so:
                        
                        conval = delta_local* P; 
                                            % h-weigted payoffs on Z
                                            % Add back Util_now' below
                        objective = (lambda * conval);
                        % Define constraint set in the form:
                        %
                        %       constraint           *   w   <= bconstraint
                        %   (|H|+(|A|-1)+1) x NZ      (NZx1)
                        %
                        % (0) Implicit:
                        %     (i) GBC: < lambda, b > >= 0. Encoded in F.
                        %     (ii)Feasible states, encoded in F:
                        %         $\lambda'\in P(a)(Q_k) \cap Q_j$.
                        %     (iii) w natural bounds. Encoded in (lb,ub).
                        % (1) Consistency: $w \in W(Q_j) :=: W(:,k')$
                        % (2) Agents' Competitive best response:
                        %     $\delta [P(a#)-P(a)] w \leq [ phi(a#)-P(a) ]$
                        %     for every $a# \neq a$.
                      
                        constraint = [  H_local;                ... (1)
                                        DevP  ];                ... (2)
                                        
                        bconstraint = [ C_kp;                   ... (1)
                                        DevA  ]                 ... (2)
                                        

                        ctype_local = repmat('U', size(bconstraint,1),1);
                        %vartype=repmat('C',NZ_local,1);
                        
                        %% GLPK: min over w's                       
                        %param.msglev=1; % Output GLPK messages     
                        %param.save=1;   % Set save options
                        sense = 1;     % Max (-1) or Min (+1) problem
                        
                        [~,fmin,status,~] = glpk( objective,...
                                                    constraint,...
                                                    bconstraint,...
                                                    lb_local,ub_local,...
                                                    ctype_local,...
                                                    vartype_local,...
                                                    sense,param_local);

                        if status == 2    
                            status_vec(n) = status; 
                            val = lambda*Util_now + fmin; % Add current payoff
                            %warning('ADMIT_OUTER_LPSET:status',...
                            %            '[2] Feasible but not optimal');
                        elseif status == 5
                            status_vec(n) = status;
                            val = lambda*Util_now + fmin;  % Add current payoff
                        else
                            val = -inf;
                            %warning('ADMIT_OUTER_LPSET:status',...
                            %                      'Solution not found');
                        end
                        
                        pival_n_temp(n) = val;

                    end % end for n
                    
                    %GLPK_SUCCESS = 5;
                    %sum_success = sum(status_vec);
                    %all_success = K*GLPK_SUCCESS;
                    
                    %if sum_success == 0
                    %j    fprintf('\nGLPK: 0.00 percent success\n');
                    %end 
                   
                    % if sum_success == all_success
                    %   fprintf('\nGLPK: Success for all samples n\n');
                    %else
                    %    fprintf('\nGLPK: %5.2f percent success n\n',...
                    %                        sum_success*100/(all_success));
                    %end
                    
                    % Max over feasible b's:
                    pival_j_temp(j) = max(pival_n_temp);
                    
                    j = j+1;
                end % End if
            end % End While j

            % min over all j:=: j(a,b): "Choosing minimizer lambda"
            %disp(pival_j_temp)
            
            % Flip -Inf to +Inf's for minizer over a below:
            pival_j_temp(pival_j_temp == -Inf) = Inf;
            pival_a_temp_local(index_a) = min(pival_j_temp); 
                                                    % c+[k,l,(a,b)]
       
    end % EndFor index_a
    
    % Minimize over a:
    pival_k =  min(pival_a_temp_local);
    
    if pival_k == Inf
        pival_k = -Inf;
    end
    
    % Clear garbage and free up parallel workers:
    %matlabpool close

end


