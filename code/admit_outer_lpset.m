function [Cnew, lambda_new, w_new ] = admit_outer_lpset(self,...
                                        PolyTriIndex,Fset,H,C,P_a,pival,...
                                                  lb,ub,vartype,epsi, iter)

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
% See also PUNISHK (or PUNISHK_LPSET optional LP version)

%% Preliminaries:

    % GLPK options
    param.msglev=0; % Output GLPK messages on workspace      
    param.save=0;   % Set save options

    [L, K] = size(C);   % No. of subgradients (L) and partition elements (K)
    NZ = size(H,2);     % Number of agent classes
    maxN = self.N;      % Max. number of -j agent states
    
    % Extract set of private agent action profiles A:
    A = self.ProfileSetA; 
    NA = self.NA;
    
    % Extract set of government action profiles B:
    B = self.ProfileSetB;
    
    % Number of action profiles (a,b) \in A x B:
    %Num_ab_Profiles = size(IndexProfileSet,1); 
    
    % Wage parameter in utility function for agents type j > 0:
    wage = self.WAGE;
    
    % Discount factor:
    delta = self.DELTA;
    delstar = self.DELSTAR; % 1 - delta

    % Preallocate memory:
    Cnew = C;
    lambda_new = zeros(1,NZ,L,K);
    w_new = zeros(1,NZ,L,K);
    
    cplus_a_temp = zeros(NA,1);
    
%% Main Loop: Admissibility w.r.t. (C)
    
fprintf('\tADMIT_OUTER_LPSET.M: ');

for k = 1 : K               % Loop: partition elements Q(k)

    matlabpool 
    fprintf('\t iter = %i | k = %i | epsi(iter-1) = %6.3g\n',...
                                                       iter, k, epsi);

    parfor l = 1 : L        % loop: subgradients h(l)

%       fprintf('\t\t\t< %i >\n',l);

        % Local worker copies of data:

        param_local = param;
        %ctype_local = ctype;
        vartype_local = vartype;
        lb_local = lb;
        ub_local = ub;

        N = maxN;       % Max. number of -j states
        NZ_local = NZ;   % Number private agents states, card(Z)
        %NA_local = NA;  % Number of private action profiles a

        delta_local = delta;        % discount factor
        delstar_local = delstar;    % 1 - delta
        wage_local = wage;          % wage level

        C_local = C;
        H_local = H;
        h = H_local(l,:);   % Current subgradient
        A_local = A;        % Set of private action profiles a
        B_local = B;        % Set of government action profiles b

        % Feasible sets F{a}{k}{i} where i := i(a,k)
        F = Fset;           % Feasible (lamda,b) sets: GBC >= 0

        % N x N x card(A) set of Markov matrices                               
        P_a_local = P_a;    
        %PolyLcons_local = xPolyLcons;
        PolyTriIndex_local = PolyTriIndex;

        % Natural bounds on total average discounted payoffs:
        % Vmin = V_MIN; % See lb
        % Vmax = V_MAX; % See ub

        % (K x 1) vector of punishment values from PUNISH_OUTER.M:
        pival_local = pival;

        % Preallocate:
        cplus_a_temp_local = cplus_a_temp;
        lambda_a_temp = zeros(1,NZ_local,NA);
        w_a_temp = lambda_a_temp;

        for index_a = 1 : NA % Loop: profiles a

            % Current P(a) matrix under a:
            P = P_a_local(:,:,index_a);
            %P_ia = repmat(P,NA_local-1,1); % Tile it NA-1 times!

            % Current action profile a:
            a = A_local(index_a,:);         % row vector
            %b = B_local(ib,:);             % row vector

            % All other a# \neq a, and induced P(a#) matrices:
            %Index_A = 1:NA;                     % indexes to { a }
            %ia_Not = setdiff(Index_A, index_a); 
                                                % All a# \neq a
            %P_ia_Not = P_a_local(:,:,ia_Not);   % All P(a#) \neq P(a), 
                                                % N x N x card(A)-1
            %P_ia_not = reshape(P_ia_Not,...
            %                    NZ_local*(NA_local-1),NZ_local);
                                            % Tile this too!
            
            %% Agent's best response:
            %
            % Current payoff from a: effort disutility
            Util_a = -phia(self,a); 
            %Util_ia = repmat(Util_a, NA_local-1, 1);
            
            % Current payoffs from all deviations from a:
            %Util_ia_Not = -phia(self, A_local(ia_Not,:) );
            
            % Agents' incentive constraints at current {k,l,a}
            %DevP = delta_local*(P_ia_not - P_ia);
            %DevA = delstar_local*(Util_ia - Util_ia_Not);
            %[nr, nc ] = size(DevA);
            %DevA = reshape(DevA', nr*nc, 1);
           

            %% 4. Solve feasible LP on each intersecting sub-state-space

            % Check # partitions intersect with current P(a)(Q_k):
            J = numel(F{index_a}{k}); % Each given P(a)[Q_k] intersects with 
                                 % J other partition elements Q_j

            % Storage for weighted-payoff level, distro and payoff vector:
            cplus_j_temp = zeros(J,1);
            lambda_j_temp = zeros(1,NZ_local,J);
            w_j_temp = zeros(1,NZ_local,J);
            
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
                    cplus_j_temp(j) = -inf;
                    j = j + 1;
                    
                else % Feasible GBC set at j

                    % Current number of feasible pairs (lambda,b) at
                    % {k,l,a}
                    N_feasible = numel(index_lambda);
                    
                    % Catch GLPK errors in current j
                    status_vec = zeros(N_feasible,1);
                    
                    cplus_n_temp = status_vec;
                    lambda_n_temp = zeros(1,NZ_local,N_feasible);
                    w_n_temp = zeros(1,NZ_local,N_feasible);
                    
                    for n = 1:N_feasible
                        
                        %% Set up current payoff U(a,b)
                        % Current state lambda:
                        lambda = set_lambda( index_lambda(n), :); 
                        
                        % Storage for final use:
                        lambda_n_temp(:,:,n) = lambda;

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
                        pival_k = pival_local(k);
                        
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
                        objective = (h * conval);
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
                        % (3) Admissibility: Govt incentive compatibility

%                         constraint = [  H;                      ... (1)
%                                         DevP;                   ... (2)
%                                        -lambda*conval  ];       ... (3)
%                                         
%                         bconstraint = [ C_kp;                   ... (1)
%                                         DevA;                   ... (2)
%                                         lambda*Util_now - pival_k ]; ... (3)
                                        
                        constraint = [  H;                      ... (1)
                                       -lambda*conval  ];       ... (3)
                                        
                        bconstraint = [ C_kp;                   ... (1)
                                        lambda*Util_now - pival_k ];... (3)

                        ctype_local = repmat('U', size(bconstraint,1),1);
                        %vartype=repmat('C',NZ_local,1);
                        
                        %% GLPK                       
                        %param.msglev=1; % Output GLPK messages     
                        %param.save=1;   % Set save options
                        sense = -1;     % Max (-1) or Min (+1) problem
                        
                        [wmin,fmin,status,~] = glpk( objective,...
                                                    constraint,...
                                                    bconstraint,...
                                                    lb_local,ub_local,...
                                                    ctype_local,...
                                                    vartype_local,...
                                                    sense,param_local);

                        if status == 2    
                            status_vec(n) = status; 
                            val = h*Util_now + fmin; % Add current payoff
                            %warning('ADMIT_OUTER_LPSET:status',...
                            %            '[2] Feasible but not optimal');
                            w_n_temp(:,:,n) = wmin;
                        elseif status == 5
                            status_vec(n) = status;
                            val = h*Util_now + fmin;  % Add current payoff
                            w_n_temp(:,:,n) = wmin;
                        else
                            val = -inf;
                            %warning('ADMIT_OUTER_LPSET:status',...
                            %                      'Solution not found');
                            w_n_temp(:,:,n) = wmin*NaN;
                        end
                        
                        cplus_n_temp(n) = val;

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
                                    
                    [ cplus_j_temp(j), index ] = max(cplus_n_temp);
                    lambda_j_temp(:,:,j) = lambda_n_temp(:,:,index);
                    w_j_temp(:,:,j) = w_n_temp(:,:,index);
                    
                    j = j+1;
                end % End if
            end % End While j

            % max over all j:=: j(a,b):
            [ c_j_temp, index_j ] = max(cplus_j_temp);
            cplus_a_temp_local(index_a) = c_j_temp; 
                                                   % c+[k,l,(a,b)]
            lambda_a_temp(:,:,index_a) = lambda_j_temp(:,:,index_j);
            w_a_temp(:,:,index_a) = w_j_temp(:,:,index_j);
            
        end % EndFor a

        [ cp_temp, index_amax ] =  max(cplus_a_temp_local);

        Cnew(l,k) = cp_temp;    % Levels at partition Q(k), direction h(l)
        lambda_new(:,:,l,k) = lambda_a_temp(:,:,index_amax);
        w_new(:,:,l,k) = w_a_temp(:,:,index_amax);
        
        
    end % EndFor l

    % Clear garbage and free up parallel workers:
    matlabpool close

end % EndFor k

fprintf('\t\tADMIT_OUTER.M: Done!\n');

end
