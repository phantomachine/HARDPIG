function C_new = admit(self, C_old, D, K, L, H, BxA, row_BxA, ...
                           minU, maxU, minw_G, iter, Options_Sol, POOLSIZE)
                            

% % ADMIT.M
% %
% % STEP 1: Consistency (part)
% % --------------------------
% % Function computes approximate index of continuation (natural) state 
% % vector as follows. Given fixed action profile of players a (N_Z x 1), 
% % and given current state lambda :=: k, k = 1,...,K, return:
% %
% %         lambda(+) = lambda P(a).
% % 
% % where P is Markov transition function, controlled by a, and K is the
% % cardinality of D, a finite set representing the space of probability 
% % measures lambda.
% %
% % STEP 2: Approximate consistency
% % -------------------------------
% % Find approximate index k(+) :=: lambda(+) in the finite set D of 
% % lambda's. This is done using the C-MEX function GETINDEX.
% %
% % =======================================================================
% %     (c) 2011-- T.Kam and R.Stauber.
% %
% % INPUT: 
% %     * self  : instance of self class KSMOD2
% %     * C_old : Previous guess or output from ADMIT, (L x NP)
% %     * D     : Finite state space, (K x N_Z)
% %     * K     : # vectors of D
% %     * L     : # subgradients, hyperplane normals
% %     * H     : Set of subgradients, (L x N_Z)
% %     * BxA   : Feasible policy correspondence, cell (K x 1)
% %     * row_BxA : #rows for each array element of BxA(k), k = 1,...,K
% %     * minU  : Lower bound on payoffs (NP x 1)
% %     * maxU  : Upper bound on payoffs (NP x 1)
% %     * minw_G : Government's min-payoff function, (1 x K)
% %     * iter  : Current iteration of ADMIT
% %     * Options_Sol : Structure containing options (see below)
% %     * POOLSIZE : # MATLAB workers to open
% %     
% % OUTPUT:
% %
% %     * C_new : Update on C_old respecting ADMISSIBILITY
% %
% % DEPENDENCIES:
% %
% %     * CompEcon Toolbox (Miranda-Fackler)
% %
% % Email: mortheus__at__gmail.com
% % =======================================================================
% % $Revision: 4.0.3 $  $Date: 2011/05/04 11:38:20 $
% %
% % See also GLPK, UVEC, MATLABPOOL, SETDIFF, LINPROG, GETINDEX

    USE_GLPK = Options_Sol.glpk;
    
    DISPLAY_DETAIL_INNER = Options_Sol.disp;

    N_Z = self.NS;        % # private states
    DELTA = self.DELTA;   % Discount factor
    np = self.NP;         % # players
    
    % GLPK settings:
    param.msglev = 1;     % Output only GLPK messages if Errors
    sense = -1;            % Minimization = 1; Maximization = -1
    
    for k = 1 : K                               % states: lambda(k) 
    
        BxA_k = BxA{k};
        rcount = row_BxA(k);
        
        lambda_k = D(k,:);
        
        U_BxA_k = (1-DELTA) * uvec( self, BxA_k, lambda_k );
        
        %fprintf('iter = %i\t k = %i\n', iter, k);
              
        %c_pp = zeros(rcount,1);  % inside j loop
        matlabpool(POOLSIZE)
        
        parfor l = 1 : L                           % subgradients: h(l)
        
            if DISPLAY_DETAIL_INNER == 1
                fprintf('iter = %i\t k = %i\t l = %i\n', iter, k, l);  
            end
				   
            % Copy H array again to parallel workers:
            Hl = H;
            
            h_l = Hl(l,:);
            
            % Copy minw_G array again to parallel workers:
            minwG = minw_G;
            
            % Copy C_old and c_pp array again to parallel workers:
            cold = C_old;
            cplus = zeros(rcount,1);
            
                % Set (payoff) optimizer lower and upper bounds for GLPK:
                lb = minU;                      % Minimal lifetime payoffs       
                ub = maxU;                      % Maximal lifetime payoffs

                % GLPK constraint type:
                ctype = repmat( 'U', L + 1, 1 );
                
                % GLPK variable type:
                vartype = repmat( 'C', np, 1 );
				   
			    % GLPK settings:
				%param.msglev = 1;     % Output only GLPK messages if Errors
				%sense = -1;            % Minimization = 1; Maximization = -1
            
            U_BxA_kj = U_BxA_k;
            BxA_kj = BxA_k;
                
            for j = 1 : rcount
                                   
                % Immediate payoffs at each fixed (b,a) - not needed in LP.
                % But this is needed when updating the current vector of
                % total payoffs:
                

                U_now = U_BxA_kj(:,j) ;

                current_payoff = h_l * U_now; % (rcount x 1)

                % Now get continuation payoffs in ADMIT:
                continue_payoff = DELTA*h_l; % Objective function in LP
                                                     % (rcount x rcount*NP)

                % Constraints: Consistency and Incentive Compatibility

                IC =  [ zeros(1,N_Z), -DELTA ];

                Confun = [ Hl ;        % Consistency: H*w <= C_old(:,k')
                           IC ];

                % Set of current guess of hyperplane levels:  

                kp = BxA_kj(j, end); % Continuation state from j := (b,a)
                Cnow = cold(:,kp);

                % maxmin payoff:
                b_eq = BxA_k( j, 1:N_Z ); % current policy vector, b

                [ ~, not_b_idx ] = setdiff( BxA_k(:,1:N_Z), ...
                                                             b_eq, 'rows');
                                                            % deviation, b'
                kp_deviate = BxA_k(not_b_idx, end); 
                                             % Set of possible continuation 
                                             % states from deviation play
                                             % (b',a'(b'))

                maxmin_vG = max( minwG( kp_deviate ) ) ;
                                               % maxmin payoff over all b'

                U_g = U_now(np) - maxmin_vG;
                                           % upper bound on govt. incentive 
                                           % compatible continuation value 

                bcon = [ Cnow   ;
                         U_g   ];


                
                % Choose w to min { continue_payoff*w | Confun*w <= bcon }:
                
                
                if USE_GLPK == 1
                    [~,fmax,~,~] = glpk( continue_payoff, ...
                                                        Confun, bcon,...
                                                        lb, ub, ...
                                                        ctype, vartype, ...
                                                           sense, param ) ;
                    
                else
                    % Use MATLAB LINPROG
                    options_lp = optimset(  'Display', 'off', ...
                                            'LargeScale', 'on', ...
                                            'Simplex', 'off',...
                                            'UseParallel','never');
                    
                    [~,fmax,~,~] = linprog(-continue_payoff,...
                                                           Confun, bcon,...
                                                                [],[],...
                                                                lb, ub,...
                                                            [],options_lp);
                    fmax = -fmax;
                end
                
                %fprintf('status = %i\n',status);
                
                % Add back h_l weighted current payoff
                cplus(j) = current_payoff + fmax;
            
				% clear xopt fmax status extra
                % clear functions
            end % j
            
            % STEP 3.3: Update operator iteration:

            C_new(l,k) = max(cplus);   
                        
        end % l
        
        matlabpool close
        
    end % k