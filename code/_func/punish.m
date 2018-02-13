function minw_G = punish(H, C, minU, maxU, ctype, vartype)

% % PUNISH.M
% %
% % Algorithm computes punishment values for the government at each state,
% % given a guess of the multifunction characterizing continuation values
% % of all players, embedded in the array describing hyperplane level sets,
% % c(l,k).
% % 
% % Solution is a linear program using GLPK-MEX (simplex method). Returns
% % -Inf if no solution exists.
% %
% % =======================================================================
% %     (c) 2011-- T.Kam and R.Stauber.
% %
% % INPUT:
% %
% %     * H         : Set of fixed search subgradients
% %     * C         : Array containing c(l,k), l = 1,...,L; k = 1,...,K
% %     * minU      : Lower bound: lifetime payoff vector 
% %     * maxU      : Upper bound: lifetime payoff vector
% %     * ctype   	: repmat( 'U', L, 1 ); Constraint type: H*w <= C(:,k)                                                        
% %     * vartype   : repmat( 'C', L, 1 ); Continuous variable type: w
% %
% % OUTPUT:
% %
% %     * minw_G : Optimal value (worst payoff set), minw_G(k), k = 1,...,K
% %     * wmin   : Optimizer(s) -- could be multiple
% %
% % Email: mortheus__at__gmail.com (T. Kam)
% % =======================================================================
% % $Revision: 4.0.3 $  $Date: 2011/05/02 11:38:20 $
% %
% % See also GLPK, MAIN

    K = size(C,2);                  % # of states, k = 1,...,K
    %L = size(C,1);                 % # of vertices, l = 1,...,L
    N = size(H,2);                  % # players, n = 1,...,N
    
    minw_G = zeros(1,K);
    statusK = zeros(1,K);
    
    param.msglev = 1;               % Output only GLPK messages if Errors
    sense = 1;                      % Minimization = 1; Maximization = -1
    
    objectif = zeros(1,N);          % Objective function in LP problem
    objectif(N) = 1;
   
    lb = minU;                      % Minimal lifetime payoffs       
    ub = maxU;                      % Maximal lifetime payoffs
    
    parfor k = 1 : K
        
        % Solve linear program at each k:     
        Cnow = C(:,k);              % Set of current hyperplane levels
        
        % Linear program - Choose w to min { objectif*w | H*w <= C(:,k) }:
        [~,fmin,status,~] = glpk( objectif, H, Cnow, ...
                                        lb, ub, ctype, vartype, ...
                                                        sense, param ) ;        
        % GLPK exit status check:                                            
        if status == 2              % 2 = feasible or 
            minw_G(k) = fmin;
            warning('PUNISH:status', ...
                        'PUNISH: GLPK-2 > Solution Feasible not Optimal');
        elseif status == 5          % 5 = optimal
            minw_G(k) = fmin;
            statusK(k) = status;
        else
            minw_G(k) = -Inf;
            warning('PUNISH:status', 'PUNISH: GLPK > Solution Not Found');
        end
               
    end
    
    GLPK_SUCCESS = 5;
    
    if sum(statusK) == K*GLPK_SUCCESS
        sprintf('\nPUNISH/GLPK: Success for all (K = %i) states!\n',K)
    end
    
end