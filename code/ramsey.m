% % RAMSEY.M
% %
% % Calculate:
% % (1) Best dynamic Ramsey commitment policy outcome for every initial
% %     aggregate state (equivalently, game state).
% %     This is MAIN.M without the need to call PUNISHMENT.M
% %
% % =======================================================================
% % (c) 2011-- T.Kam and R.Stauber.
% %
% % DEPENDENCIES:
% %
% %     * GNU Linear Programming Kit (ANSI C to MATLAB MEX)
% %     * CompEcon Toolbox (Miranda-Fackler)
% %         
% % Email: mortheus__at__gmail.com
% % =======================================================================
% % $Revision: 4.0.3 $  $Date: 2011/09/06 00:40:20 $ 
% %
% % See also AGENT, MAIN, GETINDEX, GLPK


clear all
close all
clc

%% ------------------------------------------------------------------------
%% STEP -1: MACHINE SPECIFIC SETTINGS
%% ------------------------------------------------------------------------

DISPLAY_DETAIL_INNER = 0;       % Display option for while/for loops

USE_GLPK = 1;
USE_PARFOR = 1;     POOLSIZE = 8;
USE_PROFILER = 0;

Options_Sol.glpk = USE_GLPK;
Options_Sol.pfor = USE_PARFOR;
Options_Sol.disp = DISPLAY_DETAIL_INNER;

LOAD_OLD = 'yes';               % load W_old1.mat

iter_last = 3;                % Check number in /_mat folder

UNIX = 'on';                    % UNIX users: 'on' ||  MS-Windows: 'off'



    % Add path to directory:
    if strcmp(UNIX,'on') == 1   
        addpath(strcat(pwd,'/glpkmex/'));
        addpath(strcat(pwd,'/_func/'));
        addpath(strcat(pwd,'/eq_sphere_partitions/'));
        addpath(strcat(pwd,'/_mat/'));
    else
        addpath(strcat(pwd,'\glpkmex\'));
        addpath(strcat(pwd,'\_func\'));
        addpath(strcat(pwd,'\eq_sphere_partitions\'));
        addpath(strcat(pwd,'\_mat\'));
    end
    
    if USE_PROFILER == 1
        profile on -history
    end
    
    if USE_PARFOR == 1
        if matlabpool('size') > 0
           matlabpool close force
        end
    end
    
if strcmp(LOAD_OLD,'no') == 1
    
    %% --------------------------------------------------------------------
    %% STEP 0: MODEL PARAMETERIZATION
    %% --------------------------------------------------------------------

    SIGMA = 1;                % CRRA u(.) parameter (linear = 0; log = 1)
    PHI = 0;                  % Disutility of effort parameter (linear = 0)
    M = 1;                    % "Max duration" of employment, M
    N = 2;                    % "Max duration" of unemployment, N
    WAGE = 1;                 % fixed common wage for employed
    REP_RATIO = 0.8;          % "Maximal replacement ratio"
    DELTA = 0.9615;           % common discount factor: Annual rate = 4%

    L = 200;                  % # of vertices per state
    NGRID_P = 8;              % Number of grid elements representing [0,1]

    % NOTES: Dependant parameters set in KSMOD2.M
    % -------------------------------------------
    % MBAR = WAGE*REP_RATIO;    % exogenous upperbound on benefit payout
    % DELSTAR = (1 - DELTA);    % For averaging payoffs

    %% --------------------------------------------------------------------
    %% STEP 1: CREATE INSTANCE OF CLASS: KSMOD02
    %% --------------------------------------------------------------------

    model = ksmod02(SIGMA,PHI,M,N,WAGE,REP_RATIO,DELTA,NGRID_P); 

    fprintf('\nStarting Preliminary calculations ... Please wait.\n');

    %% --------------------------------------------------------------------
    %% STEP 2: PRELIMINARIES
    %% --------------------------------------------------------------------

    N_Z = model.NS;                               % Private state dimension
    N_BGRID = size(model.BGRID,1);                % Grid on policy b(j)
                                                  % where j = 1, ...,N_Z

    np = N_Z + 1;                                 % Payoff vector dimension
                                                  % +1 element: govt payoff   

    % Fix search subgradients; Use Leopardi (2006):      

    H = subgradient(np,L);                      % H is finite subset of the 
                                                % np-dimensional sphere: 
                                                % S^(np-1). So, for each i,
                                                % sum( H(i,:).^2, 2 ) = 1.

    % Compute natural bounds on actions and payoffs:

    [V_MIN, V_MAX, ~] = payoff_bounds(model);

        % Space of private agents' action profiles: 
        [A, NA, N, A_u, A_e] = space_private_actions(model);

        % Transition probabilities from a:
        matlabpool(POOLSIZE);
        P_a = transprob_A( model, NA, A_u, A_e ); % Uses PARFOR
        matlabpool close       

    % State space D: size(D) = (K x N_Z)

    [D, K] = space_distribution(model);

    % Feasible policy correspondence B := {B(k) | k = 1,...,K}
   
    b_set = model.BGRID;       % Set of all possible b
    
    B = feasible(D, b_set);    % A cell array of dimension (K x 1)
                               % NOTE: Each b_feas{k} is a double array of
                               % dimension (nb(k) x N_Z), and nb(k) varies
                               % depending on feasibility at each state k =
                               % 1,...,K.

    % Product space of (B x A), #(B x A){k}, and continuation state space 
    % from (B x A){k}:

    [BxA, row_BxA, next_State] = transit( A, B, P_a, D );


    %% --------------------------------------------------------------------
    %% STEP 3: OUTER APPROX. OF SPE PAYOFF CORRESPONDENCE
    %% --------------------------------------------------------------------

    minU = V_MIN*ones(model.NP, 1); % Lower bound on payoff vector
    maxU = V_MAX*ones(model.NP, 1); % Upper bound on payoff vector

    % GLPK settings:
%     ctype_punish = repmat( 'U', L, 1 );     
%                                     % Variable with upper bound - i.e. LP 
%                                     % constraints of form: A(i,:)*x <= b(i)
% 
%     vartype_punish = repmat( 'C', np, 1 );   % Continuous type variables
    
    % Initial guess of hyperplane levels (payoff correspondence W_outer)
    C_old = max(max(H*repmat(V_MAX,np,1)),max(H*repmat(V_MIN,np,1)))...
                                                                *ones(L,K);
	iter_last = 0;				   
	iter = 1;

elseif strcmp(LOAD_OLD,'yes') == 1
    % Initial guess of hyperplane levels (payoff correspondence W_outer)
    %load sse_005					% W_old
	iter_last_new = iter_last;
    load(strcat('_mat/ramsey_iter',int2str(iter_last),'.mat'));
	iter = iter_last_new + 1;
end

C_new = zeros(size(C_old));           % Update array for hyperplane levels
eps_vec = zeros(1, model.MAXITER1);   % Record distances

% Initialize WHILE loop:
fprintf('\nStarting Outer Monotone Operator iteration ... Please wait.\n');


epsi = 1;

tic
while epsi >= model.TOL1 && iter <= ( iter_last + model.MAXITER1 )
    
    if DISPLAY_DETAIL_INNER ~= 1
        fprintf('Currently at iter = %i. ... Please Wait\n', iter);
    end
    
%     % STEP 3.1: Construct current guess of consistent punishment payoff:
%     minw_G = punish(H, C_old, minU, maxU, ctype_punish, vartype_punish);
   
    % STEP 3.2: ADMIT (admissibility algorithm)
    C_new = admit_ramsey(model, C_old, D, K, L, H, BxA, row_BxA, ...
                          minU, maxU, iter, Options_Sol, POOLSIZE);

    % STEP 3.3: CHECK AND STORE DISTANCE
    epsi = mdistance(C_old,C_new);  % d(Wo(iter),Wo(iter+1))  
    eps_vec(iter) = epsi;

	fprintf('%i\t d(Wo_n,Wo_n+1) = %0.5g\n', iter, epsi);
    
    % Intermediate saves for future use:
    if epsi < 1
        filename = strcat('_mat/ramsey_iter',num2str(iter),'.mat');
        save(filename)
    end
        
    % STEP 3.4: UPDATE CURRENT GUESS OF LEVELS
    C_old = C_new;                    
    
    iter = iter + 1;    
end
TotalElapsedTime = toc;
AverateTimePerOperator = TotalElapsedTime/(iter - 1 - iter_last);
save('_mat/ramsey_final')


if USE_PROFILER == 1
    profile viewer
    p = profile('info');
    profsave(p,'profile_results')
end