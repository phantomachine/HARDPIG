%% STRAT.M
%
% Construct sample pure strategies supporting initial (lambda0, w0) given
% inner approximating correspondence W^i := [ G(l,k), c(l,k) ].
%
% (c) 2012, 2013, Timothy Kam. Email: tcy.kam__at__gmail.com 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Use subject to GNU LGPL licensing terms. Cite this header and 
% author completely in subsequent re-use and modifications.
% =========================================================================
% $Revision: 5.0.0 $  $Date: 2013/09/11 12:45:20 $

close all
clear all
clc

%% SETTINGS

T = 1;

Nagrid = 5;
Nbgrid = 6;

% Small tolerance for GBC deficit/surplus
tolB = 0.15;

% GLPK options
    param.msglev=0; % Output GLPK messages on workspace      
    param.save=0;   % Set save options

% load Fixed policy result
load('_mat/ramsey_lr01.mat')

lambda_now = lambda_ramsey_ss_max;

wrss_now = vrmax;


% Load saved results from SSE
load('_mat/sse2014_iter162.mat')

    L = size(lambda_new, 3);
    K = size(lambda_new, 4);

    % lambda_new and w_new are (1 x 3 x L x K ) arrays
    vG_k = sum( lambda_new .* w_new, 2 );
    vG_k = reshape(vG_k, L, K);
    [vG_k_max, max_index] = max(vG_k,[],1);
    [vG_k_min, min_index] = min(vG_k,[],1);

    wmin = zeros(K,N_Z);
    wmax = wmin;
    lmin = wmax;
    lmax = lmin;

    for k = 1 : K
        wmin(k,:) = w_new(:,:,min_index(k),k)'; % w supporting min(vG)
        wmax(k,:) = w_new(:,:,max_index(k),k)'; % w supporting max(vG)
        lmin(k,:) = lambda_new(:,:,min_index(k),k)';
                        % Not needed: lamda supporting min(vG)
        lmax(k,:) = lambda_new(:,:,max_index(k),k)';
                        % Not needed: lamda supporting max(vG)
    end
    
%% ACTION PROFILE SETS

    % Agent Action set
    a_j = linspace(0.01, amax(model), Nagrid)'; % j's action set
    %A_j = a_j;

    agrid = a_j;
    for i = 1 : N_Z-1
        agrid = gridmake(agrid,a_j);
    end

    ProfileSetA = agrid;
    NA = size(ProfileSetA,1);


    % Government action set
    b_u = linspace(0.01,model.MBAR,Nbgrid)';
    b_e = linspace(-model.WAGE*0.5, min(-0.1,model.MBAR),...
                                                     Nbgrid)';

    %b_u = linspace(0.0,self.MBAR,self.Nbgrid)';
    %b_e = linspace(-self.WAGE*0.25, min(-0.1,self.MBAR),self.Nbgrid)';
    bsu = b_u;
    for i = 1:N-1
        bsu = gridmake(bsu,b_u);
    end

    bse = b_e;

    for i = 1:M-1
        bse = gridmake(bse,b_e);
    end

    ProfileSetB = gridmake(bsu,bse); % requires GRIDMAKE: COMPECON
    NB = size(ProfileSetB,1);

    model.ProfileSetB = ProfileSetB;
    model.ProfileSetA = ProfileSetA;
    
    parallel = 'off';
    P_a = TransProbA(model,parallel);
    
    delstar = (1-DELTA);
    

%% CONSTRUCT SAMPLE SSE PATH (INDEXES TO (a,b) in A x B
alpha = cell(T,1);
beta = cell(T,1);

% Storage
    val_a = zeros(1,NA);
    wnext_a = cell(NA,1);
    
% Initial conditions
    k_now = InPartElement(lambda_now, dt, DtriVerts);
    
    % Best SSE payoff vector over Q(k_now):
    w_now = wmax(k_now,:)';
    
    v_now = vG_k_max(k_now);
    
for t = 1 : T
   
    %% Work out (approximately) feasible set of policies b, B(lambda)
    
    Feas_temp = lambda_now * ProfileSetB';
    B_temp = ProfileSetB( (abs(Feas_temp) <= tolB), :);
    NB = size(B_temp,1);
        
    % Inner loops:
    for b_idx = 1 : NB
        
        % Variable storage:
        val_b = zeros(NB,1);
        aoptindex_b = cell(NB,1);
        wnext_b = cell(NB,1);
        
        parfor a_idx = 1 : NA 
            
            lb_local = lb;
            ub_local = ub;
            param_local = param;
            vartype_local = vartype;
            
            ProfileSetA_local = ProfileSetA;
            P_a_local = P_a;
            B_temp_local = B_temp;
            
            C_new_local = C_new;
            
            delstar_local = delstar;
            
            % Current payoff from a: effort disutility
            Util_a = -phia(model,ProfileSetA(a_idx,:)); 
            Util_ia = repmat(Util_a, NA-1, 1);
            
            % All other a# \neq a, and induced P(a#) matrices:
            Index_A = 1:NA;                     % indexes to { a }
            ia_Not = setdiff(Index_A, a_idx); 
                                                % All a# \neq a
            P_ia_Not = P_a_local(:,:,ia_Not);   % All P(a#) \neq P(a), 
                                                % N x N x card(A)-1
            P_ia_not = reshape(P_ia_Not,...
                                N_Z*(NA-1),N_Z);
                                            % Tile this too!
                                            
            % Current payoffs from all deviations from a:        
            Util_ia_Not = -phia(model, ProfileSetA_local(ia_Not,:) );
            
            % Transitions
            P_now = P_a_local(:,:,a_idx);
            P_ia = repmat(P_now,NA-1,1);
            
            % Agents' incentive constraints at current {k,l,a}
            DevP = DELTA*(P_ia_not - P_ia);
            DevA = delstar_local*(Util_ia - Util_ia_Not);
            [nr, nc ] = size(DevA);
            DevA = reshape(DevA', nr*nc, 1);
            
            % Current payoff from b:
            b_current_feasible = B_temp_local(b_idx,:);
            cons = b_current_feasible;
            cons(N+1:end) = cons(N+1:end) + model.WAGE;
            Util_b = ucons(model, cons);  

            % Total current payoff:        
            Util_now = delstar_local*( Util_a + Util_b )'; 
                                            % column vector

                                                        
            % Get continuation game state:
            
            lambda_next = lambda_now * P_now;
            
            % Check which Q(k) contains lambda_next:        
            k_next = InPartElement(lambda_next, dt, DtriVerts);
            
            % Get appropriate partitition element's correspondence:
            c_next = C_new_local(:,k_next);
            
            % Solve LP problem for maximizers w_next that enforces current
            % action-policy pair (a,b):
            
            constraint = [  H;                      %... (1)
                            DevP;                   %... (2)
                            lambda_now * DELTA*P_now   ];        %... (3)
                                        
            bconstraint1 = [ c_next;                              %...  (1)
                            DevA                      ];          %...  (2)
            bconstraint2 =  v_now - lambda_now * Util_now ;    %...  (3)
                                          
                bconstraint = [ bconstraint1; bconstraint2 ] ;
                
                
            objective = lambda_next; %ones(1,N_Z);
            
                    % Constratint types
                    ctype_local1 = repmat('U', size(bconstraint1,1),1);
                                                        % Inequality
                                                        % constraints
                    ctype_local2 = repmat('S', size(bconstraint2,1),1);
                                                        % Equality
                                                        % constraints
                    
                    ctype_local = [ ctype_local1; ctype_local2 ];
                    
                    %vartype = repmat('C',N_Z,1);

                    %% GLPK                       
                    %param.msglev=1; % Output GLPK messages     
                    %param.save=1;   % Set save options
                    sense = -1;     % Max (-1) or Min (+1) problem

                    [w_temp,fval,status,~] = glpk( objective,...
                                                constraint,...
                                                bconstraint,...
                                                lb_local,ub_local,...
                                                ctype_local,...
                                                vartype_local,...
                                                sense,param_local);
                    
                    wnext_a{a_idx} = w_temp;
                    
                    if status == 2    
                        %status_vec(n) = status; 
                        %warning('STRAT:status',...
                        %            '[2] Feasible but not optimal');
                        val_a(a_idx) = fval;
                    elseif status == 5
                        %status_vec(n) = status;
                        val_a(a_idx) = fval;
                    else
                        val_a(a_idx) = -inf;
                        %warning('STRAT:status',...
                        %                      'Solution not found');
                    end
                    
        
            
            
        end % EndFor a
        
        
        
        % Given b, maximize over the a's:
        
        [ vmax_temp, ind_temp  ] = max(val_a);
        
        val_b(b_idx) = vmax_temp; 
        aoptindex_b{b_idx} = ind_temp;
        
        wnext_ab_temp = wnext_a{ind_temp};
        wnext_b{b_idx} = wnext_ab_temp;
        
        fprintf('t = %i \t b_idx = %i \t \n',...
                                                t, b_idx);
        disp(wnext_ab_temp);
        
    end % EndFor b
    
    % Maximize over the feasible b's:
    
    [~, boptindex ] = max(val_b);
    w_next = wnext_b{boptindex};
        
    % Store indexes of (a*, b*) at stage/date t:    
    alpha{t} = aoptindex_b{boptindex};
    beta{t} = boptindex;
    
    % Update state and continuation payoffs
    lambda_now = lambda_now * P_a(:,:,alpha{t});
    w_now = w_next;
    
end % EndFor t
% 
% % Pack it all in:
% strategy.a = alpha;
% strategy.b = beta;