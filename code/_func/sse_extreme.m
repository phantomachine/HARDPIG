function vG_sse = sse_extreme(C_old, np, V_MAX, V_MIN)

% SSE_EXTREME.M
% Compute boundaries (min and max) total expected payoffs given (C_old, H)
% from running SSE routine in MAIN.M
%
% T.Kam, 2011

% Best (social) payoff:

L = size(C_old,1);
K = size(C_old,2);
NP = np;

social_payoff = zeros(NP,1);
social_payoff(end) = 1;

H = subgradient(NP,L);
Confun = H;

% GLPK settings:
param.msglev = 1;     % Output only GLPK messages if Errors

ctype = repmat( 'U', L, 1 );
vartype = repmat( 'C', NP, 1 );

vG_max = zeros(K,1);
vG_min = vG_max;

lb = V_MIN*ones(NP, 1); % Lower bound on payoff vector
ub = V_MAX*ones(NP, 1); % Upper bound on payoff vector

for k = 1 : K
    
    % Worst social payoff at state k
    
    [~,vG_min(k),~,~] = glpk( social_payoff, Confun, C_old(:,k), lb, ub,...
                                                        ctype, vartype, ...
                                                           1, param ) ;
    % Best social payoff at state k
    [~,vG_max(k),~,~] = glpk( social_payoff, Confun, C_old(:,k), lb, ub,...
                                                        ctype, vartype, ...
                                                           -1, param ) ;                                                    
end

vG_sse = [vG_min, vG_max];