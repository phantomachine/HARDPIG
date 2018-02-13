function val = uvec(self,BxA_k,lambda_k)

% % UVEC.M
% %
% % Batch compute vector of one-period payoffs for agents and government,  
% % induced by all possible (b,a) pairs in feasible action set (B x A){k}, 
% % at current state lambda_k, where k = 1,...,K.
% % 
% % =======================================================================
% %     (c) 2011-- T.Kam and R.Stauber.
% %
% % INPUT: Let N_Z = M + N 
% %
% %     * self       : KSMOD02 class
% %     * .NP        : Number of players: N_Z + 1
% %     * BxA_k      : Slice: feasible action correspondence, (B x A){k}
% %
% % OUTPUT:
% % 
% %     * val        : Numeric array, dimension (rcount*NP x 1) 
% %                    NOTE: Each b_feas{k} is a double array of dimension
% %                    (N_Z x rcount), and rcount varies depending on 
% %                    feasibility at each state k = 1,...,K.                 
% %         
% % Email: mortheus__at__gmail.com
% % =======================================================================
% % $Revision: 4.0.3 $  $Date: 2011/05/02 11:38:20 $
% %
% % See also KSMOD02, UCONS, PHIA

    n = self.N;         % Maximum unemployment state
    m = self.M;         % Maximum employment state
    wage = self.WAGE;   % Constant wage for employed

    N_Z = n + m;        % Number of individual states / small players
    
    % Number of (b,a) vectors in each BxA_k :
    rcount = size(BxA_k,1);

    % All unemployed immediate payoffs at all (b(j),a(j)), j = -N,...,-1 :
    val_u =  ucons( self, BxA_k(:,1:n) ) ...
                        - phia( self, BxA_k(:,N_Z+1:N_Z+n) ) ;

    % All employed immediate payoffs at all (b(j),a(j)), j = +1,...,+M :
    val_e =  ucons( self, wage + BxA_k(:,n+1:N_Z) ) ...
                        - phia( self, BxA_k(:,2*n+m+1:2*(N_Z)) ) ;

    % Let U = ucons(wage,b) + phia(a). So we get the array 
    % val = [ U(-N),...,U(-1),U(+1),...,U(+M) ] :
    val = [ val_u, val_e ] ; 

    % Now add the social planner's immediate payoff at state lambda_k :
    Lambda_now = repmat(lambda_k, rcount, 1); 

    valG = Lambda_now .* val;   % weigh private utilities by lambda_k
    valG = sum(valG,2);         % Utilitarian current payoff

    val = [ val, valG ]';        % Concantenate to last column of val

%     val = reshape(val', rcount*self.NP, 1); % Reshape as (rcount*NP x 1)
%                                             % vector: useful for large GLPK
%                                             % routine in ADMIT algorithm
% 

end