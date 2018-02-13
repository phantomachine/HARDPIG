function k_plus = transit(self, a, k, D)

% % TRANSIT.M
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
% %     * self  : instance of model class KSMOD2
% %     * a     : (N_Z x 1) vector of actions
% %     * k     : Current state
% %     * D     : Finite state space, (K x N_Z)
% %     
% % OUTPUT:
% %
% %     * index : index = 1,...,K 
% %
% % DEPENDENCIES:
% %
% %     * CompEcon Toolbox (Miranda-Fackler)
% %
% % Email: mortheus__at__gmail.com
% % =======================================================================
% % $Revision: 4.0.3 $  $Date: 2011/05/04 11:38:20 $

au = a( 1 : self.N );               % Decision of -j's
ae = a( self.N + 1 : end );         % Decision of +j's
 
P = self.markovmat(self,au,ae);     % Markov matrix

lambda = D(k,:);                    % Current state vector (dist of agents)

lambda_plus = lambda * P;           % Continuation state vector
                                    % (part of consistency requirement)
                                    % But lambda_plus not in D w.p.1, so...

k_plus = getindex(lambda_plus, D);  % Map back to natural number index, so:
                                    % k_plus \in {1,...,K}
                                    
end

