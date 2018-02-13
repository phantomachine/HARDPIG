function [BxA,row_BxA,next_State] = transit(A,B,P_a,D)

% % TRANSIT.M
% %
% % Let k = 1,...,K be index of current states in D -- e.g. D(k,:) is a row
% % vector of probability distribution.
% %
% % Function computes product policy-action space (B x A){k} at each k.
% % Given a pair (b,a) in B x A, function also computes the number of
% % policy action vectors, row_BxA{k}, at each k. 
% %
% % Finally, function computes the continuation state in D; equivalently,
% % finds k' in {1,...,K} arising from each (b,a) at each k.
% %
% % STEP 0:
% % --------------------------
% % Construct set (B x A){k}.
% %
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
% % INPUT: Let N_Z = M + N 
% %
% %     * A          : Combinations for all a's; Dim: NA x N_Z
% %     * B          : Combinations for all b's; Dim: ? x N_Z
% %     * P_a        : Output from transprob_A; Dim: N_Z x N_Z x NA
% %     * D          : Array of possible state vectors, K x N_Z
% % 
% % OUTPUT
% % 
% %     * BxA        : Product policy-action space (last column is integer 
% %                    index, mapping one to one to last column of P_a)
% %                    Index points to action profile a = 1,...,NA
% %     * row_BxA    : Size for all (b,a) pairs at each k
% %     * next_State : Continuation state space from set BxA{k} defined by
% %                    {k' = 1,...,K | (b,a){k} in (B x A){k}, k = 1,...,K}
% %         
% % Email: mortheus__at__gmail.com (T.Kam)
% % =======================================================================
% % $Revision: 4.0.3 $  $Date: 2011/05/02 11:38:20 $  
% %
% % See also GRIDMAKE, TRANSPROB_A

    [ K, N_Z ] = size(D);                        % # state vectors
    
    Dset = D;
    
    BxA = cell(K,1);
    row_BxA = zeros(K,1);
    next_State = cell(K,1);

    parfor k = 1 : K

        % STEP 0: feasible (b,a) at k:
        BxA{k} = gridmake( B{k}, A );     
        row_BxA(k) = size( BxA{k}, 1 );   % Measure number of vectors (b,a) 
                                          % in each discretized (B x A){k}

        PA = P_a;                         % Splice P_a data for PARFOR

        lambda_k = D(k,:);                % Get current state
        A_comb_k = BxA{k}(:,N_Z+1:end);   % Get current (b,a){k}


        % STEP 1: Transition of states:
        lambda_next = zeros(row_BxA(k), N_Z); 

        %     A_u = A_comb_k( :, 1:N );
        %     A_e = A_comb_k( :, N+1:N_Z );

        for j = 1 : row_BxA(k)   
            a_idx = A_comb_k(j,end);
            P_a_j = PA(:,:,a_idx);
            lambda_next(j,:) = lambda_k*P_a_j; 
        end

        % STEP 2: Convert to k' indexing for all continuation states:
        next_State{k} = getindex(lambda_next,Dset);
    end
end