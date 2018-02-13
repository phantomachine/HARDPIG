function B = feasible(D, polComb)

% % FEASIBLE.M
% %
% % Function computes discrete multifunction describing (finite) sets B(k) 
% % of feasible policy vectors at each state indexed by k = 1,...,K.
% % 
% % NOTE: This version assumes a balanced per-period budget policy rule.
% % So a feasible policy vector b(k), where k --> lambda(k) is
% % bijective, must satisfy:
% %
% %       lambda(k)*b(k) >= 0.
% %
% % =======================================================================
% %     (c) 2011-- T.Kam and R.Stauber.
% %
% % INPUT: Let N_Z = M + N 
% %
% %     * D          : Array of possible state vectors, K x N_Z
% %     * polComb    : Combinations for all b's: Use GRIDMAKE
% % 
% % OUTPUT
% % 
% %     * B          : A cell array of dimension (K x 1) 
% %                    NOTE: Each b_feas{k} is a double array of dimension
% %                    (nb(k) x N_Z), and nb(k) varies depending on 
% %                    feasibility at each state k = 1,...,K.                 
% %         
% % Email: mortheus__at__gmail.com
% % =======================================================================
% % $Revision: 4.0.3 $  $Date: 2011/05/02 11:38:20 $  

    K = size(D,1);  % equals K
    N_polComb = size(polComb,1);    % #{ b }
    N_Z = size(polComb,2);          % #{ D(Z) }
    
    polFeasible = -Inf*ones(K,N_polComb,N_Z);
                                    % Intermediate storage of feasible b's

    B = cell(K,1);                  % For storing feasible b correspondence
    
    for k = 1 : K
       
        lambda_k = D(k,:);
        
        for b = 1 : N_polComb          
            
            b_k_comb = polComb(b, :);
            
            budget_balance = lambda_k * b_k_comb';
            
             if budget_balance >= 0
                 polFeasible(k, b, :) = b_k_comb;
             end      
        end
        
        % Discard infeasible policies at each state k:
        bf_set_now = polFeasible(k,:,:);
        bf_set_now = bf_set_now( bf_set_now ~= -Inf ); 
        
        % Store feasible b set for each k (variable dimensions!) as cell:
        B{k} = reshape( bf_set_now, length(bf_set_now)/N_Z, N_Z);
    end
    
end

