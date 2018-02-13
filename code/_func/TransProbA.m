function P_a = TransProbA(model,parallel)    

% Set of transition probabilities P(A) from Set of A:
    
% % TransProbA.M
% %
% % Function computes transition probabilities from all action profiles a 
% % in A^N_Z, where N_Z = M + N is number of private agent states.
% % 
% % =======================================================================
% %     (c) 2011-- T.Kam and R.Stauber.
% %
% % INPUT: Let 
% %
% %     * model      : KSMOD2 class
% %     * NA         : No. Combinations for all a's: Use GRIDMAKE
% %     * A_u        : Set of all vectors [a(-1),...,a(-N)], Dim: NA x N
% %     * A_e        : Set of all vectors [a(+1),...,a(+M)], Dim: NA x N
% % 
% % OUTPUT
% % 
% %     * P_a        : A numeric array of dimension (N_Z x N_Z x NA) 
% %                                            
% % Email: mortheus__at__gmail.com
% % =======================================================================
% % $Revision: 4.0.3 $  $Date: 2011/05/02 11:38:20 $  
% %
% % See also KSMOD2, MARKOVMAT (in KSMOD2)

    N_Z = model.NS;                     % Total number of agent states j
    %M   = model.M;                      % Max(j>0)
    N   = model.N;                      % Max(j<0)
         
    A_u = model.ProfileSetA(:,1:N);     % action profiles for j < 0
    A_e = model.ProfileSetA(:,N+1:end); % action profiles for j > 0
    
    NA  = size(A_e,1);                  % Number of action profiles in A
    
    P_a = zeros(N_Z,N_Z,NA);            % Preallocation space of Markov 
                                        % matrices P(a), a \in A
    
    if strcmp(parallel,'on')
        
        %PoolSize = feature('numCores');
        %matlabpool(PoolSize);
        
        parfor a_idx = 1 : NA
            au = A_u(a_idx,:); % Note: au = [ a(-1), ..., a(-N) ]
            ae = A_e(a_idx,:); % Note: ae = [ a(+1), ..., a(+M) ]

            P_a(:,:,a_idx) = markovmat(model,au,ae);
        end
        
        %matlabpool close
    else
        for a_idx = 1 : NA
            au = A_u(a_idx,:); % Note: au = [ a(-1), ..., a(-N) ]
            ae = A_e(a_idx,:); % Note: ae = [ a(+1), ..., a(+M) ]

            P_a(:,:,a_idx) = markovmat(model,au,ae);
        end
    end
        