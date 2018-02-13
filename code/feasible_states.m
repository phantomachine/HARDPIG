function F = feasible_states(self,xPolyQjRand,xPolyIndexQjRand)
    
    % FEASIBLE_STATES.M
    %
    % INPUTS: 
    % 
    % self:
    %       Object containing model from KSMOD03
    %
    % xPolyQjRand{a}{k}:
    %       a = 1:NA, 
    %       k = 1:K,  
    % gives an (Nsim x 3) array of probability distributions randomly drawn
    % from each current partition element Q(k) of the simplex D.
    %
    % xPolyIndexQjRand{a}{k}{i}: 
    %       where 
    %       i = 1,...,I, where I = numel(xPolyRand{a}{k}), for each (a,k)
    % gives index of random vectors X(k) in Q(k) that maps to vectors
    % in subsets indexed by j(i) :=: P(a)[ Q(k) ]. This gives a useful
    % relation between each random vector and where its continuation state
    % under mapping P(a) will end up in---i.e. the partition Q(j).
    %
    % OUTPUT:
    %
    % F: A structure variable such that
    %   F{a}{k}{i}.Xi: 
    %           Xi subset of X(k) that will map to Q(j[i]) under a
    %   F{a}{k}{i}.balance:
    %           Budget balances: Xi*(self.IndexProfileSetB)'
    %   F{a}{k}{i}.rc
    %           rc = [ row, col ] indexes to nonnegative "balances".
    % Type: e.g. a = 2; k = 5; i = numel(F{a}{k}); F{a}{k}{i} to see
    % example content.
    % 
    % (c) 2012, 2013, Timothy Kam. Email: tcy.kam__at__gmail.com 
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Use subject to GNU LGPL licensing terms. Cite this header and 
    % author completely in subsequent re-use and modifications.
    % =====================================================================
    % $Revision: 5.0.0 $  $Date: 2014/01/03 12:45:20 $ 
    
    % Ad-hoc setting for tolerance for positive budget balance:
    TOL_pos = 0.4;
    
    A = self.ProfileSetA;
    B = self.ProfileSetB;

    [NA,~] = size(A);
    K = numel(xPolyQjRand{1});

    F = cell(NA,1);
    feasiblecorr = cell(K,1);

    for na = 1 : NA % Loop over current a profile to get each P(a)
        
        for k = 1 : K % loop over current state-space partition element
    
           % Random vectors in Q(k) triangle:
           X = xPolyQjRand{na}{k};
           
           I = numel(xPolyIndexQjRand{na}{k});
           
           feasible_sub = cell(1,I);
           
                for i = 1 : I % Loop over partition j that intersects k
                    
                    % Set of random vectors in Q(k) that goes into j:
                    selection = xPolyIndexQjRand{na}{k}{i};
                   
                        
                    if isempty(selection)
                        Xi = [];
                        balance = [];
                        row = [];
                        col = [];
                    else
                        
                        Xi = X(selection,:);
                        balance = -Xi*(B');
                        [row,col] = find(balance >= -0.4 & balance < TOL_pos);
                        
                        % Notes:
                        %   row: Also row number of Xi
                        %   col: Also row number of B
                        if isempty(row) || isempty(col)
                            Xi = [];
                            balance = [];
                            row = [];
                            col = [];
                        end
                    end
                    
                    feasible_sub{i}.Xi = Xi;
                    feasible_sub{i}.balance = balance;
                    feasible_sub{i}.ind_Xi = row;
                    feasible_sub{i}.ind_B = col;
                end % End i
                
           feasiblecorr{k} = feasible_sub;
           
        end % End k
        
        F{na} = feasiblecorr;
        
    end % End na
    
end

    

