function [A,NA,N,A_u,A_e] = space_private_actions(self)

        % Space of private agents' action profiles: 
        A = self.AGRID;
        NA = size(A,1);
        N = self.N;
        
        % Add index to last column of A:
        A = [ A, (1:NA)'];

        A_u = A(:,1:N);      
        A_e = A(:,N+1:end-1);