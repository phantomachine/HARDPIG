function [A,NA,N,A_u,A_e] = A_ProfileSet(self)

% A_PROFILESET.M:
% 
% Get action profile space A
% Sort by actions of
%   A_u : -ve state agents i = 1,...,N
%   A_e : +ve state agents i = N+1,...,M
%
% See also SPACE_PRIVATE_ACTIONS, KSMOD03

        % Space of private agents' action profiles: 
        A = self.ProfileSetA;
        NA = size(A,1);
        N = self.N;
        
        % DEPRECATED (1/11/2013, see SPACE_PRIVATE_ACTIONS, for this old
        % feature)
        % Add index to last column of A:
        % A = [ A, (1:NA)'];

        A_u = A(:,1:N);      
        A_e = A(:,N+1:end);