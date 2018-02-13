function U = verts_solve(H,c)

% % VERTS_SOLVE.M
% %
% % Given (L x N) matrix of directional vectors H, and (L x 1) vector of 
% % levels, c, solve for points x^l:= (x(l), ...,x(l+N)) such that 
% %
% %         H(l,:)*x^l = c(l:l+N-1),
% %
% % for all l = 1,...,L.
% %
% %     (c) 2011, T. Kam. Email: timothy.kam@anu.edu.au
% % =======================================================================
% % $Revision: 1.0.0 $  $Date: 2011/11/14 03:32:00 $ 

    [L, N] = size(H);
    
    U = zeros(L,N);
    
%     L_up = L-N+1;
%     
%     for l_idx = 1:L
%         if (l_idx < L_up)
%           l_idx_plus = l_idx + (N-1);
%           
%           X = H(l_idx:l_idx_plus,:);
%           
%           y = c(l_idx:l_idx_plus); 
%           
%         else
%           X = [   H(l_idx:L,:)          ; 
%                   H(1:N-(L-l_idx),:)     ];
%               
%           y = [   c(l_idx:L)          ; 
%                   c(1:N-(L-l_idx))     ]; 
%         end
%         
%         U(l_idx,:) = (X\y)';
%         
%     end
%     
% end

% % EXAMPLE USAGE:
% %
% % Recursively solve for ADP pairs (U_1, U_2) from two linear
% % equations:
% %
% %     [ h_1      ] * [ U_1, U_2 ]' = [ c_old(l)      ]
% %     [ h_l_plus ]                   [ c_old(l_plus) ]
% % 
% % for all l = 1,...,L-1; for final pair, set l = L and l_plus =
% % 1, to get two equations.
% 
          for l_idx = 1:L
            if (l_idx < L)
              l_idx_plus = l_idx + 1;
            else
              l_idx_plus = 1;
            end

            X = [   H(l_idx,:)          ; 
                    H(l_idx_plus,:)     ];

            y = [   c(l_idx)        ; 
                    c(l_idx_plus)   ];

            U(l_idx,:) = (mldivide(X,y))';%(invpd(X)*y)'; 

          end