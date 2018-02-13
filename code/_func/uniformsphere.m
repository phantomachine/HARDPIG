function U = uniformsphere(Npoints, M)

% UNIFORMSPHERE.M
%
% Generate Npoints number of random unit (normal) vectors that are 
% uniformly distributed on the hypershere embedded in R^M Euclidean space. 
% Input X is (Npoints x M);

        X = randn(Npoints,M);
        normX = sqrt(sum(X.^2, 2));
        normX = repmat(normX,1,M);
        U = X ./ normX;
        
end