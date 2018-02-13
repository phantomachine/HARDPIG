% % blp_randtest.m
 clear
clc

yalmip('clear')
s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s);

%% EXAMPLE 1

N_global =10;
L_global = 100;
M_global =100;


% h = randn(1,N);
% P = [ 0.2, 0.8, 0;
%       0.5, 0,  0.5;
%       0.3, 0,  0.7 ];
%   
% A = rand(L,N); 
% b = rand(L,1);
% 
% B = rand(M,N);
% c = rand(M,1);
% 
% E = rand(N,1);
% F = rand(1);
% 
% delta = 0.8;
% 
% w = sdpvar(N,1);
% l = sdpvar(N,1);
% 
% objective = h*w;%  + delta*(l')*rand(N,N)*w ;
% 
% a = zeros(1,100);

%distcomp.feature( 'LocalUseMpiexec', false ); % Speed up matlabpool (R2009)
tic
%parfor i =  1 : 1
for i = 1:1
    
    L = L_global;
    M = M_global;
    N = N_global;
    
    h = randn(1,N);
    P = [ 0.2, 0.8, 0;
          0.5, 0,  0.5;
          0.3, 0,  0.7 ];

    A = rand(L,N); 
    b = rand(L,1);

    B = rand(M,N);
    c = rand(M,1);

    E = rand(N,1);
    F = rand(1);

    delta = 0.8;

    w = sdpvar(N,1);
    l = sdpvar(N,1);

    objective = h*w;%  + delta*(l')*rand(N,N)*w ;

    a = zeros(1,100);
    
    
    constraint = [ (l')*E + delta*(l')*w <=F, A*w <= b, B*l <= c, ...
                        -1000000 <= w(:) <= 100000, 0 <= l(:) <= 1 ];

    a(i) = 1;
    
    sdp_options = sdpsettings('verbose',0,...
                                    'solver','bmibnb',...
                                    'bmibnb.roottight',0|1,...
                                    'bmibnb.numglobal', 10,...
                                    'bmibnb.uppersolver','snopt',...
                                    'bmibnb.lowersolver','glpk');


    sol = solvesdp(constraint,-objective,sdp_options);
end
toc

% 
% %% EXAMPLE 2
% %
% sdpvar x y
% F = [x^3+y^5 <= 5, y >= 0];
% F = [F, -100 <= [x y] <= 100]; % Always bounded domain
% options = sdpsettings('verbose',1,'solver','snopt');
% solvesdp(F,-x,options)

