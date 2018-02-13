% A LP example which shows all the potentials of GLPKMEX
clear;

disp('LP problem');

% c=[10,6,4]';
% a=[1,1,1;...
%    10,4,5;...
%    2,2,6];
% b=[100,600,300]';
% ctype=['U','U','U']';
% lb=[0,0,0]';
% ub=[]';
% vartype=['C','C','C']';

%c=[10,6,4,6,9,20,10,30,40,20,1,3]';

c = rand(5,1);
%c = c(1:5);

% a=[1,1,1,1,1,1,1,1,1,1,1,1;...
%    10,4,5,5,3,8,1,3,2,1,1,2;...
%    2,2,6,.4,.5,3,1,7,5,4,3,5;
%    10,4,5,5,3,3,7,1,9,0,2,4;
%    1,1,4,2,7,5,8,3,6,3,1,2;
%    1,4,3,8,6,5,0,6,9,9,0,6;
%    3,3,3,5,6,4,5,4,9,0,0,0;
%    1,1,1,1,1,1,5,3,7,8,9,0;
%    2,2,3,4,6,7,8,1,3,6,5,4;
%    1,1,3,3,4,5,6,2,8,34,2,1;
%    1,34,5,6,7,3,7,9,6,5,3,12;
%    4,4,7,8,6,3,6,7,8,9,2,4;];

a = rand(300,5);

%a = a(1:5,1:5);

%b=[100,600,300,34,290,111,32,456,7,89,23,126]';

b = rand(300,1);

%b = b(1:5);

ctype=repmat('U',300,1);%['U','U','U','U','U','U','U','U','U','U','U','U']';

%ctype = ctype(1:5);

%lb=[0,0,0,0,0,0,0,0,0,0,0,0]';

%lb = lb(1:5);

lb = zeros(5,1);

ub=[]';

vartype=repmat('C',5,1); %['C','C','C','C','C','C','C','C','C','C','C','C']';

%vartype = vartype(1:5);

% Output all GLPK messages on workspace
param.msglev=1;
% Set save options
param.save=0;
param.savefilename='SimpleLP';
param.savefiletype='fixedmps';

% profile -memory off



% for j = 1:5
%     %matlabpool open 8
%     parfor i = 1:100
    s=-1;
    [xmin,fmin,status,extra]=glpk(c,a,b,lb,ub,ctype,vartype,s,param);

    s = rand(1,10);
    S = rand(100,10);

    index = getindex(s,S);

%     end
%     %matlabpool close
% end

% profile viewer
% profsave(profile('info'), 'profiling_results')

% % Compare with linprog:
% options = optimset('Display', 'on','LargeScale', 'on', 'Simplex', 'on',...
%                     'UseParallel','never');
% tic
% [Xmin,Fmin,flag] = linprog(-c,a,b,[],[],lb,ub,[],options)
% toc

%lpsolver = param.lpsolver;
%save_pb = param.save;
%[xmin,fmin,status,extra]=glpkmex(s,c,a,b,ctype,lb,ub,vartype,param,lpsolver,save_pb)
