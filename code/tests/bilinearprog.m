function [] = bilinearprog(params, x0, w0, options)

%% BILINEARPROG.M
%
% Solve a separable bilinear program of the form:
%
% $$\min_{x | y } + r_{1}'x +  x'Cy + r_{2}'y $$
%
% s.t.
%
% $$A_{1} x \leq c_{1},  \qquad | \qquad A_{2}y \leq c_{2},$$
%
% $$x \geq \mathbf{0}, \qquad | \qquad z \geq \mathbf{0}.$$
%
% Solution is guaranteed to converge in finite steps. Solution is only
% guaranteed to be local optimum. For global optimum, we need to use either
% Konno's (1974) Cutting Plane algorithm or McCormick's Branch-and-Bound 
% algorithm (see e.g. YALMIP's BMIBNB routine).
%
% (c) 2013, T.Kam
%
% See also GLPK, LINPROG

%% 0. EXTRACT PARAMETERS
%
% Objective function:

s1 = params.s1;     % w
s2 = params.s2;     % z
C  = params.c;      % x'Cy
r1 = params.r1;     % x
r2 = params.r2;     % y

%%
% Constraint set 1 for (x,w):

A1 = params.a1;     % x
B1 = params.b1;     % w
c1 = params.c1;     % constants

%%
% Constraint set 2 for (y,z):

A2 = params.a2;     % y
B2 = params.b2;     % z
c2 = params.c2;     % constants

%%
% Convergence settings:

MAXITER = options.maxiter;
TOL = options.tolerance;





