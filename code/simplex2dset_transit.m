function PT = simplex2dset_transit(P, T)

% SIMPLEX2DSET_TRANSIT.M
%
% Let $D := \Delta(Z)$ be a 2D unit-simplex. Find set-valued image 
% of linear map, 
%
%   $$P: \Delta(Z) \rightarrow \Delta(Z);
%

% and $P$ is a stochastic/Markov matrix linear operator. Image implicitly
% defined by extreme points of $P(T)$.
%
% Example: $T \in D$ is convex polytope (e.g. triangle), and $P(T)$ is also 
% a convex-polytope subset of $\Delta (Z)$. 
% 
%
% Input:
%
% T   :    Set of extreme points of K convex subsets of $\Delta(Z)$.
%          An (N_ExtPt x K) x N_Z array.
%               * N_ExtPt:  Each $T_{k}$ has N_ExtPt (== 3 for triangles)
%                           extreme points.
%               * K      :  Number of partition elements of P
%               * N_Z    :  dim(D) (e.g. 2, for 2D simplex)
%
% P   :    Markov matrix, N_Z x N_Z array.
%
% Output:
%
% PT  :    Image of all extreme points P(T), same size as T.
%
% 
% (c) 2012, 2013, Timothy Kam. Email: tcy.kam__at__gmail.com 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% WARNING: Use subject to GNU LGPL licensing terms. Cite this header and 
% author completely in subsequent re-use and modifications.
% =======================================================================
% $Revision: 5.0.0 $  $Date: 2013/09/11 12:45:20 $ 
%
% See also SIMPLEX2DSET_INTERSECT

PT = T * P;
