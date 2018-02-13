function [ distance ] = mdistance(A,B)

% mdistance.m
%-------------------------------------------------------------------------%
%
% Input: matrices (A, B) of (L x n) dimensions each
% Output: scalar on R+
% Note: partitions on state space j = 1, ...,n
%       number of vertices, l = 1,...,L
%
%-------------------------------------------------------------------------%
% Coded by Timothy Kam. Copyleft 2007. Use and abuse subject to the spirit
% of the GNU GPL. See: http://www.fsf.org/licensing/licenses/gpl.html
%-------------------------------------------------------------------------%



% distance(W,W+) = max_j max_l | c(l,j) - c_{+}(l,j) |
    distance = max(max(abs(A - B)));