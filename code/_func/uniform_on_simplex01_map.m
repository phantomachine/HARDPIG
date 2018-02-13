function [ x, seed ] = uniform_on_simplex01_map ( m, n, seed )

%*****************************************************************************80
%
%% UNIFORM_ON_SIMPLEX01_MAP maps uniform points onto the unit simplex.
%
%  Discussion:
%
%    The surface of the unit M-dimensional simplex is the set of points
%    X(1:M) such that each X(I) is nonnegative,
%    every X(I) is no greater than 1, and
%
%    ( X(I) = 0 for some I, or sum ( X(1:M) ) = 1. )
%
%    In M dimensions, there are M sides, and one main face.
%    This code picks a point uniformly with respect to "area".
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license. 
%
%  Modified:
%
%    20 April 2013
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Reuven Rubinstein,
%    Monte Carlo Optimization, Simulation, and Sensitivity
%    of Queueing Networks,
%    Wiley, 1986.
%
%  Parameters:
%
%    Input, integer M, the dimension of the space.
%
%    Input, integer N, the number of points.
%
%    Input/output, integer SEED, a seed for the random number generator.
%
%    Output, real X(M,N), the points.
%

%
%  The construction begins by sampling M points from the
%  exponential distribution with parameter 1.
%

Nburn = n^3; % Hack by T. Kam

x = zeros(m,n + Nburn);

  for j = 1 : n + Nburn % Hack by T. Kam

    [ e, seed ] = r8vec_uniform_01 ( m, seed );

    e(1:m,1) = - log ( e(1:m,1) );

    e_sum = sum ( e(1:m,1) );

    x(1:m,j) = e(1:m,1) / e_sum;
%
%  The point X is now on the "main" face of the unit simplex.
%
%  Based on their relative areas, choose a side of the simplex,
%  or the main face.
%
    area1 = sqrt ( m );

    area2 = m;

    [ r, seed ] = r8_uniform_01 ( seed );
%
%  If we choose to move the point from the main face,
%  set a random coordinate to 0.
%
    if ( area1 / ( area1 + area2 ) < r )
      [ i, seed ] = i4_uniform_ab ( 1, m, seed );
      x(i,j) = 0.0;
    else

    end

  end
  
  % Hack by T. Kam
  y = x(:,sum(x,1) == 1);
  
  if size(y,2) >= n
    x = y(:,1:n);
  else
    x = x./repmat(sum(x,1),m,1);
    x = x(:,1:n);
  end
    

  return
end