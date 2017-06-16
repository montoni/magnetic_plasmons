function p = wavytorus( diameter, rad, Amplitude, n_wave, varargin )
%  WAVYTORUS - Faces and vertices of triangulated wavy-torus.
%
%  Usage :
%    p = wavytorus( diameter, rad, Amplitude, n_wave, varargin )
%    p = wavytorus( diameter, rad, Amplitude, n_wave, n, varargin )
%    p = wavytorus( diameter, rad, Amplitude, n_wave, n, 'triangles', varargin )
%  Input
%    diameter     :  diameter of folded cylinder
%    rad          :  radius of torus
%    Amplitude    :  amplitude of oscillations for the torus
%    n_wave       :  number of 'waves' i.e. sin(n_wave*phi)
%    n            :  number of discretization points
%    'triangles'  :  use triangles rather than quadrilaterals
%    varargin     :  additional arguments to be passed to PARTICLE
%  Output
%    p            :  faces and vertices of triangulated torus

%  extract number of discretization points
if ~isempty( varargin ) && isnumeric( varargin{ 1 } )
  [ n, varargin ] = deal( varargin{ 1 }, varargin( 2 : end ) );
  if numel( n ) == 1,  n = [ n, n ];  end
else
  n = [ 21, 21 ];
end
%  extract triangle keyword
if ~isempty( varargin ) && ischar( varargin{ 1 } )  ...
                        && strcmp( varargin{ 1 }, 'triangles' )
  [ triangles, varargin ] = deal( varargin{ 1 }, varargin( 2 : end ) );
else
  triangles = '';
end

%  grid triangulation
[ verts, faces ] = fvgrid( linspace( 0, 2 * pi, n( 1 ) ),  ...
                           linspace( 0, 2 * pi, n( 2 ) ), triangles );
%  angles
[ phi, theta ] = deal( verts( :, 1 ), verts( :, 2 ) );
%  coordinates of torus
x = ( 0.5 * diameter + Amplitude* sin(n_wave*phi) + rad  * cos( theta ) ) .* cos( phi );
y = ( 0.5 * diameter + Amplitude* sin(n_wave*phi) + rad  * cos( theta ) ) .* sin( phi );
z =                    rad * sin( theta );
% make torus
p = clean( particle( [ x, y, z ], faces, varargin{ : } ) );
