function p = trirod( diameter, height, varargin )
%  TRIROD - Faces and vertices for rod-shaped particle.
%
%  Usage :
%    p = trirod( diameter, height,                 varargin )
%    p = trirod( diameter, height, n,              varargin )
%    p = trirod( diameter, height, n, 'triangles', varargin )
%  Input
%    diameter     :  diameter of rod
%    height       :  total height (length) of rod
%    n            :  number of discretization points [ nphi, ntheta, nz ]
%    'triangles'  :  use triangles rather than quadrilaterals
%    varargin     :  additional arguments to be passed to PARTICLE
%  Output
%    p            :  faces and vertices of triangulated rod

%  extract number of discretization points
if ~isempty( varargin ) && isnumeric( varargin{ 1 } )
  [ n, varargin ] = deal( varargin{ 1 }, varargin( 2 : end ) );
  assert( numel( n ) == 3 );
else
  n = [ 15, 20, 20 ];
end

%  angles 
phi   = linspace( 0, 0.5 * pi, n( 1 ) );
theta = linspace( 0, 0.5 * pi, n( 2 ) );
%  z-values of cylinder
z = 0.5 * linspace( - 1, 1, n( 3 ) ) * ( height - diameter );

%  upper cap
cap1 = shift( trispheresegment( phi, theta, diameter, varargin{ : } ),  ...
                              [ 0, 0, 0.5 * ( height - diameter ) ] );
%  lower cap
cap2 = flip( cap1, 3 );

%  grid for cylinder discretization
[ verts, faces ] = fvgrid( phi, z, varargin{ : } );
%  cylinder coordinates
[ phi, z ] = deal( verts( :, 1 ), verts( :, 2 ) );
%  make cylinder
x = 0.5 * diameter * cos( phi );
y = 0.5 * diameter * sin( phi );


%  remove TRIANGLES argument from varargin
if ~isempty( varargin ) && ischar( varargin{ 1 } )  ...
                        && strcmp( varargin{ 1 }, 'triangles' )
  varargin = varargin( 2 : end );
end
%  cylinder particle
cyl = particle( [ x, y, z ], faces, varargin{ : } );
cyl2 = rot(cyl,90,[1,0,0]);
cyl2 = shift(cyl2,[0,-(height-diameter)/2,(height-diameter)/2]);
cyl3 = rot(cyl,90,[0,0,-1]);
cyl3 = rot(cyl3,90,[1,0,0]);
cyl3 = shift(cyl3,[0,-(height-diameter)/2,-(height-diameter)/2]);

cyl4 = rot(cyl,90,[0,-1,0]);
cyl4 = shift(cyl4,[-(height-diameter)/2,0,(height-diameter)/2]);
cyl5 = rot(cyl,90,[0,0,1]);
cyl5 = rot(cyl5,90,[0,-1,0]);
cyl5 = shift(cyl5,[-(height-diameter)/2,0,-(height-diameter)/2]);

cyl6 = flip(cyl,1);
cyl6  = shift(cyl6 ,[-(height-diameter),0,0]);

cyl7 = flip(cyl,2);
cyl7  = shift(cyl7,[0,-(height-diameter),0]);

cyl8 = flip(cyl,1);
cyl8 = flip(cyl8,2);
cyl8 = shift(cyl8,[-(height-diameter),-(height-diameter),0]);

cyl9 = flip(cyl2,1);
cyl9  = shift(cyl9,[-(height-diameter),0,0]);

cyl10 = flip(cyl3,1);
cyl10  = shift(cyl10,[-(height-diameter),0,0]);

cyl11 = flip(cyl4,2);
cyl11 = shift(cyl11,[0,-(height-diameter),0]);

cyl12 = flip(cyl5,2);
cyl12 = shift(cyl12,[0,-(height-diameter),0]);

cap3 = flip(cap1,1);
cap3 = shift(cap3,[-(height-diameter),0,0]);

cap4 = flip(cap1,2);
cap4 = shift(cap4,[0,-(height-diameter),0]);

cap5 = flip(cap1,1);
cap5 = flip(cap5,2);
cap5 = shift(cap5,[-(height-diameter),-(height-diameter),0]);

cap6 = flip(cap2,1);
cap6 = shift(cap6,[-(height-diameter),0,0]);

cap7 = flip(cap2,2);
cap7 = shift(cap7,[0,-(height-diameter),0]);

cap8 = flip(cap2,1);
cap8 = flip(cap8,2);
cap8 = shift(cap8,[-(height-diameter),-(height-diameter),0]);

val = (height-diameter)*linspace(-0.5,0.5,n(3));
[ verts_,faces_] = fvgrid( val, val);
pl = particle(verts_,faces_);

pl1 = shift(pl,[-(height-diameter)/2,-(height-diameter)/2,(height)/2]);
pl2 = flipfaces(shift(pl,[-(height-diameter)/2,-(height-diameter)/2,-(height)/2]));

pl3 = rot(pl,90,[1,0,0]);
pl3 = shift(pl3,[-(height-diameter)/2, -height+diameter/2, 0]);

pl4 = rot(pl,90,[1,0,0]);
pl4 = flipfaces(shift(pl4,[-(height-diameter)/2, diameter/2, 0]));

pl5 = rot(pl,90,[0,1,0]);
pl5 = flipfaces(shift(pl5,[-height+diameter/2,-(height-diameter)/2,  0]));

pl6 = rot(pl,90,[0,1,0]);
pl6 = shift(pl6,[diameter/2, -(height-diameter)/2, 0]);

%  compose particle

p = clean( [cap1;cap2;cap3;cap4;cap5;cap6;cap7;cap8;...
            cyl;cyl2;cyl3;cyl4;cyl5;cyl6;cyl7;cyl8;cyl9;...
            cyl10;cyl11;cyl12;pl1;pl2;pl3;pl4;pl5;pl6] );
