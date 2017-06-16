clear, clc
op = bemoptions( 'sim', 'ret', 'interp', 'curv' );

epstab  = { epsconst( 1 ), epstable( 'silver.dat' )};

rodlen_left  = [ 60, 215];
rodlen_right = [ 60, 200];
height = 20;
h=0.25;

rodedge = edgeprofile( height, 10, 'mode', '00', 'min', 0 );

rod_left  = polygon( 60, 'size', rodlen_left  ); 
rod_right = polygon( 60, 'size', rodlen_right );

radius1 = 84;
radius2 = 90;

center = 153;

% ANGLES
theta = [0, pi/3, 2*pi/3, pi, 4*pi/3, 5*pi/3];
%  hdata
hdata = struct( 'hmax', 8 );
%  extrude rod
[ p, poly ] = tripolygon( rod_left, rodedge);
[ pr, polyr ] = tripolygon( rod_right, rodedge);
% extrude substrate
%brick = tripolygon( subs, subedge);
% MOVE AND ROTATE RODS
p1 = shift(pr, [radius1*cos(theta(1))+center, radius1*sin(theta(1)), h]);
p2 = shift(rot(pr, 120), [radius1*cos(theta(3))+center, radius1*sin(theta(3)), h]);
p3 = shift(rot(pr, 60), [radius1*cos(theta(5))+center, radius1*sin(theta(5)), h]);

p4 = shift(p, [radius2*cos(theta(4))-center, radius2*sin(theta(4)), h]);
p5 = shift(rot(p, 60), [radius2*cos(theta(2))-center, radius2*sin(theta(2)), h]);
p6 = shift(rot(p, 120), [radius2*cos(theta(6))-center, radius2*sin(theta(6)), h]);

p = comparticle( epstab, {  p1, p2, p3, p4, p5, p6 },  ... 
    [ 2,1; 2,1; 2,1; 2,1; 2,1; 2,1 ], 1,2,3,4,5,6,  op );

dir = [  1, 0, 0 ];
pol = [  0, 1, 0 ];


%  points where electric fields are computed
[ x, y ] = meshgrid(  linspace( -300, 300, 301 ), linspace( -200, 200, 201 ) );
z = 0 * x + 10; %At this line I multiplied each component by 0 to get the field at 0.

%  place the points into the dielectric media
pt = compoint( p, [ x( : ), y( : ), z( : ) ], 'mindist', 0.5 );
g = compgreenret( pt, p );

units;
ene = [1.32,1.38];
enei = eV2nm./ene;
bem = bemsolver( p, op );
exc = planewave( pol, dir, op );

sca = zeros( numel( enei ), size( dir, 1 ) );
ext = zeros( numel( enei ), size( dir, 1 ) );

for ien = 1 : length( enei )
	    sig = bem \ exc( p, enei( ien ) );

field=g.field(sig);
Alicex=pt(real(field.e(:,1)));
Alicey=pt(real(field.e(:,2)));
Alicez=pt(real(field.e(:,3)));
Alice=sqrt(Alicex.^2+Alicey.^2+Alicez.^2);
Alice = reshape(Alice,size(x));
Bobz=pt(real(field.h(:,3)));
Bob = reshape(Bobz,size(x));

save(strcat('200E',num2str(ene(ien))),'Alice','-ascii')
save(strcat('200B',num2str(ene(ien))),'Bob','-ascii')
ext( ien, : ) = exc.ext( sig );
sca( ien, : ) = exc.sca( sig );
end

