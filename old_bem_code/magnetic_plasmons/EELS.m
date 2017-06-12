addpath(genpath('/data/home/dmasiell/code/MNPBEM14'));
op = bemoptions( 'sim', 'ret', 'interp', 'curv', 'eels.refine', 2 );

 %  table of dielectric functions
epstab  = { epsconst( 1 ), epstable( 'silver.dat' )};

%  X-LENGTH, Y-LENGTH OF RODS
rodlen_left  = [ 60, 215];
rodlen_right = [ 60, 215];
height = 20;
h=0.25;

% HEIGHT OF RODS
rodedge = edgeprofile( height, 10, 'mode', '00', 'min', 0 );

% MAKE ROD POLYGON
rod_left  = polygon( 30, 'size', rodlen_left  ); 
rod_right = polygon( 30, 'size', rodlen_right );

% RADIUS OF ONE RING
radius1 = 90;
radius2 = 100;

% LOCATION OF CENTER OF ONE RING
center = 160;
del = 17*pi/180;
distx = 19;
disty = 14;

% ANGLES
theta = [0, pi/3, 2*pi/3, pi, 4*pi/3, 5*pi/3];

hdata = struct( 'hmax', 8 );

[ p, poly ] = tripolygon( rod_left, rodedge);
[ pr, polyr ] = tripolygon( rod_right, rodedge);

p4 = shift(p, [radius2*cos(theta(4))-center+distx, radius2*sin(theta(4)), h]);
p5 = shift(rot(p, 60-del*180/pi), [radius2*cos(theta(2))-center, radius2*sin(theta(2))+disty, h]);
p6 = shift(rot(p, 120+del*180/pi), [radius2*cos(theta(6))-center, radius2*sin(theta(6))-disty, h]);

p = comparticle( epstab, { p4, p5, p6},  ... 
 [ 2,1; 2,1; 2,1 ], 1,2,3, op );
plot(p);
%%
units;
[ width, vel ] = deal( 0.2, eelsbase.ene2vel( 200e3 ) );
imp = [-259,120;-234,127];
enei = eV2nm./linspace( NUM,NUM,1);
ene = transpose(eV2nm./enei);

bem = bemsolver( p, op );
exc = electronbeam( p, imp, width, vel, op );
[ psurf, pbulk ] = deal( zeros( size( imp, 1 ), numel( enei ) ) );

for ien = 1 : length( enei )
  sig = bem \ exc( enei( ien ) );
  [ psurf( :,ien ), pbulk( :,ien ) ] = exc.loss( sig );
end

fid = fopen(['Spectrum'],'wt');
for j = 1:length(ene)
          fprintf(fid,' %g', ene(j));
          fprintf(fid,' %g', psurf(1,j));
          fprintf(fid,' %g', psurf(2,j));
          fprintf(fid, '\n');
end
fclose(fid)
