%%  initialization
%  options for BEM simulation
clear, clc, clf

%  options for BEM simulation
op = bemoptions( 'sim', 'ret', 'interp', 'curv' );

 %  table of dielectric functions
epstab  = { epsconst( 1 ), epstable( 'silver.dat' ), epsconst( 2 ),  };
%  X-LENGTH, Y-LENGTH OF RODS


for i=0:5
a=1+i/4;    
w0=25;
width=w0*a;
scale=3;
le=scale*width;
h0=2;

rodlen_left  = [ w0, le];
rodlen_right = [ w0, le];
height = 14;
h=0;

% HEIGHT OF RODS
rodedge = edgeprofile( height, 10, 'mode', '00', 'min', 0 );

% X-LENGTH, Y-LENGTH OF SUBSTRATE
sublen = [ 900, 700];

% MAKE ROD POLYGON
rod_left  = polygon( 25, 'size', rodlen_left  ); 
rod_right = polygon( 25, 'size', rodlen_right );


% MAKE ROD POLYGON
% p  = trirod( width, le, [12,12,12]  ); 
% pr = trirod( width, le, [12,12,12]  );
% pr= rot(pr,90,[0,1,0]);


% RADIUS OF ONE RING
space=a;
radius1 = (space+le/2+width/2)/sqrt(3)-h0;
radius2 = (space+le/2+width/2)/sqrt(3)-h0;

% LOCATION OF CENTER OF ONE RING
center = 0;

% ANGLES
%theta = [0, pi/3, 2*pi/3, pi, 4*pi/3, 5*pi/3];

% MOVE AND ROTATE RODS
% p1 = shift(rot(pr,90)  , [radius1*cos(theta(1))+center, radius1*sin(theta(1)), h]);
% p2 = shift(rot(pr, 210), [radius1*cos(theta(3))+center, radius1*sin(theta(3)), h]);
% p3 = shift(rot(pr, 150), [radius1*cos(theta(5))+center, radius1*sin(theta(5)), h]);
% 
% p4 = shift(rot(p,90), [radius2*cos(theta(4))-center, radius2*sin(theta(4)), h]);
% p5 = shift(rot(p, 150), [radius2*cos(theta(2))-center, radius2*sin(theta(2)), h]);
% p6 = shift(rot(p, 210), [radius2*cos(theta(6))-center, radius2*sin(theta(6)), h]);

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

 p = comparticle( epstab, {  p1, p2, p3 },  ... 
  [ 2,1; 2,1; 2,1 ], 1,2,3, op );

plot(p)

%p = comparticle( epstab, {  p1, p2, p3 },  ... 
 %   [ 2,1; 2,1; 2,1 ], 1,2,3,  op );
%figure(1)
%plot(p, 'EdgeColor','b')
%  EELS excitation     
units;
%  width of electron beam and electron velocity 
[ width, vel ] = deal( 0.1, eelsbase.ene2vel( 200e3 ) );
%  impact parameters
%imp = [ 315, 95 ];
  imp = [ -(1.5)*radius1, -width/2 ];
%  loss energies in eV
enei = eV2nm./linspace( 1, 3.5, 91);
ene = transpose(eV2nm./enei);
%  BEM simulation
%  BEM solver
bem = bemsolver( p, op );
 
%  electron beam excitation
exc = electronbeam( p, imp, width, vel, op );
%  surface and bulk loss
[ psurf, pbulk ] = deal( zeros( size( imp, 1 ), numel( enei ) ) );
%
multiWaitbar( 'BEM solver', 0, 'Color', 'g', 'CanCancel', 'on' );
%  loop over energies
tic 
%parpool(6)
for ien = 1 : length( enei )
  %tic
  %  surface charges
  sig = bem \ exc( enei( ien ) );
  %  EELS losses
  [ psurf( :,ien ), pbulk( :,ien ) ] = exc.loss( sig );
  %eV2nm./enei(ien)
  %psurf( :, ien )
  %toc
  multiWaitbar( 'BEM solver', ien / numel( enei ) );
end

fid = fopen(strcat('Spectrum_',num2str(a)),'wt');
for j = 1:length(ene)
          fprintf(fid,' %g', ene(j));
          fprintf(fid,' %g', psurf(1,j));
          fprintf(fid, '\n');
end
fclose(fid)

multiWaitbar( 'CloseAll' );
plot( eV2nm./enei, psurf ); 
toc
end