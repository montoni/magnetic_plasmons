%  initialization
clear
close all
clc

%  options for BEM simulation
op = bemoptions( 'sim', 'ret', 'interp', 'curv', 'eels.refine', 2 );
 
 %  table of dielectric functions
epstab  = { epsconst( 1 ), epstable( 'Au_tn.dat' ), epsconst( 1.5^2 )  };

%  X-LENGTH, Y-LENGTH OF RODS
rodlen_left  = [ 25, 50];
rodlen_right = [ 25, 50];
height = 15;
h=0.5;

% HEIGHT OF RODS
rodedge = edgeprofile( height, 10, 'mode', '01', 'min', 0 );

% X-LENGTH, Y-LENGTH OF SUBSTRATE
sublen = [ 900, 600];

% HEIGHT OF SUBSTRATE
subedge = edgeprofile( 35, 10, 'mode', '01', 'min', -35 );

% MAKE ROD POLYGON
rod_left = polygon( 20, 'size', rodlen_left ); 
rod_right =polygon( 20, 'size', rodlen_right );

% RADIUS OF ONE RING
radius1 = 100;
radius2 = 120;

% LOCATION OF CENTER OF ONE RING
center = 200;

% MAKE SUBSTRATE BRICK POLYGON
%subs = polygon( 4, 'size', sublen)
cr  = 15 ;
pcube2 = asymcube2(2*cr, sublen(1), sublen(2), 35, [4, 4, 40, 40, 20 ] );
pcube2 = shift(pcube2,[sublen(1)/2,sublen(1)/2,-height/2]);

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
p4 = shift(p, [radius2*cos(theta(4))-center, radius1*sin(theta(4)), h]);
p5 = shift(rot(p, 60), [radius2*cos(theta(2))-center, radius1*sin(theta(2)), h]);
p6 = shift(rot(p, 120), [radius2*cos(theta(6))-center, radius1*sin(theta(6)), h]);

%  make particle pup1, plo1, pup2, plo2, 2,1; 2, 3; 2,1; 2,3; [3,4],
% p = comparticle( epstab, {  p1, p2, p3, p4, p5, p6, pcube2},... 
%     [ 2,1 ; 2,1; 2,1; 2,1; 2,1; 2,1; 3,1], 1,2,3,4,5,6,7,  op );
p = comparticle( epstab, { p }, [ 2,1 ], 1 , op );


plot(p, 'EdgeColor','b') 
%%  EELS excitation     
units;
%  width of electron beam and electron velocity 
[ width, vel ] = deal( 0.2, eelsbase.ene2vel( 200e3 ) );
%  impact parameters
%imp = [ 315, 95 ];
  imp = [ 0 , -60 ];
%  loss energies in eV
enei = eV2nm./linspace( 1.25, 5.0, 200);
ene = transpose(eV2nm./enei);

%%  BEM simulation
%  BEM solver
bem = bemsolver( p, op );
 
%  electron beam excitation
exc = electronbeam( p, imp, width, vel, op );
%  surface and bulk loss
[ psurf, pbulk ] = deal( zeros( size( imp, 1 ), numel( enei ) ) );
 %%
multiWaitbar( 'BEM solver', 0, 'Color', 'g', 'CanCancel', 'on' );
%  loop over energies
tic
for ien = 1 : length( enei )
  %  surface charges
  sig = bem \ exc( enei( ien ) );
  %  EELS losses
  [ psurf( :,ien ), pbulk( :,ien ) ] = exc.loss( sig );
  multiWaitbar( 'BEM solver', ien / numel( enei ) );
end
%  close waitbar
multiWaitbar( 'CloseAll' );
toc
%%
plot( transpose(ene), psurf );