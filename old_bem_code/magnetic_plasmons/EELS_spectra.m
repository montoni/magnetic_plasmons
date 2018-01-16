%%  initialization
%  options for BEM simulation
clear, clc, clf

%  options for BEM simulation
op = bemoptions( 'sim', 'ret', 'interp', 'curv' );

 %  table of dielectric functions
epstab  = { epsconst( 1 ), epsdrude( 'Ag' ), epsconst( 1 ),  };
%  X-LENGTH, Y-LENGTH OF RODS
rodlen_left  = [ 10, 20];
rodlen_right = [ 10, 20];
height = 20;
h=0.25;

% HEIGHT OF RODS
rodedge = edgeprofile( height, 10, 'mode', '00', 'min', 0 );

% X-LENGTH, Y-LENGTH OF SUBSTRATE
sublen = [ 900, 700];

% MAKE ROD POLYGON
rod_left  = polygon( 35, 'size', rodlen_left  ); 
rod_right = polygon( 35, 'size', rodlen_right );

% RADIUS OF ONE RING
radius1 = 12.5;
radius2 = 12.5;

% LOCATION OF CENTER OF ONE RING
center = 19 ;

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
p1 = shift(pr, [radius1*cos(theta(1))-center, radius1*sin(theta(1)), h]);
p2 = shift(rot(pr, 120), [radius1*cos(theta(3))-center, radius1*sin(theta(3)), h]);
p3 = shift(rot(pr, 60), [radius1*cos(theta(5))-center, radius1*sin(theta(5)), h]);

p4 = shift(p, [radius2*cos(theta(4))+center, radius2*sin(theta(4)), h]);
p5 = shift(rot(p, 60), [radius2*cos(theta(2))+center, radius2*sin(theta(2)), h]);
p6 = shift(rot(p, 120), [radius2*cos(theta(6))+center, radius2*sin(theta(6)), h]);

% MAKE SUBSTRATE BRICK POLYGON
%subs = polygon( 4, 'size', sublen)
cr  = 15 ;
pcube2 = asymcube2(2*cr, sublen(1), sublen(2), 15, [6, 6, 30, 30, 18 ] );
pcube2 = shift(pcube2,[sublen(1)/2-cr,sublen(1)/2-cr,-height/2]);


 p = comparticle( epstab, {  p1, p2, p3, p4, p5, p6},  ... 
  [ 2,1; 2,1; 2,1; 2,1; 2,1; 2,1 ], 1,2,3,4,5,6,  op );


%plot(p)

%p = comparticle( epstab, {  p1, p2, p3 },  ... 
 %   [ 2,1; 2,1; 2,1 ], 1,2,3,  op );
figure(1)
%plot(p)
%hold on
%scatter([7, 0],[12,15])
%%  EELS excitation     
units;
%  width of electron beam and electron velocity 
[ width, vel ] = deal( 0.2, eelsbase.ene2vel( 200e3 ) );
%  impact parameters
%imp = [ 315, 95 ];
  imp = [ 20, 0 ; 20, 10 ];
%  loss energies in eV
enei = eV2nm./linspace( 1, 3.9, 101);
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
tic
%  loop over energies
%parpool(10)
for ien = 1 : length( enei )
  %  surface charges
  sig = bem \ exc( enei( ien ) );
  %  EELS losses
  [ psurf( :,ien ), pbulk( :,ien ) ] = exc.loss( sig );
  eV2nm./enei(ien);
  psurf( :, ien );
  multiWaitbar( 'BEM solver', ien / numel( enei ) );
end
%  close waitbar
multiWaitbar( 'CloseAll' );
toc
%%
figure(2)
plot( eV2nm./enei, psurf ); 