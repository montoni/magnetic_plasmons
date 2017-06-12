%%  initialization
%  options for BEM simulation
clear, clc, clf

%  options for BEM simulation
op = bemoptions( 'sim', 'ret', 'interp', 'curv' );

 %  table of dielectric functions
epstab  = { epsconst( 1 ), epstable( 'silver.dat' ), epsconst( 2 ),  };
%  X-LENGTH, Y-LENGTH OF RODS
width=30;
scale=2.5;
le=scale*width;
h=0;

% MAKE ROD POLYGON
p  = trirod( width, le, [12,12,12]  ); 
p= rot(p,90,[0,1,0]);
pr = trirod( width, le, [12,12,12]  );
pr= rot(pr,90,[0,1,0]);


% X-LENGTH, Y-LENGTH OF SUBSTRATE
sublen = [ 450, 300];

% RADIUS OF ONE RING
radius1 = 50;
radius2 = 50;

% LOCATION OF CENTER OF ONE RING
center = 75;

% ANGLES
theta = [0, pi/3, 2*pi/3, pi, 4*pi/3, 5*pi/3];
% extrude substrate
%brick = tripolygon( subs, subedge);
% MOVE AND ROTATE RODS
p1 = shift(rot(pr,90)  , [radius1*cos(theta(1))+center, radius1*sin(theta(1)), h]);
p2 = shift(rot(pr, 210), [radius1*cos(theta(3))+center, radius1*sin(theta(3)), h]);
p3 = shift(rot(pr, 150), [radius1*cos(theta(5))+center, radius1*sin(theta(5)), h]);

p4 = shift(rot(p,90), [radius2*cos(theta(4))-center, radius2*sin(theta(4)), h]);
p5 = shift(rot(p, 150), [radius2*cos(theta(2))-center, radius2*sin(theta(2)), h]);
p6 = shift(rot(p, 210), [radius2*cos(theta(6))-center, radius2*sin(theta(6)), h]);

% MAKE SUBSTRATE BRICK POLYGON
%subs = polygon( 4, 'size', sublen)
cr  = 15 ;
a=2.75;
b=1.15;
c=0.9;
%a=2;b=1;c=1;
pcube2 = asymcube2(2*cr, sublen(1), sublen(2), 30, [10, 10, 20, 20, 20 ] );
pcube2 = shift(pcube2,[sublen(1)/a-b*cr,sublen(1)/a-c*cr, -width/2-cr]);


 p = comparticle( epstab, {  p1, p2, p3, p4, p5,p6, pcube2 },  ... 
  [ 2,1; 2,1; 2,1; 2,1; 2,1; 2,1;3,1 ], 1,2,3,4,5,6,7, op );

%plot(p)

%p = comparticle( epstab, {  p1, p2, p3 },  ... 
 %   [ 2,1; 2,1; 2,1 ], 1,2,3,  op );
figure(1)
%plot(p,'EdgeColor', 'b')
%%  EELS excitation     
units;
%  width of electron beam and electron velocity 
[ width, vel ] = deal( 0.1, eelsbase.ene2vel( 200e3 ) );
%  impact parameters
%imp = [ 315, 95 ];
  imp = [ 0, 23 ];
%  loss energies in eV
enei = eV2nm./linspace( 0.8, 3, 91);
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
%  close waitbar
multiWaitbar( 'CloseAll' );
toc
%%
plot( eV2nm./enei, psurf ); 