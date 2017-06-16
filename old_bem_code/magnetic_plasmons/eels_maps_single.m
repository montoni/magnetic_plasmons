%  initialization
clear
close all
clc

%  options for BEM simulation
op = bemoptions( 'sim', 'ret', 'interp', 'curv', 'eels.refine', 2 );

 %  table of dielectric functions
epstab  = { epsconst( 1 ), epsdrude( 'Ag' ), epsconst( 1.5^2 ),  };

%  X-LENGTH, Y-LENGTH OF RODS
rodlen_left  = [ 75, 220];
rodlen_right = [ 75, 190];
height = 20;
h=0.5;

% HEIGHT OF RODS
rodedge = edgeprofile( height, 10, 'mode', '01', 'min', 0 );

% X-LENGTH, Y-LENGTH OF SUBSTRATE
sublen = [ 900, 600];

% HEIGHT OF SUBSTRATE
subedge = edgeprofile( 35, 10, 'mode', '01', 'min', -35 );

% MAKE ROD POLYGON
rod_left = polygon( 30, 'size', rodlen_left ); 
rod_right =polygon( 30, 'size', rodlen_right );

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
  %   [ 2,1 ; 2,1; 2,1; 2,1; 2,1; 2,1; 3,1], 1,2,3,4,5,6,7,  op );
p = comparticle( epstab, { p4},  ... 
[ 2,1 ], 1 , op );


%plot(p)

%p = comparticle( epstab, {  p1, p2, p3, p4, p5, p6 },  ... 
 %   [ 2,1; 2,1; 2,1; 2,1; 2,1; 2,1 ], 1,2,3,4,5,6,  op );

plot(p, 'EdgeColor','b') 
%%  EELS excitation                 
%  width of electron beam and electron velocity 
[ width, vel ] = deal( 0.2, eelsbase.ene2vel( 200e3 ) );
%  impact parameters (triangle corner, middle, and midpoint of edge)
[ x, y ] = meshgrid( linspace(  -400, 400, 201 ), linspace( -200, 0, 101 ) );
%  impact parameters
imp = [ x( : ), y( : ) ];
 
%  loss energies in eV
ene = [ 0.3, 1.57, 1.4 ];%linspace( 2, 3.8, 81 );
 
%  convert energies to nm
units;  enei = eV2nm ./ ene;
%%  BEM simulation
%  BEM solver
tic
bem = bemsolver( p, op );
 
%  electron beam excitation
exc = electronbeam( p, imp, width, vel, op );
%  surface and bulk loss
[ psurf, pbulk ] = deal( zeros( size( imp, 1 ), numel( enei ) ) );
 
multiWaitbar( 'BEM solver', 0, 'Color', 'g', 'CanCancel', 'on' );
%  loop over energies
for ien = 1 : length( ene )
  %  surface charges
  sig = bem \ exc( enei( ien ) );
  %  EELS losses
  [ psurf( :,ien ), pbulk( :,ien ) ] = exc.loss( sig );
  
Alice  = reshape(psurf(:,ien),size(x)) + reshape(pbulk(:,ien),size(x));
Alice  = Alice';
Alice2 = fliplr(Alice);
Alice  = [Alice , Alice2];
 
 save(strcat('Map_single_',num2str(ene(ien))),'Alice','-ascii')
  
  multiWaitbar( 'BEM solver', ien / numel( enei ) );
end
%  close waitbar
multiWaitbar( 'CloseAll' );
toc
%%
fid = fopen('EELS_Z_XY-Z','wt');
for i = 1:length(ene)
    fprintf(fid, ' %g', ene(i));
    fprintf(fid, ' %g', psurf(1,i) + pbulk(1,i));
    fprintf(fid, ' %g', psurf(2,i) + pbulk(2,i));
    fprintf(fid,'\n');
end
%%
cmapEfield = [ 0.0417         0         0
    0.0833         0         0
    0.1250         0         0
    0.1667         0         0
    0.2083         0         0
    0.2500         0         0
    0.2917         0         0
    0.3333         0         0
    0.3750         0         0
    0.4167         0         0
    0.4583         0         0
    0.5000         0         0
    0.5417         0         0
    0.5833         0         0
    0.6250         0         0
    0.6667         0         0
    0.7083         0         0
    0.7500         0         0
    0.7917         0         0
    0.8333         0         0
    0.8750         0         0
    0.9167         0         0
    0.9583         0         0
    1.0000         0         0
    1.0000    0.0417         0
    1.0000    0.0833         0
    1.0000    0.1250         0
    1.0000    0.1667         0
    1.0000    0.2083         0
    1.0000    0.2500         0
    1.0000    0.2917         0
    1.0000    0.3333         0
    1.0000    0.3750         0
    1.0000    0.4167         0
    1.0000    0.4583         0
    1.0000    0.5000         0
    1.0000    0.5417         0
    1.0000    0.5833         0
    1.0000    0.6250         0
    1.0000    0.6667         0
    1.0000    0.7083         0
    1.0000    0.7500         0
    1.0000    0.7917         0
    1.0000    0.8333         0
    1.0000    0.8750         0
    1.0000    0.9167         0
    1.0000    0.9583         0
    1.0000    1.0000         0
    1.0000    1.0000    0.0625
    1.0000    1.0000    0.1250
    1.0000    1.0000    0.1875
    1.0000    1.0000    0.2500
    1.0000    1.0000    0.3125
    1.0000    1.0000    0.3750
    1.0000    1.0000    0.4375
    1.0000    1.0000    0.5000
    1.0000    1.0000    0.5625
    1.0000    1.0000    0.6250
    1.0000    1.0000    0.6875
    1.0000    1.0000    0.7500
    1.0000    1.0000    0.8125
    1.0000    1.0000    0.8750
    1.0000    1.0000    0.9375
    1.0000    1.0000    1.0000];

figure(2)
colormap(cmapEfield)
EELSmap=load('Map_single_1.57'); 
pcolor(EELSmap);
caxis([0 0.01])
shading interp
shading flat