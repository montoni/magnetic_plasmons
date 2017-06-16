%  initialization
clear, clc
%  options for BEM simulation
op = bemoptions( 'sim', 'ret', 'interp', 'curv', 'eels.refine', 2 );
 %  table of dielectric functions
epstab  = { epsconst( 1 ), epstable( 'gold.dat' ) };

height = 16;
%  X-LENGTH, Y-LENGTH OF RODS
rodlen_left   = [ 10, 20 ];
rodlen_right  = [ 10, 20 ];
% Length of Substrate
sublen = [1000,1000];

% MAKE ROD POLYGON
rod_left   = polygon( 32, 'size', rodlen_left ,'dir', -1 ); 
rod_right  = polygon( 32, 'size', rodlen_right ,'dir', -1 ); 

subPoly   = polygon(  4, 'size', sublen); 

rodedge = edgeprofile( height, 10, 'mode', '00', 'min', 0 );

p1 = rod_left;  p2 = rod_left;  p3 = rod_left;
p4 = rod_right; p5 = rod_right; p6 = rod_right;

center =  7.5;
radius1 = 10.5;
radius2 = radius1;
theta = [0, pi/3, 2*pi/3, pi, 4*pi/3, 5*pi/3];

for j = 1:length(p1.pos)
p1.pos(j,:) = [cos(theta(1)), sin(theta(1)); -sin(theta (1)),cos(theta(1))]*p1.pos(j,:)'...
            - [2*radius1*cos(theta(1))+center, radius1*sin(theta(1))]';
p2.pos(j,:) = [cos(theta(3)), sin(theta(3)); -sin(theta(3)),cos(theta(3))]*p2.pos(j,:)'...
            + [radius1*cos(theta(3))-center, radius1*sin(theta(3))]';
p3.pos(j,:) = [cos(theta(5)), sin(theta(5)); -sin(theta(5)),cos(theta(5))]*p3.pos(j,:)'...
            + [radius1*cos(theta(5))-center, radius1*sin(theta(5))]';

p4.pos(j,:) = [cos(theta(4)), sin(theta(4)); -sin(theta(4)),cos(theta(4))]*p4.pos(j,:)'...
            - [2*radius2*cos(theta(4))-center, radius2*sin(theta(4))]';
p5.pos(j,:) = [cos(theta(2)), sin(theta(2)); -sin(theta(2)),cos(theta(2))]*p5.pos(j,:)'...
            + [radius2*cos(theta(2))+center, radius2*sin(theta(2))]';
p6.pos(j,:) = [cos(theta(6)), sin(theta(6)); -sin(theta(6)),cos(theta(6))]*p6.pos(j,:)'...
            + [radius2*cos(theta(6))+center, radius2*sin(theta(6))]';                
end

poly = [subPoly,p1,p2,p3,p4,p5,p6];
p = tripolygon(poly, rodedge );
p = comparticle( epstab, {p}, [ 2,1 ], 1, op );
clf 
plot(p,'EdgeColor', 'b')
%plot(p)
%%  EELS excitation     
units;
%  width of electron beam and electron velocity 
[ width, vel ] = deal( 0.2, eelsbase.ene2vel( 60e3 ) );
%  impact parameters
%imp = [ 315, 95 ];
imp = [ 0, 0 ; 0, 5 ; 36, 12 ; 100, 100];
%  loss energies in eV
enei = eV2nm./linspace( 1, 3.5, 101);
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
parpool(10)
parfor ien = 1 : length( enei )
  %  surface charges
  sig = bem \ exc( enei( ien ) );
  %  EELS losses
  [ psurf( :,ien ), pbulk( :,ien ) ] = exc.loss( sig );
  psurf( :, ien );
  multiWaitbar( 'BEM solver', ien / numel( enei ) );
end
%  close waitbar
multiWaitbar( 'CloseAll' );
toc
%%
figure(2)
plot( eV2nm./enei, psurf ); 