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
rod_left   = polygon( 35, 'size', rodlen_left ,'dir', -1 ); 
rod_right  = polygon( 35, 'size', rodlen_right ,'dir', -1 ); 

subPoly   = polygon(  4, 'size', sublen); 

rodedge = edgeprofile( height, 10, 'mode', '00', 'min', 0 );

p1 = rod_left  ;  p2 = rod_left  ;  p3 = rod_left;
p4 = rod_right ;  p5 = rod_right ;  p6 = rod_right;

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
%plot(p,'EdgeColor', 'b')
plot(p)
%%  EELS excitation                 
%  width of electron beam and electron velocity 
[ width, vel ] = deal( 0.1, eelsbase.ene2vel( 200e3 ) );
%  impact parameters (triangle corner, middle, and midpoint of edge)
[ x, y ] = meshgrid( linspace(  -70, 70, 101 ), linspace( 70, 70, 101 ) );
%  impact parameters
imp = [ x( : ), y( : ) ];
 
%  loss energies in eV
%ene = [1.25, 1.45 , 1.54 ];%linspace( 2, 3.8, 81 );
ene = [ 1.913, 2.375, 2.65, 3.25, 3.58 ] ;  
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
%Alice  = Alice';
%Alice2 = fliplr(Alice);
%Alice  = [Alice , Alice2];
%Alice = [flipud(Alice3) ; Alice3];

 
 save(strcat('ribbed_torus',num2str(ene(ien))),'Alice','-ascii')
  
  multiWaitbar( 'BEM solver', ien / numel( enei ) );
end
%  close waitbar
multiWaitbar( 'CloseAll' );
toc
%%
cmapEfield = [
    0.0417         0         0
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

figure(3)
colormap(cmapEfield)
EELSmap=load('Map_holes_10X20_1micronX1micron1.975'); 
pcolor(EELSmap);
caxis([0.015 0.0185])
colorbar
shading interp
shading flat