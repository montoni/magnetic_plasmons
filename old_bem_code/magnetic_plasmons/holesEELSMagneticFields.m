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
sublen = [300,300];

% MAKE ROD POLYGON
rod_left   = polygon( 30, 'size', rodlen_left ,'dir', -1 ); 
rod_right  = polygon( 30, 'size', rodlen_right ,'dir', -1 ); 

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
units;
%  width of electron beam and electron velocity 
[ width, vel ] = deal( 0.2, eelsbase.ene2vel( 200e3 ) );
%  impact parameters
%imp = [ 315, 95 ];
  imp = [ 36, 12 ]; 
%  loss energies in eV
ene = [ 1.9, 2.17 ];
enei = eV2nm./ene;

%%
%  impact parameters (triangle corner, middle, and midpoint of edge)
[ x, y ] = meshgrid( linspace(  -35, 35, 101 ), linspace( -35, 35, 101 ) );
z = 0 * x; %At this line I multiplied each component by 0 to get the field at 0.

%  place the points into the dielectric media
pt = compoint( p, [ x( : ), y( : ), z( : ) ], 'mindist', 0.1 );
%%  BEM simulation
%  BEM solver
bem = bemsolver( p, op );
g = compgreenret( pt, p );
%  electron beam excitation
exc = electronbeam( p, imp, width, vel, op );
%  surface and bulk loss
[ psurf, pbulk ] = deal( zeros( size( imp, 1 ), numel( enei ) ) );
%%
multiWaitbar( 'BEM solver', 0, 'Color', 'g', 'CanCancel', 'on' );
%  loop over wavelengths
for ien = 1 : length( enei )
    clear Alice
    tic
  %  surface charges
  sig = bem \ exc( enei( ien ) );
  
    field=g.field(sig);
    Alicex=pt(real(field.e(:,1)));
    Alicey=pt(real(field.e(:,2)));
    Alicez=pt(real(field.e(:,3)));
    Alice=sqrt(Alicex.^2+Alicey.^2+Alicez.^2);
    
    Alice = reshape(Alice,size(x));
    
    %Alice  = Alice';
    %Alice2 = fliplr(Alice);
    %Alice  = [Alice , Alice2];
  %  Alice = [flipud(Alice3) ; Alice3];
    
    Bob=pt(real(field.h(:,3)));
    
    Bob = reshape(Bob,size(x));
    
    %Bob  = Bob';
    %Bob2 = fliplr( Bob);
    %Bob  = [ Bob ,  Bob2];
  %  Bob = [flipud( Bob3) ;  Bob3];

    
    save(strcat('Efield_holes_EELS_',num2str(eV2nm./enei(ien))),'Alice','-ascii')
    save(strcat('Bfield_holes_EELS_',num2str(eV2nm./enei(ien))),'Bob','-ascii')
    
    multiWaitbar( 'BEM solver', ien / numel( enei ) );
    
    clear field bem2 sig2 E Ez Alicex Alicey Alicez
   
  multiWaitbar( 'BEM solver', ien / numel( enei ) );
  toc
end

%  close waitbar
multiWaitbar( 'CloseAll' );