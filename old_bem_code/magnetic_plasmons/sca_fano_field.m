%%  initialization
%  options for BEM simulation
clear, clc, clf

%  options for BEM simulation
op = bemoptions( 'sim', 'ret', 'interp', 'curv' );

 %  table of dielectric functions
epstab  = { epsconst( 1 ), epstable( 'silver.dat' ), epsconst( 1.55^2 ),  };

%  X-LENGTH, Y-LENGTH OF RODS
rodlen_left  = [ 70, 225];
rodlen_right = [ 70, 190];
height = 20;
h=0.25;

% HEIGHT OF RODS
rodedge = edgeprofile( height, 10, 'mode', '00', 'min', 0 );

% X-LENGTH, Y-LENGTH OF SUBSTRATE
sublen = [ 900, 700];

% MAKE ROD POLYGON
rod_left  = polygon( 40, 'size', rodlen_left  ); 
rod_right = polygon( 40, 'size', rodlen_right );

% RADIUS OF ONE RING
radius1 = 100;
radius2 = 110;

% LOCATION OF CENTER OF ONE RING
center = 170;

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

% MAKE SUBSTRATE BRICK POLYGON
%subs = polygon( 4, 'size', sublen)
cr  = 15 ;
pcube2 = asymcube2(2*cr, sublen(1), sublen(2), 15, [6, 6, 30, 30, 18 ] );
pcube2 = shift(pcube2,[sublen(1)/2-cr,sublen(1)/2-cr,-height/2]);


p = comparticle( epstab, {  p1, p2, p3, p4, p5, p6},  ... 
 [ 2,1; 2,1; 2,1; 2,1; 2,1; 2,1 ], 1,2,3,4,5,6,  op );


%plot(p)

%p = comparticle( epstab, {  p1, p2, p3, p4, p5, p6 },  ... 
 %   [ 2,1; 2,1; 2,1; 2,1; 2,1; 2,1 ], 1,2,3,4,5,6,  op );
figure(1)
plot(p, 'EdgeColor','b')
%% 
%  TM mode, excitation from above
dir = [ 1, 0, 0 ];%; 0, 0, 1 ; 0, 1, 0 ; 0, 1, 0 ];
pol = [ 0, 1, 0 ];%; 1, 0, 0 ; 1, 0, 0 ; 1, 0, 0 ];
%  photon wavelengths
ene = [ 1.375, 1.525 ] ;  
%  convert energies to nm
units;  enei = eV2nm ./ ene;
%% Next Green Functions are calculated.
%  impact parameters (triangle corner, middle, and midpoint of edge)
[ x, y ] = meshgrid( linspace(  -400, 400, 201 ), linspace( -200, 200, 101 ) );
z = 0 * x; %At this line I multiplied each component by 0 to get the field at 0.

%  place the points into the dielectric media
pt = compoint( p, [ x( : ), y( : ), z( : ) ], 'mindist', 0.1 );
%%  BEM solver
%  initialize BEM solver
bem = bemsolver( p, op );
g = compgreenret( pt, p );
%  initialize plane wave excitation
exc = planewave( pol, dir, op );
%  scattering cross section
sca = zeros( numel( enei ), size( dir, 1 ) );
ext = zeros( numel( enei ), size( dir, 1 ) );
%%
multiWaitbar( 'BEM solver', 0, 'Color', 'g', 'CanCancel', 'on' );
%  loop over wavelengths
for ien = 1 : length( enei )
    clear Alice
    tic
  %  surface charges
  sig = bem \ exc( p, enei( ien ) );
  
    field=g.field(sig);
    Alicex=pt(real(field.e(:,1)));
    Alicey=pt(real(field.e(:,2)));
    Alicez=pt(real(field.e(:,3)));
    Alice=sqrt(Alicex.^2+Alicey.^2+Alicez.^2);
    
    Alice = reshape(Alice,size(x));
    
    %Alice  = Alice';
    %Alice2 = fliplr(Alice);
    %Alice3  = [Alice , Alice2];

    Bob=pt(real(field.h(:,3)));
    Bob = reshape(Bob,size(x));
    
    %Bob  = Bob';
    %Bob2 = fliplr( Bob );
    %Bob  = [ Bob ,  Bob2];

    save(strcat('PW_Efield_fano_kxPy_',num2str(eV2nm./enei(ien))),'Alice','-ascii')
    save(strcat('PW_Bfield_fano_kxPy_',num2str(eV2nm./enei(ien))),'Bob','-ascii')

    multiWaitbar( 'BEM solver', ien / numel( enei ) );
    
    clear field bem2 sig2 E Ez Alicex Alicey Alicez
   
  multiWaitbar( 'BEM solver', ien / numel( enei ) );
  toc
end
%  close waitbar
multiWaitbar( 'CloseAll' );