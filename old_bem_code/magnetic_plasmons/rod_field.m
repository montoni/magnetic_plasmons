%%  initialization
%  options for BEM simulation
clear, clc, clf

%  options for BEM simulation
op = bemoptions( 'sim', 'ret', 'interp', 'curv' );

 %  table of dielectric functions
epstab  = { epsconst( 1 ), epstable( 'silver.dat' ), epsconst( 1.55^2 ),  };

%  X-LENGTH, Y-LENGTH OF RODS
rodlen_left  = [ 60, 215];
rodlen_right = [ 60, 215];
height = 20;
h=0.25;

% HEIGHT OF RODS
rodedge = edgeprofile( height, 10, 'mode', '00', 'min', 0 );

% X-LENGTH, Y-LENGTH OF SUBSTRATE
sublen = [ 900, 700];

% MAKE ROD POLYGON
rod_left  = polygon( 20, 'size', rodlen_left  ); 
rod_right = polygon( 20, 'size', rodlen_right );

% RADIUS OF ONE RING
radius1 = 90;
radius2 = 90;

% LOCATION OF CENTER OF ONE RING
center = 160;

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
dir = [ 0, 0, 1 ]%; 0, 0, 1 ; 0, 1, 0 ; 0, 1, 0 ];
pol = [ 0, 1, 0 ]%; 1, 0, 0 ; 1, 0, 0 ; 1, 0, 0 ];
%  photon wavelengths
units;
%enei = eV2nm./linspace( 1, 2.5, 101 );
%ene = transpose(eV2nm./enei);
%% Next Green Functions are calculated.
%  impact parameters (triangle corner, middle, and midpoint of edge)
[ x, y ] = meshgrid( linspace(  0, 400, 101 ), linspace( -200, 0, 101 ) );
z = 0 * x; %At this line I multiplied each component by 0 to get the field at 0.

%  place the points into the dielectric media
pt = compoint( p, [ x( : ), y( : ), z( : ) ], 'mindist', 0.1 );
%%  BEM solver
%  initialize BEM solver
bem = bemsolver( p, op );
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
    
    Alice  = Alice';
    Alice2 = fliplr(Alice);
    Alice3  = [Alice , Alice2];
    Alice = [flipud(Alice3) ; Alice3];

    
    save(strcat('PW_Efield_',num2str(eV2nm./enei(ien))),'Alice','-ascii')
    
    multiWaitbar( 'BEM solver', ien / numel( enei ) );
    
    clear field bem2 sig2 E Ez Alicex Alicey Alicez
   
  multiWaitbar( 'BEM solver', ien / numel( enei ) );
  toc
end


%  close waitbar
multiWaitbar( 'CloseAll' );


fileID = fopen('PW_sca_cube.txt','w');
fprintf(fileID,'%6s %12s\n','E','sca');
fprintf(fileID,'%6.2f %12.8f\n',sca);
fclose(fileID);
%%
figure(4)
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

for j=1:length(enei)
E=load(strcat('PW_Efield_',num2str(eV2nm./enei(j))),'-ascii'); 
colormap(cmapEfield)
pcolor(E')
axis equal
shading flat
shading interp  
colorbar 
title('E-Field mag')     
title(eV2nm/enei(j));    
caxis([0 40])
drawnow  
%pause
end