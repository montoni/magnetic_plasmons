addpath(genpath('/data/home/dmasiell/code/MNPBEM14'));
clear all
close all
clc
op1 = bemoptions( 'sim', 'ret', 'interp', 'curv' );

 %  table of dielectric functions
epstab  = { epsconst( 1 ), epstable('silver_normal.dat') };
% epstable( 'agsubs.dat' ),
numPart = 3;
theta = linspace(0,2*pi,numPart+1);
phi = 2*pi/(numPart*2);
for diam = 60;
%diam=5;
separation = 1.5*diam;
RAD = separation/(2*sin(phi));

p1 = shift(trisphere(144, diam), [RAD*cos(theta(1)), RAD*sin(theta(1)), 0]);
p2 = shift(trisphere(144, diam), [RAD*cos(theta(2)), RAD*sin(theta(2)), 0]);
p3 = shift(trisphere(144, diam), [RAD*cos(theta(3)), RAD*sin(theta(3)), 0]);
%p4 = shift(trisphere(60, diam), [RAD*cos(theta(4)), RAD*sin(theta(4)), 0]);
%p5 = shift(trisphere(60, diam), [RAD*cos(theta(5)), RAD*sin(theta(5)), 0]);
%p6 = shift(trisphere(60, diam), [RAD*cos(theta(6)), RAD*sin(theta(6)), 0]);

p = comparticle( epstab, { p1, p2, p3 },  ... 
		 [ 2,1;2,1;2,1], 1,2,3, op1 );
plot(p)

units;
enei = eV2nm./linspace( 3.25, 3.75, 26);
ene = transpose(eV2nm./enei);
bem = bemsolver( p, op1 );
 
pt = compoint( p, [ 0, 0, 0 ] );
%  dipole excitation
dip = magdipoleret( pt, [ 0, 1, 0 ], op1 );
%  initialize total and radiative scattering rate
[ sca ] = deal( zeros( numel( enei ), 1 ) );
for ien = 1 : length( enei )
%  surface charge
  sig = bem \ dip( p, enei( ien ) );
  %  total and radiative decay rate
  [ sca( ien ) ] = dip.scattering( sig );
  
  multiWaitbar( 'BEM solver', ien / numel( enei ) );
end
%  close waitbar
multiWaitbar( 'CloseAll' );
end
%%
figure
plot( ene, sca, 'o-' ); 

xlabel( 'Energy (eV)' );
ylabel( 'Scattering' );
% fid = fopen(strcat('hexamer_',num2str(diam/2)),'wt');
% for i = 1:length(ene)
%     fprintf(fid, ' %g', ene(i));
%     fprintf(fid, ' %g', psurf(1,i));
%     fprintf(fid,'\n');
% end



