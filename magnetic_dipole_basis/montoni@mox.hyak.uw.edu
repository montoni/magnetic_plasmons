%addpath(genpath('/data/home/dmasiell/code/MNPBEM14'));
op2 = bemoptions( 'sim', 'ret', 'interp', 'curv' );

epstab = {epsconst( 1 ),epstable('silver_smallgamma.dat')};

for sep = 1:10;


radius2=30;
r = 2*radius2;
diameter2 = 2*radius2;

ps = trisphere(144,diameter2);

ex = sqrt(3)/2 * r;
ey = (1/2) * r;
ez = sqrt(2/3) * r;

p1 = shift(ps,[-ex-ey*sep,0,0]);
p2 = shift(ps,[-ey*sep,-ey,0]);
p3 = shift(ps,[ey*sep,-ey,0]);
p4 = shift(ps,[ex+ey*sep,0,0]);
p5 = shift(ps,[ey*sep,ey,0]);
p6 = shift(ps,[-ey*sep,ey,0]);
p7 = shift(ps,[-ex/3 - ey*sep,0,ez]);
p8 = shift(ps,[ex/3 + ey*sep,0,ez]);

Ptot= {p1,p2,p3,p4,p5,p6,p7,p8};

p=comparticle(epstab,Ptot,[2,1;2,1;2,1;2,1;2,1;2,1;2,1;2,1],1,2,3,4,5,6,7,8,op2);

units;
enei = eV2nm./linspace( 2.7, 3.7, 101 );
bem = bemsolver( p, op1 );
units
%  plane wave excitation
exc = planewave( [ 0, 1, 0; 1, 0, 0; 0, 1, 0], [ 1, 0, 0; 0, 0, 1; 0, 0, 1], op1 );
%exc2 = planewave( [ 0, 1, 0 ], [ 0, 0, 1], op1 );
%  light wavelength in vacuum

%  allocate scattering and extinction cross sections
sca = zeros( length( enei ), 3 );
ext = zeros( length( enei ), 3 );

for ien = 1 : length( enei )
  %  surface charge
  sig = bem \ exc( p, enei( ien ) );
  %sig2 = bem \ exc2( p, enei(ien));
  %  scattering and extinction cross sections
  
  
  sca( ien, : ) = exc.sca( sig );
  ext( ien, : ) = exc.ext( sig );
  multiWaitbar( 'BEM solver', ien / numel( enei ) );
end
%  close waitbar


ene = linspace(2.7,3.7,101);
fid = fopen('optical_prisms/spectrum_prism_180','wt');
fprintf(fid, ' %s', 'Energy(eV)     Ext/Sca_Zpol     Ext/Sca_Ypol     Ext/Sca_Xpol');
fprintf(fid, '\n');
for i = 1:length(ene)
    fprintf(fid, ' %g', ene(i));
    fprintf(fid, ' %g', ext(i,1));
    fprintf(fid, ' %g', sca(i,1));
    fprintf(fid, ' %g', ext(i,2));
    fprintf(fid, ' %g', sca(i,2));
    fprintf(fid, ' %g', ext(i,3));
    fprintf(fid, ' %g', sca(i,3));
    fprintf(fid,'\n');
end
end