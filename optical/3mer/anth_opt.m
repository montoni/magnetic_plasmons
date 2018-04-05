addpath(genpath('/usr/lusers/montoni/MNPBEM14'));
op2 = bemoptions( 'sim', 'ret', 'interp', 'curv' );

epstab = {epsconst( 1 ),epstable('silver_normal.dat')};

for num = [1,5,10,15,20,25,30];
radius2 = num;
val = 5;

r = 3*radius2;
diameter2 = 2*radius2;

ps = trisphere(144,diameter2);

ex = sqrt(3)/2 * r;
ey = (1/2) * r;

p1 = shift(ps,[-2*ex,2*ey,0]);
p2 = shift(ps,[-3*ex,ey,0]);
p3 = shift(ps,[-3*ex,-ey,0]);
p4 = shift(ps,[-2*ex,-2*ey,0]);
p5 = shift(ps,[-ex,-ey,0]);
p6 = shift(ps,[0,-2*ey,0]);
p7 = shift(ps,[ex,-ey,0]);
p8 = shift(ps,[2*ex,-2*ey,0]);
p9 = shift(ps,[3*ex,-ey,0]);
p10 = shift(ps,[3*ex,ey,0]);
p11 = shift(ps,[2*ex,2*ey,0]);
p12 = shift(ps,[ex,ey,0]);
p13 = shift(ps,[0,2*ey,0]);
p14 = shift(ps,[-ex,ey,0]);

Ptot= {p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14};

pret=comparticle(epstab,Ptot,[2,1;2,1;2,1;2,1;2,1;2,1;2,1;2,1;2,1;2,1;2,1;2,1;2,1;2,1],1,2,3,4,5,6,7,8,9,10,11,12,13,14,op2);

units;
enei = eV2nm./linspace( 2.8, 3.8, 101 );
%  set up BEM solver
bem = bemsolver( pret, op2 );
units
%  plane wave excitation
exc = planewave( [ 1, 0, 0; 0, 1, 0], [ 0, 1, 0; 0, 0, 1], op2 );
%  light wavelength in vacuum
%  allocate scattering and extinction cross sections
sca = zeros( length( enei ), 2 );
ext = zeros( length( enei ), 2 );
%  loop over wavelengths

for ien = 1 : length( enei )
  %  surface charge
            sig = bem \ exc( pret, enei( ien ) );
%sig2 = bem \ exc2( p, enei(ien));
  %  scattering and extinction cross section
  sca( ien, : ) = exc.sca( sig );
ext( ien, : ) = exc.ext( sig );
end
fid = fopen(strcat('anth_optical_nm',num2str(num)),'wt')
  fprintf(fid, ' %s', 'Energy(eV)     Ext_Z     Sca_Z     Ext_X     Sca_X');
fprintf(fid, '\n');
for i = 1:length(enei)
          fprintf(fid, ' %g', enei(i));
fprintf(fid, ' %g', ext(i,1));
fprintf(fid, ' %g', sca(i,1));
fprintf(fid, ' %g', ext(i,2));
fprintf(fid, ' %g', sca(i,2));
fprintf(fid,'\n');
end
end
