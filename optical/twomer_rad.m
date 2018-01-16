%addpath(genpath('/usr/lusers/montoni/MNPBEM14'));
op2 = bemoptions( 'sim', 'ret', 'interp', 'curv' );

epstab = {epsconst( 1 ),epstable('silver_normal.dat')};
num = [15];
%for sep = 2:12
radius2 = num;
r = 3*radius2;
diameter2 = 2*radius2;

ps = trisphere(144,diameter2);

ex = sqrt(3)/2 * r;
ey = (1/2) * r;

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

Ptot= {p5,p6,p7,p8,p9,p10,p11,p12,p13,p14};

pret=comparticle(epstab,Ptot,[2,1;2,1;2,1;2,1;2,1;2,1;2,1;2,1;2,1;2,1],1,2,3,4,5,6,7,8,9,10,op2);

units;
enei = eV2nm./3.52;
%  set up BEM solver
bem = bemsolver( pret, op2 );
units
%  plane wave excitation
exc = planewave( [ 0, 1, 0], [ 1, 0, 0], op2 );
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

theta2 = reshape( linspace( 0.00001, 2*pi - 0.001 + 0*pi/100, 301 ), [], 1 );
%  directions for emission
dir = [ 0 * theta2, cos( theta2 ),  sin( theta2 ) ];
%  set up spectrum object
spec = spectrum( dir, op2 );

%  farfield radiation
fpl = farfield( spec, sig );
%  norm of Poynting vector
spl = vecnorm( 0.5 * real( cross( fpl.e, conj( fpl.h ), 2 ) ) );
figure()
%  plot electric field
%imagesc( x( : ), z( : ), log10( ee1 ) );  hold on
%plot( [ min( x( : ) ), max( x( : ) ) ], [ 0, 0 ], 'w--' );

%  Cartesian coordinates of Poynting vector
[ sx, sy ] = pol2cart( theta2, 8 * spl / max( spl ) );
%  overlay with Poynting vector
plot( sx, sy, 'b-', 'LineWidth', 1.5 );
xlabel('Y-Direction')
ylabel('Z-Direction')
% fid = fopen(strcat('twomer_optical_sep',num2str(sep-2)),'wt')
% fprintf(fid, ' %s', 'Energy(eV)     Ext_yx     Sca_yx     Ext_xz     Sca_xz');
% fprintf(fid, '\n');
% for i = 1:length(enei)
% 	  fprintf(fid, ' %g', enei(i));
% fprintf(fid, ' %g', ext(i,1));
% fprintf(fid, ' %g', sca(i,1));
% fprintf(fid, ' %g', ext(i,2));
% fprintf(fid, ' %g', sca(i,2));
% fprintf(fid,'\n');
%end
%end
