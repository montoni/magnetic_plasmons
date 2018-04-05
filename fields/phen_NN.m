addpath(genpath('/usr/lusers/montoni/MNPBEM14'));
op2 = bemoptions( 'sim', 'ret', 'interp', 'curv' );

epstab = {epsconst( 1 ),epstable('silver_022717.dat')};

for radius2 = [1,5,10,15,20,25,30];

val = 5;

r = 3*radius2;
diameter2 = 2*radius2;

ps = trisphere(144,diameter2);

ex = sqrt(3)/2 * r;
ey = (1/2) * r;

p1 = shift(ps,[-ex,2*ey,0]);
p2 = shift(ps,[-2*ex,ey,0]);
p3 = shift(ps,[-2*ex,-ey,0]);
p4 = shift(ps,[-ex,-2*ey,0]);
p5 = shift(ps,[0,-ey,0]);
p6 = shift(ps,[ex,-2*ey,0]);
p7 = shift(ps,[2*ex,-ey,0]);
p8 = shift(ps,[-ex,4*ey,0]);
p9 = shift(ps,[0,5*ey,0]);
p10 = shift(ps,[ex,4*ey,0]);
p12 = shift(ps,[2*ex,ey,0]);
p13 = shift(ps,[ex,2*ey,0]);
p14 = shift(ps,[0,ey,0]);

Ptot= {p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p12,p13,p14};

p=comparticle(epstab,Ptot,[2,1;2,1;2,1;2,1;2,1;2,1;2,1;2,1;2,1;2,1;2,1;2,1;2,1],1,2,3,4,5,6,7,8,9,10,11,12,13,op2);

units;

exc = planewave( [ 0, 1, 0], [ 1, 0, 0], op2 );
%  light wavelength in vacuum
eV = linspace(3.0,3.8,41);
enei = eV2nm./linspace( 3.0,3.8,41);
ene = transpose(eV2nm./enei);

%  allocate scattering and extinction cross sections
sca = zeros( length( enei ), 2 );
ext = zeros( length( enei ), 2 );


bem = bemsolver( p, op2 );

% Next Green Functions are calculated.  
Size=75; %number of total points
SizeX=14*radius2; %total dimension, going negative to positve
  SizeY=17*radius2; 

%  points where electric fields are computed
[ x, y ] = meshgrid( (SizeX-1)/2 * linspace( - 1, 1, Size ),(SizeY-1)/2 * linspace( - 1, 1, Size ) ); z = 0 * x; %At this line I multiplied each component by 0 to get the field at 0.

%  place the points into the dielectric media
pt = compoint( p, [ x( : ), y( : ), z( : ) ], 'mindist', 0.1*radius2 );

%  Green function between points and particleMo
g = compgreenret( pt, p );

multiWaitbar( 'BEM solver', 0, 'Color', 'g', 'CanCancel', 'on' );
%  loop over wavelengths
for ien = 1 : length( enei )
  %  surface charge
	    sig = bem \ exc(p, enei( ien ) );
sca( ien, : ) = exc.sca( sig );
ext( ien, : ) = exc.ext( sig );
field=g.field(sig);
Bob=pt(real(field.h(:,3)));
 
for k=1:max(length(y));
for i=1:max(length(x));
Bz(i,k)=Bob((k-1)*max(length(x))+i);
    end
 end
    save(strcat('NN_phenBfield_',num2str(radius2),'_nm_',num2str(eV(ien))),'Bz','-ascii')

  clear Bob field sig Bz
end

      save(strcat('NN_phensca_',num2str(radius2),'_nm'),'sca','-ascii')

end
