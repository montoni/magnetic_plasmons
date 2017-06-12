addpath(genpath('/data/home/dmasiell/code/MNPBEM14'));
op2 = bemoptions( 'sim', 'ret', 'interp', 'curv', 'eels.refine', 2  );

epstab = {epsconst( 1 ),epstable('silver_smallgamma.dat')};

for radius2 = [15,20,25,30];

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

[ width, vel ] = deal( 0.01, eelsbase.ene2vel( 200e3 ) );

%imp = [0, ey + 0.1*radius2
imp = [ 1.1*radius2,5*ey];

enei = eV2nm./linspace( 3.0,3.4,201);
ene = transpose(eV2nm./enei);

bem = bemsolver( p, op2 );

exc = electronbeam( p , imp, width, vel, op2 );

% Next Green Functions are calculated.  
Size=75; %number of total points
SizeX=14*radius2; %total dimension, going negative to positve
  SizeY=17*radius2; 

%  points where electric fields are computed
[ x, y ] = meshgrid( (SizeX-1)/2 * linspace( - 1, 1, Size ),(SizeY-1)/2 * linspace( - 1, 1, Size ) ); z = 0 * x; %At this line I multiplied each component by 0 to get the field at 0.

%  place the points into the dielectric media
pt = compoint( p, [ x( : ), y( : ), z( : ) ], 'mindist', 1.5 );

%  Green function between points and particleMo
g = compgreenret( pt, p );

%  loop over wavelengths
for ien = 1 : length( enei )
  %  surface charge
	    sig = bem \ exc( enei( ien ) );
field=g.field(sig);
Bob=pt(real(field.h(:,3)));
 
for k=1:max(length(y));
for i=1:max(length(x));
Bz(i,k)=Bob((k-1)*max(length(x))+i);
    end
 end
    save(strcat('NN_Bfield_',num2str(radius2),'_nm_',num2str(enei(ien))),'Bz','-ascii')


  

  clear Bob field bem2 sig Bz
end

end
