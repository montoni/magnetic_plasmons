addpath(genpath('/usr/lusers/montoni/MNPBEM14'));
op2 = bemoptions( 'sim', 'ret', 'interp', 'curv', 'eels.refine', 2  );

epstab = {epsconst( 1 ),epstable('silver_normal.dat')};
for num = [1,5,10,15,20,25,30];
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

[ width, vel ] = deal( 0.01, eelsbase.ene2vel( 200e3 ) );

imp = [ex, ey+1.1*radius2; ex-1.1*radius2, ey];

enei = eV2nm./linspace( 3.0,3.65,401);
ene = transpose(eV2nm./enei);

bemr = bemsolver( pret, op2 );

excr = electronbeam( pret , imp, width, vel, op2 );

[ psurfret , pbulkret  ] = deal( zeros( size( imp, 1 ), numel( enei ) ) );

for ien = 1 : length( enei )
  sigr = bemr \ excr( enei( ien ) );
  [ psurfret( :,ien ), pbulkret( :,ien ) ] = excr.loss( sigr );
end

fid = fopen([strcat('nm',num2str(num),'_twomer')],'wt');
for j = 1:length(ene)
          fprintf(fid,' %g', ene(j));
          fprintf(fid,' %g', psurfret(1,j));
          fprintf(fid,' %g', psurfret(2,j));
          fprintf(fid, '\n');
end
fclose(fid)
end
