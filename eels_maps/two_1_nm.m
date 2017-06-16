addpath(genpath('/usr/lusers/montoni/MNPBEM14'));
op2 = bemoptions( 'sim', 'ret', 'interp', 'curv', 'eels.refine', 2  );

epstab = {epsconst( 1 ),epstable('silver_normal.dat')};
num = 1 %,15,30];
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
enei = eV2nm./linspace( 3.6,3.65,11);
ene = transpose(eV2nm./enei);

%  mesh for electron beams
[ x, y ] = meshgrid( linspace( ex, 4*ex, 31 ), linspace( -3*ey, 3*ey, 61 ) );
%  impact parameters
impact = [ x( : ), y( : ) ];
%  width of electron beam and electron velocity 
[ width, vel ] = deal( 0.2, eelsbase.ene2vel( 200e3 ) );

%  BEM simulation
%  BEM solver
bem = bemsolver( pret, op2 );
%  EELS excitation
exc = electronbeam( pret, impact, width, vel, op2 );
%  electron energy loss probabilities
[ psurf, pbulk ] = deal( zeros( size( impact, 1 ), length( enei ) ) );

multiWaitbar( 'BEM solver', 0, 'Color', 'g', 'CanCancel', 'on' );
%  loop over energies
for ien = 1 : length( enei )
  %  surface charges
	    sig = bem \ exc( enei( ien ) );
  %  EELS losses
  [ psurf( :, ien ), pbulk( :, ien ) ] = exc.loss( sig );
  
fid = fopen(strcat('1_nm_EELS_map_',num2str(ene(ien))),'wt');
prob = reshape( psurf( :, ien ) + pbulk( :, ien ), size( x ) );
probnew = [fliplr( prob( 2 : end, : ) ); prob];
for numcol = 1:length(probnew(1,:))
	       for numrow = 1:length(probnew(:,1))
			      fprintf(fid, ' %g', probnew(numrow,numcol));
    end
    fprintf(fid,'\n');
end
  clear prob probnew
end

