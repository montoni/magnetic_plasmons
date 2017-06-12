%  initialization
clear
close all
clc

%  options for BEM simulation
op = bemoptions( 'sim', 'ret', 'interp', 'curv', 'eels.refine', 2 );

 %  table of dielectric functions
epstab  = { epsconst( 1 ), epsdrude( 'Ag' ), epsconst( 1.5^2 )  };

% nanosphere 
radius1 = 37.5;
radius2 = 1;

%for LoopCount = 1:15
height1 = 220;
r = height1;
diameter1 =2*radius1;
prod = trirod(diameter1,height1,[20 20 20]);
prod = rot(prod,90,[0,1,0]); 
%psphere = trisphere(500, 10);
 %  make particle
p = comparticle( epstab, { prod }, [ 2, 1 ], 1, op );

%plot(p, 'EdgeColor','b') 
%%  EELS excitation     
units;
%  width of electron beam and electron velocity 
[ width, vel ] = deal( 0.2, eelsbase.ene2vel( 200e3 ) );
%  impact parameters
%imp = [ 315, 95 ];
  imp = [ -500 , 0 ];
%  loss energies in eV
enei = eV2nm./linspace( 0.01, 2.5, 60);
ene = transpose(eV2nm./enei);

%%  BEM simulation
%  BEM solver
bem = bemsolver( p, op );
 
%  electron beam excitation
exc = electronbeam( p, imp, width, vel, op );
%  surface and bulk loss
[ psurf, pbulk ] = deal( zeros( size( imp, 1 ), numel( enei ) ) );
 %%
multiWaitbar( 'BEM solver', 0, 'Color', 'g', 'CanCancel', 'on' );
%  loop over energies
tic
for ien = 1 : length( enei )
    tic
  %  surface charges
  sig = bem \ exc( enei( ien ) );
  %  EELS losses
  [ psurf( :,ien ), pbulk( :,ien ) ] = exc.loss( sig );
  eV2nm./enei(ien)
  psurf( :, ien )
  toc
  multiWaitbar( 'BEM solver', ien / numel( enei ) );
end
%  close waitbar
multiWaitbar( 'CloseAll' );
toc
%%
plot( transpose(ene), psurf ); 