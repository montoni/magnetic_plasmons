%%  initialization
%  options for BEM simulation
clear, clc

%  options for BEM simulation
op1 = bemoptions( 'sim', 'ret', 'interp', 'curv' );


 %  table of dielectric functions
epstab  = { epsconst( 1 ), epsdrude( 'Ag' ), epsconst( 1.5^2 )  };

% nanosphere 
radius1 = 37.5;
radius2 = 1;

%for LoopCount = 1:15
height1 = 0;
r = height1;
diameter1 =2*radius1;
prod = trirod(diameter1,height1,[15 15 15]);
prod = rot(prod,90,[0,1,0]); 
%psphere = trisphere();
 %  make particle
p = comparticle( epstab, { prod }, [ 2, 1 ], 1, op1 );

plot(p, 'EdgeColor','b') 
%% 
%  TM mode, excitation from above
dir = [ 0, 0, 1 ];
pol = [ 1, 0, 0 ];
%  photon wavelengths
units;
enei = eV2nm./linspace( 0.1, 2.5, 90 );
ene = transpose(eV2nm./enei);
%%  BEM solver
%  initialize BEM solver
bem = bemsolver( p, op1 );
%  initialize plane wave excitation
exc = planewave( pol, dir, op1 );
%  scattering cross section
sca = zeros( numel( enei ), size( dir, 1 ) );
ext = zeros( numel( enei ), size( dir, 1 ) );
%%
multiWaitbar( 'BEM solver', 0, 'Color', 'g', 'CanCancel', 'on' );
%  loop over wavelengths
for ien = 1 : length( enei )
  %  surface charges
  sig = bem \ exc( p, enei( ien ) );
  multiWaitbar( 'BEM solver', ien / numel( enei ) );
  
  %  scattering cross section
  ext( ien, : ) = exc.ext( sig );
  sca( ien, : ) = exc.sca( sig );
 
  multiWaitbar( 'BEM solver', ien / numel( enei ) );
end


%  close waitbar
multiWaitbar( 'CloseAll' );

%%  final plot
figure(3)
plot( ene, sca(:,1) );

xlabel( 'Energy (eV)' );
ylabel( 'Scattering cross section (nm^2)' );