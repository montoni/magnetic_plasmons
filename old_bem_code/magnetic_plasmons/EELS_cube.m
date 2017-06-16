%%  initialization
%  options for BEM simulation
clear, clc, clf

%  options for BEM simulation
op1 = bemoptions( 'sim', 'ret', 'interp', 'curv' );

 %  table of dielectric functions
epstab  = { epsconst( 1 ), epstable( 'silver.dat' ), epsconst( 3 )  };
diam = 20;
h=0.1;
th = 3;


% X-LENGTH, Y-LENGTH OF SUBSTRATE
sublen = [ 100, 100 ];


% MAKE SUBSTRATE BRICK POLYGON
%subs = polygon( 4, 'size', sublen)
pcube2 = asymcube2(th, sublen(1), sublen(2), th, [10, 10, 55, 55, 20 ] );
pcube2 = shift(pcube2,[sublen(1)/2-th/2,sublen(1)/2-th/2,-th/2]);

psphere = shift(trisphere(1000, diam), [0, 0, diam/2+h]);

%p = comparticle( epstab, { psphere }, [ 2,1 ], 1, op1 );

 p = comparticle( epstab, { psphere, pcube2 },  ... 
  [ 2,1 ; 3,1], 1,2, op1 );

%plot(p)

%p = comparticle( epstab, {  p1, p2, p3 },  ... 
 %   [ 2,1; 2,1; 2,1 ], 1,2,3,  op );
figure(1)
plot(p,'EdgeColor','b')
%plot(p)
%%  EELS excitation     
units;
%  width of electron beam and electron velocity 
[ width, vel ] = deal( 0.01, eelsbase.ene2vel( 200e3 ) );
%  impact parameters
%imp = [ 315, 95 ];
  imp = [ 0, diam+2 ];
%  loss energies in eV
enei = eV2nm./linspace( 3.2, 4.1, 51);
ene = transpose(eV2nm./enei);
%%  BEM simulation
%  BEM solver
bem = bemsolver( p, op1 );

%  electron beam excitation
exc = electronbeam( p, imp, width, vel, op1 );
%  surface and bulk loss
[ psurf, pbulk ] = deal( zeros( size( imp, 1 ), numel( enei ) ) );
 %%
multiWaitbar( 'BEM solver', 0, 'Color', 'g', 'CanCancel', 'on' );
%  loop over energies
tic 
parpool(11)
parfor ien = 1 : length( enei )
  %tic
  %  surface charges
  sig = bem \ exc( enei( ien ) );
  %  EELS losses
  [ psurf( :,ien ), pbulk( :,ien ) ] = exc.loss( sig );
  %eV2nm./enei(ien)
  %psurf( :, ien )
  %toc
  multiWaitbar( 'BEM solver', ien / numel( enei ) );
end
%  close waitbar
multiWaitbar( 'CloseAll' );
toc
%%
plot( eV2nm./enei, psurf ); 