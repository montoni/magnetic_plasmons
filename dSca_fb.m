%%  initialization

%  options for BEM simulation
clear, clc
op = bemoptions( 'sim', 'ret', 'interp', 'curv' );
epstab  = { epsconst( 1^2 ), epstable( 'silver.dat' )};

diameter = 60;
sep = 3;
radius = (diameter+sep)*cos(pi/6)*2/3;
pr = trisphere(128,diameter);

theta = [0, pi/3, 2*pi/3, pi, 4*pi/3, 5*pi/3];

 p1 = shift(pr, [radius*cos(theta(1)), radius*sin(theta(1)), 0]);
 p2 = shift(pr, [radius*cos(theta(3)), radius*sin(theta(3)), 0]);
 p3 = shift(pr, [radius*cos(theta(5)), radius*sin(theta(5)), 0]);

ptot = {p1,p2,p3};
p = comparticle( epstab, ptot, [ 2,1; 2,1; 2,1 ], 1,2,3, op );

clf
plot(p)

%%  BEM simulation
%  set up BEM solver
bem = bemsolver( p, op );
units;

pol = [ 1, 0, 0 ];
dir = [ 0, 1, 0 ];
%  plane wave excitation
exc = planewave( pol, dir, op );
%  light wavelength in vacuum
enei = eV2nm./linspace( 2.2, 3.7, 101 );

%  allocate scattering and extinction cross sections
sca   = zeros( length( enei ), 1 );
sca_f = zeros( length( enei ), 1 );
sca_b = zeros( length( enei ), 1 );

%%
multiWaitbar( 'BEM solver', 0, 'Color', 'g', 'CanCancel', 'on' );
%  loop over wavelengths
for ien = 1 : length( enei )
  %  surface charge
  sig = bem \ exc( p, enei( ien ) );

  %  scattering and extinction cross sections
  [ sca( ien, : ), dsca1 ] = exc.sca( sig );
    
  %The differential scattering function computes scattering on surface
  %elements of a sphere at infinity. We are going to grab the coordinates
  %for the surface elements.
  coords = dsca1.p.nvec; 
  
  % Next, let's find the column of coordinates that corresponds to the
  % direction of propogation.
  col = find(dir==1);
    
for j = 1:length(coords)
    % The following if statement determines if the coordinates of the
    % surface elements of the sphere at infinity are positive or negative
    % for the direction of propogation
    if coords(j, col ) >0 
        % If greater that zero, back scattering (check if this is true, or
        % if the labels somehow are flipped)
        sca_b(ien) = sca_b(ien) + dsca1.dsca(j)*dsca1.p.area(j);
    else
        % Otherwise we say it is forward scattering.
        sca_f(ien) = sca_f(ien) + dsca1.dsca(j)*dsca1.p.area(j);
    end
end
  
  multiWaitbar( 'BEM solver', ien / numel( enei ) );
end
%  close waitbar
multiWaitbar( 'CloseAll' );

%%
ene = eV2nm./enei;

figure(2)
clf(2)
plot( ene, sca,'k',ene,sca_f,'b',ene,sca_b,'r','LineWidth',2); 

xlabel( 'Wavelength (nm)' );
ylabel( 'Scattering cross section (nm^2)' );

