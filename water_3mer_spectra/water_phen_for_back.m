addpath(genpath('/usr/lusers/montoni/MNPBEM14'));
op2 = bemoptions( 'sim', 'ret', 'interp', 'curv' );

epstab = {epsconst( 1.77 ),epstable('silver_normal.dat')};

for radius2 = [1,5,10,15,20,25,30];
num = radius2
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
p8 = shift(ps,[-2*ex,4*ey,0]);
p9 = shift(ps,[-ex,5*ey,0]);
p10 = shift(ps,[0,4*ey,0]);
p12 = shift(ps,[ex,ey,0]);
p13 = shift(ps,[0,2*ey,0]);
p14 = shift(ps,[-ex,ey,0]);

Ptot= {p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p12,p13,p14};

pret=comparticle(epstab,Ptot,[2,1;2,1;2,1;2,1;2,1;2,1;2,1;2,1;2,1;2,1;2,1;2,1;2,1],1,2,3,4,5,6,7,8,9,10,11,12,13,op2);

units;
enei = eV2nm./linspace( 2.8, 3.8, 101 );
%  set up BEM solver
bem = bemsolver( pret, op2 );

units
pol = [0,1,0];
dir = [1,0,0];
%  plane wave excitation
exc = planewave( [ 0, 1, 0], [1, 0, 0], op2);%[ 0, 1, 0; 0, 0, 1], op2 );
%  light wavelength in vacuum
%  allocate scattering and extinction cross sections

sca   = zeros( length( enei ), 1 );
sca_left = zeros( length( enei ), 1 );
sca_right = zeros( length( enei ), 1 );

%  loop over wavelengths
for ien = 1 : length( enei )
  %  surface charge
	    sig = bem \ exc( pret, enei( ien ) );

  %  scattering and extinction cross sections
  [ sca( ien, : ), dsca1 ] = exc.sca( sig );
    
  %The differential scattering function computes scattering on surface
  %elements of a sphere at infinity. We are going to grab the coordinates
  %for the surface elements.
  coords = dsca1.p.nvec; 
  
% Next, let us find the column of coordinates that corresponds to the
  % direction of propogation.
  col = find(dir==1);
    
for j = 1:length(coords)
    % The following if statement determines if the coordinates of the
    % surface elements of the sphere at infinity are positive or negative
    % for the direction of propogation
    if coords(j, col ) >0 
        % If greater that zero, back scattering (check if this is true, or
        % if the labels somehow are flipped)
        sca_left(ien) = sca_left(ien) + dsca1.dsca(j)*dsca1.p.area(j);
    else
        % Otherwise we say it is forward scattering.
        sca_right(ien) = sca_right(ien) + dsca1.dsca(j)*dsca1.p.area(j);
    end
end
  

end

fid = fopen(strcat('water_phen_forback_sca_nm',num2str(num)),'wt')
fprintf(fid, ' %s', 'Energy(eV)     Sca_L     Sca_R');
fprintf(fid, '\n');
for i = 1:length(enei)
	  fprintf(fid, ' %g', enei(i));
fprintf(fid, ' %g', sca_left(i));
fprintf(fid, ' %g', sca_right(i));
fprintf(fid,'\n');
end
end
