addpath(genpath('/usr/lusers/montoni/MNPBEM14'));
clear all
close all
clc
op2 = bemoptions( 'sim', 'ret', 'interp', 'curv' );

epstab = {epsconst( 1 ),epstable('silver_normal.dat')};
for num = [1,5,10,15,20,25,30];
sep = [3];
radius2 = num;
r = sep*radius2;
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

imp = [ex, ey+1.1*radius2; ex+1.1*radius2, ey; 3*ex+1.1*radius2, 0];

enei = eV2nm./linspace( 3,4,101);
ene = linspace(3,4,101);

bemr = bemsolver( pret, op2 );

excr = eelsret( pret , imp, width, vel, op2 );

[ psurfret, pbulkret ] = deal( zeros( size( imp, 1 ), numel( enei ) ) );
cathodo   = zeros( size(imp,1),numel( enei ));
NS_cat_left = zeros( 1,numel( enei ) );
NS_cat_right = zeros( 1,numel( enei ) );
NN_cat_left = zeros( 1,numel( enei ) );
NN_cat_right = zeros( 1,numel( enei ) );
both_cat_left = zeros( 1,numel( enei ) );
both_cat_right = zeros( 1,numel( enei ) );

dir = [1,0,0];
pol = [0,1,0];
multiWaitbar( 'BEM solver', 0, 'Color', 'g', 'CanCancel', 'on' );
for ien = 1 : length( enei )
	    sigr = bemr \ excr( enei( ien ) );
[ psurfret( :,ien ), pbulkret( :,ien ) ] = excr.loss( sigr );


% cathodoluminescence
[ cathodo( :, ien ), dcat1 ] = excr.rad( sigr );
    
  %The differential scattering function computes scattering on surface
  %elements of a sphere at infinity. We are going to grab the coordinates
  %for the surface elements.
  coords = dcat1.p.nvec; 
  
% Next, let us find the column of coordinates that corresponds to the
  % direction of propogation.
  col = find(pol ==0 & dir == 0);
    
for j = 1:length(coords)
    % The following if statement determines if the coordinates of the
    % surface elements of the sphere at infinity are positive or negative
    % for the direction of propogation
	  if coords(j, col ) > .9 
        % If greater that zero, back scattering (check if this is true, or
        % if the labels somehow are flipped)
        NS_cat_left(ien) = NS_cat_left(ien) + dcat1.dprad.struct.dsca(j)*dcat1.p.area(j);
        NN_cat_left(ien) = NN_cat_left(ien) + dcat1.dprad.struct.dsca(j+length(coords))*dcat1.p.area(j);
        both_cat_left(ien) = both_cat_left(ien) + dcat1.dprad.struct.dsca(j+length(coords)+length(coords))*dcat1.p.area(j);
    else
        % Otherwise we say it is forward scattering.
        NS_cat_right(ien) = NS_cat_right(ien) + dcat1.dprad.struct.dsca(j)*dcat1.p.area(j);
        NN_cat_right(ien) = NN_cat_right(ien) + dcat1.dprad.struct.dsca(j+length(coords))*dcat1.p.area(j);
        both_cat_right(ien) = both_cat_right(ien) + dcat1.dprad.struct.dsca(j+length(coords)+length(coords))*dcat1.p.area(j);
      end
    
end
  
multiWaitbar( 'BEM solver', ien / numel( enei ) );
end
multiWaitbar( 'CloseAll' );

fid = fopen([strcat('nm',num2str(num),'_2mer_CL_updo_broad')],'wt');
fprintf(fid, ' %s', 'Energy(eV)     NS_eels     NN_eels     NS_cat     NN_cat     NS_CL_left     NS_CL_right     NN_CL_left     NN_CL_right     both_left     both_right');
fprintf(fid, '\n');
for j = 1:length(ene)
          fprintf(fid,' %g', ene(j));
fprintf(fid,' %g', psurfret(1,j));
fprintf(fid,' %g', psurfret(2,j));
fprintf(fid,' %g', cathodo(1,j));
fprintf(fid,' %g', cathodo(2,j));
fprintf(fid,' %g', NS_cat_left(1,j));
fprintf(fid,' %g', NS_cat_right(1,j));
fprintf(fid,' %g', NN_cat_left(1,j));
fprintf(fid,' %g', NN_cat_right(1,j));
fprintf(fid,' %g', both_cat_left(1,j));
fprintf(fid,' %g', both_cat_right(1,j));
fprintf(fid, '\n');
end
fclose(fid)
end

