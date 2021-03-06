addpath(genpath('/usr/lusers/montoni/MNPBEM14'));
clear all
close all
clc
op2 = bemoptions( 'sim', 'ret', 'interp', 'curv' );

epstab = {epsconst( 1 ),epstable('silver_normal.dat')};
for num = [1,5,10,15,20,25,30];
sep = [2.5];
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

dir1 = [sin(pi/4)*cos(0), sin(pi/4)*sin(0), cos(pi/4)];
dir2 = [sin(pi/4)*cos(pi/2), sin(pi/4)*sin(pi/2), cos(pi/4)];
dir3 = [sin(pi/4)*cos(pi), sin(pi/4)*sin(pi), cos(pi/4)];
dir4 = [sin(pi/4)*cos(3*pi/2), sin(pi/4)*sin(3*pi/2), cos(pi/4)];
dir5 = [sin(0)*cos(0), sin(0)*sin(0), cos(0)];

enei = eV2nm./linspace( 2.5,3.8,231);
ene = linspace(2.5,3.8,231);

bemr = bemsolver( pret, op2 );

excr = eelsret( pret , imp, width, vel, op2 );

[ psurfret, pbulkret ] = deal( zeros( size( imp, 1 ), numel( enei ) ) );
cathodo   = zeros( size(imp,1),numel( enei ));
NS_cat = zeros( 5,numel( enei ) );
NN_cat = zeros( 5,numel( enei ) );
x_cat = zeros( 5,numel( enei ) );

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
  tol = [.075,.075,.075];
for j = 1:length(coords)
    % The following if statement determines if the coordinates of the
    % surface elements of the sphere at infinity are positive or negative
    % for the direction of propogation
	  if abs(coords(j,: ) - dir1) < tol
	     NS_cat(1,ien) = NS_cat(1,ien) + dcat1.dprad.struct.dsca(j)*dcat1.p.area(j);
             NN_cat(1,ien) = NN_cat(1,ien) + dcat1.dprad.struct.dsca(j+length(coords))*dcat1.p.area(j);
             x_cat(1,ien) = x_cat(1,ien) + dcat1.dprad.struct.dsca(j+length(coords)+length(coords))*dcat1.p.area(j); 

          elseif abs(coords(j,: ) - dir2) < tol
	     NS_cat(2,ien) = NS_cat(2,ien) + dcat1.dprad.struct.dsca(j)*dcat1.p.area(j);
             NN_cat(2,ien) = NN_cat(2,ien) + dcat1.dprad.struct.dsca(j+length(coords))*dcat1.p.area(j); 
             x_cat(2,ien) = x_cat(2,ien) + dcat1.dprad.struct.dsca(j+length(coords)+length(coords))*dcat1.p.area(j);

          elseif abs(coords(j,: ) - dir3) < tol
             NS_cat(3,ien) = NS_cat(3,ien) + dcat1.dprad.struct.dsca(j)*dcat1.p.area(j);
             NN_cat(3,ien) = NN_cat(3,ien) + dcat1.dprad.struct.dsca(j+length(coords))*dcat1.p.area(j);
             x_cat(3,ien) = x_cat(3,ien) + dcat1.dprad.struct.dsca(j+length(coords)+length(coords))*dcat1.p.area(j);

          elseif abs(coords(j,: ) - dir4) < tol
             NS_cat(4,ien) = NS_cat(4,ien) + dcat1.dprad.struct.dsca(j)*dcat1.p.area(j);
             NN_cat(4,ien) = NN_cat(4,ien) + dcat1.dprad.struct.dsca(j+length(coords))*dcat1.p.area(j);
             x_cat(4,ien) = x_cat(4,ien) + dcat1.dprad.struct.dsca(j+length(coords)+length(coords))*dcat1.p.area(j);

          elseif abs(coords(j,: ) - dir5) < tol
             NS_cat(5,ien) = NS_cat(5,ien) + dcat1.dprad.struct.dsca(j)*dcat1.p.area(j);
             NN_cat(5,ien) = NN_cat(5,ien) + dcat1.dprad.struct.dsca(j+length(coords))*dcat1.p.area(j);
             x_cat(5,ien) = x_cat(5,ien) + dcat1.dprad.struct.dsca(j+length(coords)+length(coords))*dcat1.p.area(j);

end
    
end
  
multiWaitbar( 'BEM solver', ien / numel( enei ) );
end
multiWaitbar( 'CloseAll' );

fid = fopen([strcat('nm',num2str(num),'_2mer_CL_NS_close')],'wt');
fprintf(fid, ' %s', 'Energy(eV)     NS_101     NS_0-11     NS_-101     NS_011     NS_000');
fprintf(fid, '\n');
for j = 1:length(ene)
          fprintf(fid,' %g', ene(j));
fprintf(fid,' %g', NS_cat(1,j));
fprintf(fid,' %g', NS_cat(2,j));
fprintf(fid,' %g', NS_cat(3,j));
fprintf(fid,' %g', NS_cat(4,j));
fprintf(fid,' %g', NS_cat(5,j));
fprintf(fid, '\n');
end
fclose(fid)

fid = fopen([strcat('nm',num2str(num),'_2mer_CL_NN_close')],'wt');
fprintf(fid, ' %s', 'Energy(eV)     NN_101     NN_0-11     NN_-101     NN_011     NN_000');
fprintf(fid, '\n');
for j = 1:length(ene)
          fprintf(fid,' %g', ene(j));
fprintf(fid,' %g', NN_cat(1,j));
fprintf(fid,' %g', NN_cat(2,j));
fprintf(fid,' %g', NN_cat(3,j));
fprintf(fid,' %g', NN_cat(4,j));
fprintf(fid,' %g', NN_cat(5,j));
fprintf(fid, '\n');
end
fclose(fid)

fid = fopen([strcat('nm',num2str(num),'_2mer_CL_x_close')],'wt');
fprintf(fid, ' %s', 'Energy(eV)     x_101     x_0-11     x_-101     x_011     x_000');
fprintf(fid, '\n');
for j = 1:length(ene)
          fprintf(fid,' %g', ene(j));
fprintf(fid,' %g', x_cat(1,j));
fprintf(fid,' %g', x_cat(2,j));
fprintf(fid,' %g', x_cat(3,j));
fprintf(fid,' %g', x_cat(4,j));
fprintf(fid,' %g', x_cat(5,j));
fprintf(fid, '\n');
end
fclose(fid)


end

