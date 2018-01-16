addpath(genpath('/usr/lusers/montoni/MNPBEM14'));
clear all
close all
clc
op2 = bemoptions( 'sim', 'ret', 'interp', 'curv' );

epstab = {epsconst( 1 ),epstable('silver_normal.dat')};
for num = [10,20,30,40,50];

numPart = 6;
theta = linspace(0,2*pi,numPart+1);
phi = 2*pi/(numPart*2);
diam = 60;
%diam=5;
separation = 1.5*diam;
RAD = separation/(2*sin(phi));
second_layer = sqrt(3)*separation;
third_layer = 2*separation;
fourth_layer = sqrt(7)*separation;
sphere = trisphere(144,diam);

p0 = sphere;
p1 = shift(sphere, [RAD*cos(theta(1)), RAD*sin(theta(1)), 0]);
p2 = shift(sphere, [RAD*cos(theta(2)), RAD*sin(theta(2)), 0]);
p3 = shift(sphere, [RAD*cos(theta(3)), RAD*sin(theta(3)), 0]);
p4 = shift(sphere, [RAD*cos(theta(4)), RAD*sin(theta(4)), 0]);
p5 = shift(sphere, [RAD*cos(theta(5)), RAD*sin(theta(5)), 0]);
p6 = shift(sphere, [RAD*cos(theta(6)), RAD*sin(theta(6)), 0]);
p7 = shift(sphere, [second_layer*cos(theta(1)+phi), second_layer*sin(theta(1)+phi), 0]);
p8 = shift(sphere, [second_layer*cos(theta(2)+phi), second_layer*sin(theta(2)+phi), 0]);
p9 = shift(sphere, [second_layer*cos(theta(3)+phi), second_layer*sin(theta(3)+phi), 0]);
p10 = shift(sphere, [second_layer*cos(theta(4)+phi), second_layer*sin(theta(4)+phi), 0]);
p11 = shift(sphere, [second_layer*cos(theta(5)+phi), second_layer*sin(theta(5)+phi), 0]);
p12 = shift(sphere, [second_layer*cos(theta(6)+phi), second_layer*sin(theta(6)+phi), 0]);
p13 = shift(sphere, [third_layer*cos(theta(1)), third_layer*sin(theta(1)), 0]);
p14 = shift(sphere, [third_layer*cos(theta(2)), third_layer*sin(theta(2)), 0]);
p15 = shift(sphere, [third_layer*cos(theta(3)), third_layer*sin(theta(3)), 0]);
p16 = shift(sphere, [third_layer*cos(theta(4)), third_layer*sin(theta(4)), 0]);
p17 = shift(sphere, [third_layer*cos(theta(5)), third_layer*sin(theta(5)), 0]);
p18 = shift(sphere, [third_layer*cos(theta(6)), third_layer*sin(theta(6)), 0]);
p19 = shift(sphere, [fourth_layer*cos(theta(1)+phi/2), fourth_layer*sin(theta(1)+phi/2), 0]);
p20 = shift(sphere, [fourth_layer*cos(theta(2)+phi/2), fourth_layer*sin(theta(2)+phi/2), 0]);
p21 = shift(sphere, [fourth_layer*cos(theta(3)+phi/2), fourth_layer*sin(theta(3)+phi/2), 0]);
p22 = shift(sphere, [fourth_layer*cos(theta(4)+phi/2), fourth_layer*sin(theta(4)+phi/2), 0]);
p23 = shift(sphere, [fourth_layer*cos(theta(5)+phi/2), fourth_layer*sin(theta(5)+phi/2), 0]);
p24 = shift(sphere, [fourth_layer*cos(theta(6)+phi/2), fourth_layer*sin(theta(6)+phi/2), 0]);
p25 = shift(sphere, [fourth_layer*cos(theta(1)+3*phi/2), fourth_layer*sin(theta(1)+3*phi/2), 0]);
p26 = shift(sphere, [fourth_layer*cos(theta(2)+3*phi/2), fourth_layer*sin(theta(2)+3*phi/2), 0]);
p27 = shift(sphere, [fourth_layer*cos(theta(3)+3*phi/2), fourth_layer*sin(theta(3)+3*phi/2), 0]);
p28 = shift(sphere, [fourth_layer*cos(theta(4)+3*phi/2), fourth_layer*sin(theta(4)+3*phi/2), 0]);
p29 = shift(sphere, [fourth_layer*cos(theta(5)+3*phi/2), fourth_layer*sin(theta(5)+3*phi/2), 0]);
p30 = shift(sphere, [fourth_layer*cos(theta(6)+3*phi/2), fourth_layer*sin(theta(6)+3*phi/2), 0]);

p = comparticle( epstab, { p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30 },  ...
                 [ 2,1; 2,1; 2,1; 2,1; 2,1; 2,1; 2,1; 2,1; 2,1; 2,1; 2,1; 2,1; 2,1; 2,1; 2,1; 2,1; 2,1; 2,1; 2,1; 2,1; 2,1; 2,1; 2,1; 2,1; 2,1; 2,1; 2,1; 2,1; 2,1; 2,1; 2,1], 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31, op2 );

units;

[ width, vel ] = deal( 0.01, eelsbase.ene2vel( 200e3 ) );

imp = [0, fourth_layer+1.05*(diam/2); 1.05*(diam/2), fourth_layer];

dir1 = [sin(pi/4)*cos(0), sin(pi/4)*sin(0), cos(pi/4)];
dir2 = [sin(pi/4)*cos(pi/2), sin(pi/4)*sin(pi/2), cos(pi/4)];
dir3 = [sin(pi/4)*cos(pi), sin(pi/4)*sin(pi), cos(pi/4)];
dir4 = [sin(pi/4)*cos(3*pi/2), sin(pi/4)*sin(3*pi/2), cos(pi/4)];
dir5 = [sin(0)*cos(0), sin(0)*sin(0), cos(0)];

enei = eV2nm./linspace( 2.5,3.8,231);
ene = linspace(2.5,3.8,231);

bemr = bemsolver( p, op2 );

excr = eelsret( p , imp, width, vel, op2 );

[ psurfret, pbulkret ] = deal( zeros( size( imp, 1 ), numel( enei ) ) );
cathodo   = zeros( size(imp,1),numel( enei ));
NS_cat = zeros( 5,numel( enei ) );
NN_cat = zeros( 5,numel( enei ) );


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


          elseif abs(coords(j,: ) - dir2) < tol
	     NS_cat(2,ien) = NS_cat(2,ien) + dcat1.dprad.struct.dsca(j)*dcat1.p.area(j);
             NN_cat(2,ien) = NN_cat(2,ien) + dcat1.dprad.struct.dsca(j+length(coords))*dcat1.p.area(j); 


          elseif abs(coords(j,: ) - dir3) < tol
             NS_cat(3,ien) = NS_cat(3,ien) + dcat1.dprad.struct.dsca(j)*dcat1.p.area(j);
             NN_cat(3,ien) = NN_cat(3,ien) + dcat1.dprad.struct.dsca(j+length(coords))*dcat1.p.area(j);


          elseif abs(coords(j,: ) - dir4) < tol
             NS_cat(4,ien) = NS_cat(4,ien) + dcat1.dprad.struct.dsca(j)*dcat1.p.area(j);
             NN_cat(4,ien) = NN_cat(4,ien) + dcat1.dprad.struct.dsca(j+length(coords))*dcat1.p.area(j);


          elseif abs(coords(j,: ) - dir5) < tol
             NS_cat(5,ien) = NS_cat(5,ien) + dcat1.dprad.struct.dsca(j)*dcat1.p.area(j);
             NN_cat(5,ien) = NN_cat(5,ien) + dcat1.dprad.struct.dsca(j+length(coords))*dcat1.p.area(j);


end
    
end
  
multiWaitbar( 'BEM solver', ien / numel( enei ) );
end
multiWaitbar( 'CloseAll' );

fid = fopen([strcat('nm',num2str(num),'_31_CL_NS')],'wt');
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

fid = fopen([strcat('nm',num2str(num),'_31_CL_NN')],'wt');
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


end

