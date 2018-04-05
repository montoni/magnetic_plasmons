addpath(genpath('/usr/lusers/montoni/MNPBEM14'));
op2 = bemoptions( 'sim', 'ret', 'interp', 'curv' );
units;
% table of dielectric functions
epsmetal=epstable('silver_022717.dat');
epstab = {epsconst( 1),epsmetal};
val = 5;
for rad = 11:30
diameter = 2*rad;
height = 2*(2.5*rad);
dist = height;
n = 15;
rod = rot(trirod(diameter,height,[n n n],'triangles'),90,[0,1,0]);
p1 = shift(rot(rod,150,[0,0,1]),[-sqrt(3)*dist-dist/4,dist/2,0]);
p2 = shift(rot(rod,30,[0,0,1]),[-sqrt(3)*dist-dist/4,-dist/2,0]);
p3 = shift(rot(rod,90,[0,0,1]),[-dist/4,0,0]);
p4 = shift(rot(rod,90,[0,0,1]),[dist/4,0,0]);
p5 = shift(rot(rod,150,[0,0,1]),[sqrt(3)*dist+dist/4,-dist/2,0]);
p6 = shift(rot(rod,30,[0,0,1]),[sqrt(3)*dist+dist/4,dist/2,0]);
Ptot= {p1,p2,p3,p4,p5,p6};
p=comparticle(epstab,Ptot,[2,1;2,1;2,1;2,1;2,1;2,1],1,2,3,4,5,6,op2);
exc = planewaveret( [ 1, 0, 0 ], [ 0, 1, 0 ] );
enei = eV2nm./linspace( 2, 2.7, 71 );
bem = bemret( p,op2 );
sca = zeros( length( enei ), 1 );
ext = zeros( length( enei ), 1 );
for ien = 1 : length( enei )
            sig = bem \ exc( p, enei( ien ) );
sca( ien, : ) = exc.sca( sig );
ext( ien, : ) = exc.ext( sig );
  count = ien
end

A(:,1) = eV2nm./enei';
A(:,2) = ext;
A(:,3) = sca;
fid = fopen(strcat('closer_rods_NN_',num2str(rad),'nm'),'wt');
           for j = 1:length(enei)
           fprintf(fid,' %g', A(j,1));
           fprintf(fid,' %g', A(j,2));
           fprintf(fid,' %g', A(j,3));
           fprintf(fid, '\n');
           end
fclose(fid)
end                              
