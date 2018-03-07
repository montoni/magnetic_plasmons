addpath(genpath('/usr/lusers/montoni/MNPBEM14'));
op2 = bemoptions( 'sim', 'ret', 'interp', 'curv' );
units;
% table of dielectric functions
epsmetal=epstable('silver_normal.dat');
epstab = {epsconst( 1),epsmetal};
val = 5;
%for rad = 1:30
rad = 15;
diameter = 2*rad;
height = 2*(8*rad);
dist = .75*height;
offset = 3*rad/2;
n = 15;
rod = rot(trirod(diameter,height,[n n n],'triangles'),90,[0,1,0]);

p1 = shift(rot(rod,30,[0,0,1]),[-sqrt(3)*dist/2-offset,dist/2,0]);
p2 = shift(rot(rod,150,[0,0,1]),[-sqrt(3)*dist/2-offset,-dist/2,0]);
p3 = shift(rot(rod,90,[0,0,1]),[-offset,0,0]);
p4 = shift(rot(rod,90,[0,0,1]),[offset,0,0]);
p5 = shift(rot(rod,30,[0,0,1]),[sqrt(3)*dist/2+offset,-dist/2,0]);
p6 = shift(rot(rod,150,[0,0,1]),[sqrt(3)*dist/2+offset,dist/2,0]);
Ptot= {p1,p2,p3,p4,p5,p6};
p=comparticle(epstab,Ptot,[2,1;2,1;2,1;2,1;2,1;2,1],1,2,3,4,5,6,op2);
plot(p)
%%
exc = planewaveret( [ 1, 0, 0 ], [ 0, 1, 0 ] );
enei = eV2nm./linspace( 2.25, 2.75, 51 );
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
fid = fopen(strcat('rods_NN_',num2str(rad),'nm'),'wt');
           for j = 1:length(enei)
           fprintf(fid,' %g', A(j,1));
           fprintf(fid,' %g', A(j,2));
           fprintf(fid,' %g', A(j,3));
           fprintf(fid, '\n');
           end
fclose(fid)
%end                              
