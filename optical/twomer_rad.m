%addpath(genpath('/usr/lusers/montoni/MNPBEM14'));
op2 = bemoptions( 'sim', 'ret', 'interp', 'curv' );

epstab = {epsconst( 1 ),epstable('silver_normal.dat')};
num = 20;
%for sep = linspace(2,4,11);

numPart = 6;
theta = linspace(0,2*pi,numPart+1);
phi = 2*pi/(numPart*2);
diam = 55;
%diam=5;
separation = 1.02*diam;
RAD = separation/(2*sin(phi));
second_layer = sqrt(3)*separation;
sphere = trisphere(256,diam);
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

p = comparticle( epstab, { p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12 },  ...
                 [ 2,1; 2,1; 2,1; 2,1; 2,1; 2,1; 2,1; 2,1; 2,1; 2,1; 2,1; 2,1; 2,1], 1,2,3,4,5,6,7,8,9,10,11,12,13, op2 );



plot(p)
% p = comparticle(epstab,{ps},[2,1],1,op2);
dir = [ 0, -1, 0  ];
pol = [ 1, 0, 0  ];
units;
%%
enei = eV2nm./linspace( 2.8,3.8,51 );
%enei = eV2nm./3.05;
ene = transpose(eV2nm./enei);
bem = bemsolver( p, op2 );
exc = planewave( pol, dir, op2 );
sca  = zeros( numel( enei ), size( dir, 1 ) );
sca2 = zeros( numel( enei ), size( dir, 1 ) );
ext = zeros( numel( enei ), size( dir, 1 ) );
multiWaitbar( 'BEM solver', 0, 'Color', 'g', 'CanCancel', 'on' );
theta = linspace(0,pi/2,181);
phi = linspace(0,2*pi,361);
for ien = 1 : length( enei )
  ien
  sig = bem \ exc( p, enei( ien ) );
  sigma = sig.sig2;
   r = [sqrt((p.pos(1:length(p.pos)/10,1)).^2 + ...
            (p.pos(1:length(p.pos)/10,2)).^2 + ...
            (p.pos(1:length(p.pos)/10,3)).^2) ; ...
            sqrt((p.pos(length(p.pos)/10+1:length(p.pos)*2/10,1)).^2 + ...
                 (p.pos(length(p.pos)/10+1:length(p.pos)*2/10,2)).^2 + ...
                 (p.pos(length(p.pos)/10+1:length(p.pos)*2/10,3)).^2) ; ...
            sqrt((p.pos(length(p.pos)*2/10+1:length(p.pos)*3/10,1)).^2 + ...
                 (p.pos(length(p.pos)*2/10+1:length(p.pos)*3/10,2)).^2 + ...
                 (p.pos(length(p.pos)*2/10+1:length(p.pos)*3/10,3)).^2) ; ...
            sqrt((p.pos(length(p.pos)*3/10+1:length(p.pos)*4/10,1)).^2 + ...
                 (p.pos(length(p.pos)*3/10+1:length(p.pos)*4/10,2)).^2 + ...
                 (p.pos(length(p.pos)*3/10+1:length(p.pos)*4/10,3)).^2) ; ...
            sqrt((p.pos(length(p.pos)*4/10+1:length(p.pos)*5/10,1)).^2 + ...
                 (p.pos(length(p.pos)*4/10+1:length(p.pos)*5/10,2)).^2 + ...
                 (p.pos(length(p.pos)*4/10+1:length(p.pos)*5/10,3)).^2) ; ...
            sqrt((p.pos(length(p.pos)*5/10+1:length(p.pos)*6/10,1)).^2 + ...
                 (p.pos(length(p.pos)*5/10+1:length(p.pos)*6/10,2)).^2 + ...
                 (p.pos(length(p.pos)*5/10+1:length(p.pos)*6/10,3)).^2) ; ...
            sqrt((p.pos(length(p.pos)*6/10+1:length(p.pos)*7/10,1)).^2 + ...
                 (p.pos(length(p.pos)*6/10+1:length(p.pos)*7/10,2)).^2 + ...
                 (p.pos(length(p.pos)*6/10+1:length(p.pos)*7/10,3)).^2) ; ...
            sqrt((p.pos(length(p.pos)*7/10+1:length(p.pos)*8/10,1)).^2 + ...
                 (p.pos(length(p.pos)*7/10+1:length(p.pos)*8/10,2)).^2 + ...
                 (p.pos(length(p.pos)*7/10+1:length(p.pos)*8/10,3)).^2); ...
            sqrt((p.pos(length(p.pos)*8/10+1:length(p.pos)*9/10,1)).^2 + ...
                 (p.pos(length(p.pos)*8/10+1:length(p.pos)*9/10,2)).^2 + ...
                 (p.pos(length(p.pos)*8/10+1:length(p.pos)*9/10,3)).^2) ; ...
            sqrt((p.pos(length(p.pos)*9/10+1:length(p.pos),1)).^2 + ...
                 (p.pos(length(p.pos)*9/10+1:length(p.pos),2)).^2 + ...
                 (p.pos(length(p.pos)*9/10+1:length(p.pos),3)).^2)] ; 
         q = p.area.*(sigma);
         Q = [q,q,q];
         r_ = [r,r,r];
         R = p.pos./r_;
         d = sum(Q.*R);
         for j = 1:length(theta)
             for k = 1:length(phi)
                 n = [sin(theta(j))*cos(phi(k)), sin(theta(j))*sin(phi(k)),cos(theta(j))];
                 n = n/norm(n);
             sca2(ien) = sca2(ien) +  1/(4*pi*3e8^3)*dot(cross(n,-(eV2nm./enei(ien))^2*d),...
                                                         cross(n,-(eV2nm./enei(ien))^2*d));
             end
         end
  ext( ien, : ) = exc.ext( sig );
  sca( ien, : ) = exc.sca( sig );
  multiWaitbar( 'BEM solver', ien / numel( enei ) );
  clear Efx Efy Efz Bfx Bfy Bfz EFieldf BFieldf
  clear Ebx Eby Ebz Bbx Bby Bbz EFieldb BFieldb
end

 

multiWaitbar( 'CloseAll' );

 



plot(eV2nm./enei,sca2./max(sca2),eV2nm./enei,sca./max(sca))
legend('Backward','Everywhere')

%%
units;
enei = eV2nm./3.3;
%  set up BEM solver
bem = bemsolver( p, op2 );
units
%  plane wave excitation
exc = planewave( pol, dir, op2 );
%  light wavelength in vacuum
%  allocate scattering and extinction cross sections
sca = zeros( length( enei ), 2 );
ext = zeros( length( enei ), 2 );
%  loop over wavelengths

for ien = 1 : length( enei )
  %  surface charge
	    sig = bem \ exc( p, enei( ien ) );
%sig2 = bem \ exc2( p, enei(ien));
  %  scattering and extinction cross section
  sca( ien, : ) = exc.sca( sig );
ext( ien, : ) = exc.ext( sig );
end

theta = linspace( 0, pi, 181 );
phi = linspace(0, 2*pi, 181);

[theta2, phi2] = meshgrid(theta, phi);
%  directions for emission
dir = surface.pos/250;%[ sin(theta2(:)).*cos( phi2(:) ), sin(theta2(:)) .* sin(phi2(:)), cos( theta2(:) ) ];
%  set up spectrum object
spec = spectrum( dir, op2 );
exc
%  farfield radiation
fpl = farfield( spec, sig );
%  norm of Poynting vector
spl = vecnorm( 0.5 * real( cross( fpl.e, conj( fpl.h ), 2 ) ) );
%poynt_angle = [theta2(:),phi2(:),spl];
plot(surface,spl)
%figure()
%  plot electric field
%imagesc( x( : ), z( : ), log10( ee1 ) );  hold on
%plot( [ min( x( : ) ), max( x( : ) ) ], [ 0, 0 ], 'w--' );

%  Cartesian coordinates of Poynting vector
%[ sx, sy ] = pol2cart( theta2, spl / max( spl ) );
%  overlay with Poynting vector
%clf
%plot( sx, sy, 'b-', 'LineWidth', 1.5 );
%xlabel('Y-Direction')
%ylabel('Z-Direction')
%fid = fopen(strcat('pol_30_sep',num2str(sep-2)),'wt');
%fprintf(fid, ' %g', sx);
%fprintf(fid, ' %g', sy);

%end
light
