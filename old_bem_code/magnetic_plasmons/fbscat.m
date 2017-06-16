op = bemoptions( 'sim', 'ret', 'interp', 'curv' );

 %  table of dielectric functions
epstab  = { epsconst( 1 ), epstable( 'silver.dat' )};

%  X-LENGTH, Y-LENGTH OF RODS
rodlen_left  = [ 70, 225];
rodlen_right = [ 70, 225];
height = 25;
h=0.25;

% HEIGHT OF RODS
rodedge = edgeprofile( height, 10, 'mode', '00', 'min', 0 );

% MAKE ROD POLYGON
rod_left  = polygon( 40, 'size', rodlen_left  ); 
rod_right = polygon( 40, 'size', rodlen_right );

% RADIUS OF ONE RING
radius1 = 100;
radius2 = 100;

% LOCATION OF CENTER OF ONE RING
center = 160;
del = 0;
disty = 0;

% ANGLES
theta = [0, pi/3-del, 2*pi/3, pi, 4*pi/3, 5*pi/3+del];

hdata = struct( 'hmax', 8 );

[ p, poly ] = tripolygon( rod_left, rodedge);
[ pr, polyr ] = tripolygon( rod_right, rodedge);

p4 = shift(p, [radius2*cos(theta(4))-center, radius2*sin(theta(4)), h]);
p5 = shift(rot(p, 60), [radius2*cos(theta(2))-center, radius2*sin(theta(2))+disty, h]);
p6 = shift(rot(p, 120), [radius2*cos(theta(6))-center, radius2*sin(theta(6))-disty, h]);

p = comparticle( epstab, { p4, p5, p6},  ... 
 [ 2,1; 2,1; 2,1 ], 1,2,3, op );
plot(p)
%%
dir = [ 0, 1, 0  ];
pol = [ 1, 0, 0  ];

rad = 1e6;
[theta,phi] = meshgrid( linspace(0,pi,81),linspace(0,pi,161));
xf = rad*sin(theta).*cos(phi);
yf = rad*sin(theta).*sin(phi);
zf = rad*cos(theta);

[theta,phi] = meshgrid( linspace(0,pi,81),linspace(pi,2*pi,161));
xb = rad*sin(theta).*cos(phi);
yb = rad*sin(theta).*sin(phi);
zb = rad*cos(theta);

ptf = compoint( p, [ xf(:), yf(:), zf(:)], 'mindist', 0.5);
ptb = compoint( p, [ xb(:), yb(:), zb(:)], 'mindist', 0.5);

gf = compgreenret(ptf, p);
gb = compgreenret(ptb, p);


% for j = 1:length(x(:,1))
%     hold on
% scatter3(xb(j,:),yb(j,:),zb(j,:))
% xlabel('x')
% ylabel('y')
% zlabel('z')
% end

units;
enei = eV2nm./linspace( 1,2,101 );
%enei = eV2nm./3.05;
ene = transpose(eV2nm./enei);

bem = bemsolver( p, op );
exc = planewave( pol, dir, op );

sca = zeros( numel( enei ), size( dir, 1 ) );
ext = zeros( numel( enei ), size( dir, 1 ) );
Powerf = zeros( numel( enei ), size( dir, 1 ) );
Powerb = zeros( numel( enei ), size( dir, 1 ) );


for ien = 1 : length( enei )
  sig = bem \ exc( p, enei( ien ) );
  fieldf = gf.field(sig);
  fieldb = gb.field(sig);
  
  Efx = ptf((fieldf.e(:,1)));
  Efy = ptf((fieldf.e(:,2)));
  Efz = ptf((fieldf.e(:,3)));
  
  Bfx = ptf((fieldf.h(:,1)));
  Bfy = ptf((fieldf.h(:,2)));
  Bfz = ptf((fieldf.h(:,3)));
  
  EFieldf = [Efx,Efy,Efz];
  BFieldf = [Bfx,Bfy,Bfz];

  Ebx = ptb((fieldb.e(:,1)));
  Eby = ptb((fieldb.e(:,2)));
  Ebz = ptb((fieldb.e(:,3)));
  
  Bbx = ptb((fieldb.h(:,1)));
  Bby = ptb((fieldb.h(:,2)));
  Bbz = ptb((fieldb.h(:,3)));
  
  EFieldb = [Ebx,Eby,Ebz];
  BFieldb = [Bbx,Bby,Bbz];
  
  for j = 1:length(Efx);
  Powerf(ien) = Powerf(ien) + norm(cross(EFieldf(j,:),BFieldf(j,:)));
  Powerb(ien) = Powerb(ien) + norm(cross(EFieldb(j,:),BFieldb(j,:)));
  end
  
  ext( ien, : ) = exc.ext( sig );
  sca( ien, : ) = exc.sca( sig );
  
  clear Efx Efy Efz Bfx Bfy Bfz EFieldf BFieldf
  clear Ebx Eby Ebz Bbx Bby Bbz EFieldb BFieldb
end

fid = fopen(['SpectrumS-rods-new'],'wt');
for j = 1:length(ene)
          fprintf(fid,' %g', ene(j));
          fprintf(fid,' %g', sca(j));
          fprintf(fid,' %g', Powerf(j));
          fprintf(fid,' %g', Powerb(j));
          fprintf(fid, '\n');
end
fclose(fid)


