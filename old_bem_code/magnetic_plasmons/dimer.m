%%  initialization
%  options for BEM simulation
clear, clc

%  options for BEM simulation
op1 = bemoptions( 'sim', 'stat', 'interp', 'curv' );
op2 = bemoptions( 'sim', 'ret', 'interp', 'curv' );


epsmetal = epsdrude('Ag');
epsmetal2 = epsdrude('Au');
epstab = {epsconst( 2 ),epsmetal,epsmetal2};
% nanosphere 
radius1 = 1;
radius2 = 10;

val = 5;

%for LoopCount = 1:15
height1 = 50;
r = radius2+12;
diameter1 =2*radius1;
diameter2 = 2*radius2;
p = trirod(diameter1,height1,[15 15 15]);
ps = trisphere(256,diameter2);


theta = linspace(pi/6,2*pi+pi/6,7);

p1 = rot(ps,90,[cos(theta(1)),sin(theta(1)),0]);
p1 = shift(p1,[r*cos(theta(1)), r*sin(theta(1)), 0]);

p1s = shift(ps,[r*cos(theta(6)), r*sin(theta(6)), 0]);
%p1s = shift(p1s,[0, r/2, 0]);

p2s = shift(ps,[r*cos(theta(7)), r*sin(theta(7)), 0]);
%p2s = shift(p2s,[0, -r/2, 0]);

p2 = rot(ps,90,[cos(theta(2)),sin(theta(2)),0]);
p2 = shift(p2,[r*cos(theta(2)), r*sin(theta(2)), 0]);
p3 = rot(ps,90,[cos(theta(3)),sin(theta(3)),0]);
p3 = shift(p3,[r*cos(theta(3)), r*sin(theta(3)), 0]);
p4 = rot(ps,90,[cos(theta(4)),sin(theta(4)),0]);
p4 = shift(p4,[r*cos(theta(4)), r*sin(theta(4)), 0]);
p5 = rot(ps,90,[cos(theta(5)),sin(theta(5)),0]);
p5 = shift(p5,[r*cos(theta(5)), r*sin(theta(5)), 0]);



p6 = rot(ps,90,[cos(-theta(2)),sin(-theta(2)),0]);
p6 = shift(p6,[sqrt(3)*r, 0, 0]);
p6 = shift(p6,[-r*cos(-theta(2)), -r*sin(-theta(2)), 0]);
 
p7 = rot(ps,90,[cos(-theta(3)),sin(-theta(3)),0]);
p7 = shift(p7,[sqrt(3)*r, 0, 0]);
p7 = shift(p7,[-r*cos(-theta(3)), -r*sin(-theta(3)), 0]);


p8 = rot(ps,90,[cos(-theta(4)),sin(-theta(4)),0]);
p8 = shift(p8,[sqrt(3)*r, 0, 0]);
p8 = shift(p8,[-r*cos(-theta(4)), -r*sin(-theta(4)), 0]);
 
p9 = rot(ps,90,[cos(-theta(5)),sin(-theta(5)),0]);
p9 = shift(p9,[sqrt(3)*r, 0, 0]);
p9 = shift(p9,[-r*cos(-theta(5)), -r*sin(-theta(5)), 0]);




R = sqrt((height1/2+4).^2 + r^2);
phi = acos(r/R);
gamma = theta(3) + phi;

pbeam = shift(trisphere(32,1),[R*cos(gamma), R*sin(gamma),0]);


Ptot= {p1s,p2s,p2,p3,p4,p5,p6,p7,p8,p9};

p=comparticle(epstab,Ptot,[2,1;2,1;2,1;2,1;2,1;2,1;2,1;2,1;2,1;2,1],1,2,3,4,5,6,7,8,9,10,op1);
% p = rotate(p, 
clf
plot(p)
xlabel('x');
ylabel('y');
zlabel('z');
view([0,90])

%Update name
%name = 'Dimer';


%% 
%  TM mode, excitation from above
dir = [ 0, 0, 1 ; 0, 1, 0 ; 1, 0, 0 ];
pol = [ 0, 1, 0 ; 1, 0, 0 ; 0, 1, 0 ];
%  photon wavelengths
units;
enei = eV2nm./linspace( 1, 3.5, 101 );
ene = transpose(eV2nm./enei);
%%  BEM solver
%  initialize BEM solver
bem = bemsolver( p, op1 );
%  initialize plane wave excitation
exc = planewave( pol, dir, op1 );
%  scattering cross section
sca = zeros( numel( enei ), size( dir, 1 ) );
ext = zeros( numel( enei ), size( dir, 1 ) );
%%
multiWaitbar( 'BEM solver', 0, 'Color', 'g', 'CanCancel', 'on' );
%  loop over wavelengths
for ien = 1 : length( enei )
  %  surface charges
  sig = bem \ exc( p, enei( ien ) );
  multiWaitbar( 'BEM solver', ien / numel( enei ) );
  
  %  scattering cross section
  ext( ien, : ) = exc.ext( sig );
  sca( ien, : ) = exc.sca( sig );
 
  multiWaitbar( 'BEM solver', ien / numel( enei ) );
end


%  close waitbar
multiWaitbar( 'CloseAll' );

%%  final plot
figure(3)
plot( ene, sca );

xlabel( 'Energy (eV)' );
ylabel( 'Scattering cross section (nm^2)' );
