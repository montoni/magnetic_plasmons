clear, clc

M1 = load('Dimer_Map_1.33');
M2 = load('Dimer_Map_1.48');

M1 = [M1;flipud(M1)];

M2 = [M2;flipud(M2)];

figure(1)
pcolor(M1)
shading flat
shading interp
axis equal 
axis off

figure(2)
pcolor(M2)
shading flat
shading interp
axis equal 
axis off

figure(3)
pcolor(B1x)
shading flat
shading interp
axis equal 
axis off

figure(4)
pcolor(B2x)
shading flat
shading interp
axis equal 
axis off

figure(5)
pcolor(E1x)
shading flat
shading interp
axis equal 
axis off

figure(6)
pcolor(E2x)
shading flat
shading interp
axis equal 
axis off