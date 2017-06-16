clear all
close all
clc
%  final plot
%  x and y limits
ene = linspace(3.4,3.5,21);
for ien = 1:length(ene)
eels = load(strcat('15_nm_EELS_map_',num2str(ene(ien))));
eels = [eels(:,1:60);eels(:,62:121)];
if any(eels < .5*max(max(eels)))
    eels(eels < .5*max(max(eels))) = 0;
end
pcolor(eels)
colormap('hot')
axis equal
axis tight
shading flat
shading interp
caxis([0 max(max(eels))])
title(ene(ien))
colorbar
drawnow
pause
end
