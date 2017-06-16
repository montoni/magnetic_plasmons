clear all
close all
clc

units;
ene = linspace(3.1,3.5,201);
enei = eV2nm./ene;

for i = 1:length(ene)
    field = load(strcat('NS_Bfield_30_nm_',num2str(enei(i))));
    ene(i)
    clf
    pcolor(field)
    colorbar
    colormap('jet')
    shading flat
    shading interp
    pause
end