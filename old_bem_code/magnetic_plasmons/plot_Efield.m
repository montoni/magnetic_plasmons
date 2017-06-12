cmapEfield = [ 
    0.0417         0         0
    0.0833         0         0
    0.1250         0         0
    0.1667         0         0
    0.2083         0         0
    0.2500         0         0
    0.2917         0         0
    0.3333         0         0
    0.3750         0         0
    0.4167         0         0
    0.4583         0         0
    0.5000         0         0
    0.5417         0         0
    0.5833         0         0
    0.6250         0         0
    0.6667         0         0
    0.7083         0         0
    0.7500         0         0
    0.7917         0         0
    0.8333         0         0
    0.8750         0         0
    0.9167         0         0
    0.9583         0         0
    1.0000         0         0
    1.0000    0.0417         0
    1.0000    0.0833         0
    1.0000    0.1250         0
    1.0000    0.1667         0
    1.0000    0.2083         0
    1.0000    0.2500         0
    1.0000    0.2917         0
    1.0000    0.3333         0
    1.0000    0.3750         0
    1.0000    0.4167         0
    1.0000    0.4583         0
    1.0000    0.5000         0
    1.0000    0.5417         0
    1.0000    0.5833         0
    1.0000    0.6250         0
    1.0000    0.6667         0
    1.0000    0.7083         0
    1.0000    0.7500         0
    1.0000    0.7917         0
    1.0000    0.8333         0
    1.0000    0.8750         0
    1.0000    0.9167         0
    1.0000    0.9583         0
    1.0000    1.0000         0
    1.0000    1.0000    0.0625
    1.0000    1.0000    0.1250
    1.0000    1.0000    0.1875
    1.0000    1.0000    0.2500
    1.0000    1.0000    0.3125
    1.0000    1.0000    0.3750
    1.0000    1.0000    0.4375
    1.0000    1.0000    0.5000
    1.0000    1.0000    0.5625
    1.0000    1.0000    0.6250
    1.0000    1.0000    0.6875
    1.0000    1.0000    0.7500
    1.0000    1.0000    0.8125
    1.0000    1.0000    0.8750
    1.0000    1.0000    0.9375
    1.0000    1.0000    1.0000];

figure(2)
colormap('jet')
%EELSmap=load('PW_Efield_kxPy_1.35'); 
%EELSmap=load('PW_Efield_unita_kxPy_1.33'); 
%EELSmap=load('PW_Efield_fano_kxPy_1.34'); 
EELSmap=load('Efield_holes_EELS_2.17'); 
pcolor(EELSmap);
colorbar
caxis( [ 0 0.25 ] )
shading interp
shading flat
axis equal

% figure(4)
% colormap('jet')
% %EELSmap=load('PW_Efield_kxPy_1.35'); 
% %EELSmap=load('PW_Efield_unita_kxPy_1.33'); 
% %EELSmap=load('PW_Efield_fano_kxPy_1.34'); 
% EELSmap=load('185B1.48'); 
% pcolor(EELSmap);
% colorbar
% caxis( [ 0 10 ] )
% shading interp
% shading flat
% axis equal
% 
% figure(6)
% colormap(cmapEfield)
% %EELSmap=load('PW_Efield_kxPy_1.35'); 
% EELSmap=load('PW_Efield_unita_kxPy_1.33'); 
% %EELSmap=load('PW_Efield_fano_kxPy_1.34'); 
% %EELSmap=load('PW_Efield_fano_kyPx_1.34'); 
% pcolor(EELSmap);
% caxis( [ 0 70 ] )
% shading interp
% shading flat
% axis equal