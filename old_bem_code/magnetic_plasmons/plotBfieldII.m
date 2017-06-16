cmapBfield =[ 0.8000         0    0.2000
    1.0000    0.0159         0
    1.0000    0.0510    0.0357
    1.0000    0.0862    0.0714
    1.0000    0.1213    0.1071
    1.0000    0.1565    0.1429
    1.0000    0.1916    0.1786
    1.0000    0.2268    0.2143
    1.0000    0.2619    0.2500
    1.0000    0.2971    0.2857
    1.0000    0.3322    0.3214
    1.0000    0.3673    0.3571
    1.0000    0.4025    0.3929
    1.0000    0.4376    0.4286
    1.0000    0.4728    0.4643
    1.0000    0.5079    0.5000
    1.0000    0.5431    0.5357
    1.0000    0.5782    0.5714
    1.0000    0.6134    0.6071
    1.0000    0.6485    0.6429
    1.0000    0.6837    0.6786
    1.0000    0.7188    0.7143
    1.0000    0.7540    0.7500
    1.0000    0.7891    0.7857
    1.0000    0.8243    0.8214
    1.0000    0.8594    0.8571
    1.0000    0.8946    0.8929
    1.0000    0.9297    0.9286
    1.0000    0.9649    0.9643
    1.0000    1.0000    1.0000
    1.0000    1.0000    1.0000
    1.0000    1.0000    1.0000
    1.0000    1.0000    1.0000
    1.0000    1.0000    1.0000
    0.9655    0.9655    1.0000
    0.9310    0.9310    1.0000
    0.8966    0.8966    1.0000
    0.8621    0.8621    1.0000
    0.8276    0.8276    1.0000
    0.7931    0.7931    1.0000
    0.7586    0.7586    1.0000
    0.7241    0.7241    1.0000
    0.6897    0.6897    1.0000
    0.6552    0.6552    1.0000
    0.6207    0.6207    1.0000
    0.5862    0.5862    1.0000
    0.5517    0.5517    1.0000
    0.5172    0.5172    1.0000
    0.4828    0.4828    1.0000
    0.4483    0.4483    1.0000
    0.4138    0.4138    1.0000
    0.3793    0.3793    1.0000
    0.3448    0.3448    1.0000
    0.3103    0.3103    1.0000
    0.2759    0.2759    1.0000
    0.2414    0.2414    1.0000
    0.2069    0.2069    1.0000
    0.1724    0.1724    1.0000
    0.1379    0.1379    1.0000
    0.1034    0.1034    1.0000
    0.0690    0.0690    1.0000
    0.0345    0.0345    1.0000
         0         0    1.0000
    0.0784    0.1686    0.5490];

figure(1)
colormap(cmapBfield)
%EELSmap=load('PW_Bfield_kyPx_1.33'); 
%EELSmap=load('PW_Bfield_unita_kxPy_1.33'); 
EELSmap=load('PW_Bfield_fano_kxPy_1.375'); 
pcolor(EELSmap);
colorbar
caxis( [-5 5] )
shading interp
shading flat
axis equal

figure(2)
colormap(cmapBfield)
%EELSmap=load('PW_Bfield_kyPx_1.33'); 
%EELSmap=load('PW_Bfield_unita_kxPy_1.33'); 
EELSmap=load('185B1.48'); 
pcolor(EELSmap);
colorbar
caxis( [-20 20] )
shading interp
shading flat
axis equal

% figure(2)
% colormap(cmapBfield)
% %EELSmap=load('PW_Bfield_kyPx_1.33'); 
% EELSmap=load('PW_Bfield_unita_kxPy_1.33'); 
% %EELSmap=load('PW_Bfield_fano_kxPy_1.34'); 
% pcolor(EELSmap);
% caxis( [-5 5] )
% shading interp
% shading flat
% axis equal
% 
% figure(3)
% colormap(cmapBfield)
% %EELSmap=load('PW_Bfield_kyPx_1.33'); 
% EELSmap=load('PW_Bfield_kyPx_1.33'); 
% %EELSmap=load('PW_Bfield_fano_kxPy_1.34'); 
% pcolor(EELSmap);
% caxis( [-5 5] )
% shading interp
% shading flat
% axis equal
% 
% figure(4)
% colormap(cmapBfield)
% %EELSmap=load('PW_Bfield_kyPx_1.33'); 
% EELSmap=load('PW_Bfield_kxPy_1.35'); 
% %EELSmap=load('PW_Bfield_fano_kxPy_1.34'); 
% pcolor(EELSmap);
% caxis( [-10 10] )
% shading interp
% shading flat
% axis equal
% 
% figure(6)
% colormap(cmapBfield)
% %EELSmap=load('PW_Bfield_kyPx_1.33'); 
% %EELSmap=load('PW_Bfield_unita_kxPy_1.33'); 
% EELSmap=load('PW_Bfield_fano_kxPy_1.43'); 
% pcolor(EELSmap);
% caxis( [-8 8] )
% shading interp
% shading flat
% axis equal
%%
save('map.dat','map','-ascii');