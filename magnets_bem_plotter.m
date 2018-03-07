clear all
close all
clc

num = 1;
%A = importdata(strcat('drude_polarization_plots/twomer_spectrum_nm',num2str(num)),' ',1);
%D = importdata(strcat('3mer/sharp_twomer_optical_nm',num2str(num)),' ',1);
%E = importdata(strcat('3mer/sharp_twomer_eels_nm',num2str(num)),' ',1);
clf
theta = linspace(0,2*pi,361);
eV = linspace(2.8,3.8,101);
for ene = linspace(min(eV),max(eV),101);
     B = load(strcat('water_dP_domega/water_beam_rad_',num2str(num),'_ene_',num2str(ene)));
     %B = load(strcat('diff_sca_water_yz/water_yz_rad_',num2str(num),'_ene_',num2str(ene)));
     [sx,sy] = pol2cart( theta(:), B  );
     clf
%      figure(1)
%      plot(eV,A.data(:,3)./max(A.data(:,3)))
%      hold on
%      plot([ene,ene],[0,1],'Color','black')
%      legend('Both Polarizations')
%      hold off
%     
    figure(2)
    plot(sx,sy)
    title(ene)
    xlabel('X-Direction')
    ylabel('Z-Direction')
    pause

end
%%

for num = [1,5,10,15,20,25,30];
    
    A = importdata(strcat('static_optical/vac_stat_twomer_fb_sca_nm',num2str(num)),' ',1);
    B = importdata(strcat('static_optical/vac_stat_twomer_leri_sca_',num2str(num)),' ',1);
    C = importdata(strcat('static_optical/vac_stat_twomer_updo_sca_',num2str(num)),' ',1);
    
    eV = linspace(2.8,3.8,101);
    clf
    %plot(eV,A.data(:,2),eV,A.data(:,3))

    plot(eV,B.data(:,2),eV,C.data(:,2))
    legend('Left','Up')
    pause
end

%%
num = 15;
ene = linspace(2.8,3.8,51);
surface = trisphere(256,1);
for i = 1:length(ene)
    ene(i)
    A = load(strcat('poynt_rad_',num2str(num),'_ene_',num2str(ene(i))));
    clf
    plot(surface,A)
    delete(findall(gcf,'Type','light'))
    pause
end
%%
num = 20;
theta = linspace(0,2*pi,361);
for sep = linspace(0.2,3,15);
    A = load(strcat('anth_beaming/poyntyz_anth_',num2str(num),'_sep_',num2str(sep)));
    [sx,sy] = pol2cart( theta(:),A  );
     clf
%      figure(1)
%      plot(eV,A.data(:,3)./max(A.data(:,3)))
%      hold on
%      plot([ene,ene],[0,1],'Color','black')
%      legend('Both Polarizations')
%      hold off
%     
    figure(2)
    plot(sx,sy)
    title(sep)
    xlabel('X-Direction')
    ylabel('Z-Direction')
    pause
end

%%
for num = [15,30]
    A = importdata(strcat('CL_points/nm',num2str(num),'_2mer_CL_NS_thin_new'),' ',1);
    B = importdata(strcat('CL_points/nm',num2str(num),'_2mer_CL_NN_thin_new'),' ',1);
    C = importdata(strcat('CL_points/nm',num2str(num),'_2mer_CL_x_thin_new'),' ',1);
    
    subplot(3,1,1)
    plot(A.data(:,1),A.data(:,2),A.data(:,1),A.data(:,3))
    legend('NS xz','NS yz')

    subplot(3,1,2)
    plot(B.data(:,1),B.data(:,4),B.data(:,1),B.data(:,5))
    legend('NN xz','NN yz')

    subplot(3,1,3)
    plot(C.data(:,1),C.data(:,2),C.data(:,1),C.data(:,3))
    legend('x xz','x yz')
    pause
end

%%
for num = [15,30]
    A = importdata(strcat('CL_points/nm',num2str(num),'_anth_CL_NSN'),' ',1);
    B = importdata(strcat('CL_points/nm',num2str(num),'_anth_CL_NNN'),' ',1);
    C = importdata(strcat('CL_points/nm',num2str(num),'_anth_CL_x'),' ',1);
    D = importdata(strcat('CL_points/nm',num2str(num),'_anth_CL_N_S'),' ',1);
    
    subplot(4,1,1)
    plot(A.data(:,1),A.data(:,2))
    legend('NSN xz','NSN yz')

    subplot(4,1,2)
    plot(B.data(:,1),B.data(:,5))
    legend('NNN xz','NNN yz')

    subplot(4,1,3)
    plot(C.data(:,1),C.data(:,2),C.data(:,1),C.data(:,3))
    legend('x xz','x yz')
    
    subplot(4,1,4)
    plot(D.data(:,1),D.data(:,2))
    legend('N_S xz','N_S yz')
    
    %figure()
    %blarg = B.data(:,5) - A.data(:,2) - C.data(:,3);
    %plot(B.data(:,1),blarg)
    
    
    
    pause
end

%%
for num = [15,30]
    A = importdata(strcat('CL_points/nm',num2str(num),'_phen_CL_NS1'),' ',1);
    B = importdata(strcat('CL_points/nm',num2str(num),'_phen_CL_NS2'),' ',1);
    C = importdata(strcat('CL_points/nm',num2str(num),'_phen_CL_NN'),' ',1);
    
    subplot(3,1,1)
    plot(A.data(:,1),A.data(:,2),A.data(:,1),A.data(:,3))
    legend('NS xz','NS yz')

    subplot(3,1,2)
    plot(B.data(:,1),B.data(:,2),B.data(:,1),B.data(:,3))
    legend('NN xz','NN yz')

    subplot(3,1,3)
    plot(C.data(:,1),C.data(:,2),C.data(:,1),C.data(:,3))
    legend('x xz','x yz')
    pause
end
%%
clear all
close all
clc

for num = [40]
A = importdata(strcat('greybush_2017/nm',num2str(num),'_19_CL_NN'),' ',1);
B = importdata(strcat('greybush_2017/nm',num2str(num),'_19_CL_NS'),' ',1);
plot(A.data(:,1),A.data(:,2),B.data(:,1),B.data(:,2))
legend('NN','NS')
pause
end

%%
clear all
close all
clc
units
for num = [1,5,10,15,20,25,30]
A = importdata(strcat('diff_hemi_spectra/yx_vac_twomer_points_sca_',num2str(num)),' ',1);
B = importdata(strcat('diff_hemi_spectra/yz_vac_twomer_points_sca_',num2str(num)),' ',1);

plot(eV2nm./A.data(:,1),A.data(:,4)./max(A.data(:,4)),eV2nm./B.data(:,1),B.data(:,4)./max(B.data(:,4)))
legend('up','left')
pause
end

%%
clear all
close all
clc
units
for num = [20,30,40,50,60,70,80,90,100]
A = importdata(strcat('prism_twomer/2mer_test_scattering',num2str(num)),' ',1);

plot(A.data(:,1),A.data(:,3),A.data(:,1),A.data(:,5))
legend('EyBz','EyBx','ExBz','ExBy')
pause
end

%%
clear all
close all
clc

A = load('40nm_19_CL_map_3.13');
B = load('40nm_19_EELS_map_3.13');

AA = [fliplr(A),A];
BB = [fliplr(B),B];

CL = [flipud(AA);AA];
EELS = [flipud(BB);BB];

figure(1)
pcolor(CL./max(max(CL)))
caxis([0,1])
shading flat
shading interp 
figure(2)
pcolor(EELS./max(max(EELS)))
caxis([0,1])
shading flat
shading interp