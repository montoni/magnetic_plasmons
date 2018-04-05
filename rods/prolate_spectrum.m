%    computes the scattering cross section for different light wavelengths
%    using the full Maxwell equations, and compares the results with Mie
%    theory.
%
%  Runtime on my computer:  7.4 sec.

%% Initialization
clear all
close all
clc
units;
%  options for BEM simulation
op = bemoptions( 'sim', 'ret', 'interp', 'curv' );
for rad = 20:30;
%  table of dielectric functions
epstab  = { epsconst( 1 ), epstable( 'silver_normal.dat' )  };
minor = rad;
major = 5*rad;

rod = prolate(major,minor);

p = comparticle( epstab, { rod},  ...
                   [  2, 1; ], 1, op );

%plot(p, 'EdgeColor','b')
ene = linspace(1,4,151);
enei = eV2nm./ene;


%  set up BEM solver
bem = bemsolver( p, op );
units
%  plane wave excitation
exc = planewave( [ 0, 0, 1], [ 0, 1, 0], op );
%exc2 = planewave( [ 0, 1, 0 ], [ 0, 0, 1], op1 );
%  light wavelength in vacuum

%  allocate scattering and extinction cross sections
sca = zeros( length( enei ), 2 );
ext = zeros( length( enei ), 2 );


multiWaitbar( 'BEM solver', 0, 'Color', 'g', 'CanCancel', 'on' );
%  loop over wavelengths

for ien = 1 : length( enei )
  %  surface charge
  sig = bem \ exc( p, enei( ien ) );
  %sig2 = bem \ exc2( p, enei(ien));
  %  scattering and extinction cross sections
  
  
  sca( ien, : ) = exc.sca( sig );
  ext( ien, : ) = exc.ext( sig );
  multiWaitbar( 'BEM solver', ien / numel( enei ) );
end
%  close waitbar
multiWaitbar( 'CloseAll' );



plot(  eV2nm./enei, ext(:,1)-sca(:,1), eV2nm./enei, sca(:,1), 'linewidth', 2 );

xlabel( 'Wavelength (nm)' );
ylabel( 'cross section (nm^2)' );
legend( 'Absorption', 'Scattering')
%save('ext_Pt.dat','ext','-ascii')


fid = fopen(strcat('prolate_5_to_1_spec_nm',num2str(rad)),'wt');
%fprintf(fid, ' %s', 'Energy(eV)     Ext_long     Sca_long');
for i = 1:length(ene)
    fprintf(fid, ' %g', ene(i));
    fprintf(fid, ' %g', ext(i,1));
    fprintf(fid, ' %g', sca(i,1));
    fprintf(fid, ' %g', ext(i,1)-sca(i,1));
    fprintf(fid,'\n');
end

end

%%
%for num = 1:30;
    %A = load(strcat('prolate_spec_nm',num2str(num)));
    %plot(A(:,1),A(:,4))
    %pause
%end