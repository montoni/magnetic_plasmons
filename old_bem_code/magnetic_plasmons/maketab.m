clear, clc
 
h = 4.135e-15;
c = 2.99e8;
 
A = load('Ag_dr.dat');
%B = load('Al_Im.txt');
 
%ref = [h*c./(1e-6*A(:,1)),A(:,2),B(:,2)];
%ref = [A(:,1),sqrt((sqrt(A(:,2).^2+A(:,3).^2)+A(:,2))/2)...
 %   ,sqrt((sqrt(A(:,2).^2+A(:,3).^2)-A(:,2))/2)];
ref2 = [A(:,1),sqrt((sqrt(A(:,2).^2+A(:,3).^2)+A(:,2))/2)...
    ,sqrt((sqrt(A(:,2).^2+A(:,3).^2)-A(:,2))/2)];
%  
% for j = 1:length(ref)
%     ref2(length(ref)+1-j,:) = ref(j,:);
% end
%  
plot(ref2(:,1),ref2(:,2),ref2(:,1),ref2(:,3))

fid = fopen(['Ag_Drude_b.dat'],'wt');
for j = 1:length(ref2)
          fprintf(fid,' %g', ref2(j,1));
          fprintf(fid,' %g', ref2(j,2));
          fprintf(fid,' %g', ref2(j,3));
          fprintf(fid, '\n');
end

fclose(fid)