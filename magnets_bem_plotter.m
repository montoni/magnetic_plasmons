clear all
close all
clc

for num = 1:10
    A = load(strcat('optical/twomer_optical_sep',num2str(num)));
    clf
    eV = 1240./(A(:,1))
    plot(eV,A(:,3),eV,A(:,5))
    legend('xy' , 'yz')
    num
    pause
end
