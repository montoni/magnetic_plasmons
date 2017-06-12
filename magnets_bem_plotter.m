clear all
close all
clc

for num = [1,15,30]
    A = load(strcat('optical/twomer_optical_nm',num2str(num)));
    clf
    eV = 1240./(A(:,1))
    plot(eV,A(:,2)-A(:,3),eV,A(:,4)-A(:,5))
    legend('xy' , 'yz')
    num
    pause
end
