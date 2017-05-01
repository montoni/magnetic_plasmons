clear all
close all
clc

for num = [1,5,10,15,20,25,30]
    A = load(strcat('coarse_bem_data/phenalene_nm',num2str(num)));
    clf
    plot(A(:,1),A(:,2),A(:,1),A(:,3))
    legend('NS' , 'NN')
    num
    pause
end