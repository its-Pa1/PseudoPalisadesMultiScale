% This matlab scripts read the DW data from ncdf file and save DWs for a
% given domain and plots fractional anisotropy of DW
clear all;
clc;
close all;
%%
ncdisp('data.ncdf') %it will display all the variables stored and their info
W_raw  = ncread('data.ncdf','DW'); % reading th raw data
%% Access the data on a domain and store it
% if domain range is outside the domain of Dw, Dw is zero

x = 0:.5:250; % x-domain
y = 0:.5:300; % y
index_z = 100;% fixed z i.e a slice
test = zeros(2,2,length(x),length(y));% stores Dw at each point of a domain
for i = 1:length(x)
    for j = 1:length(y)
        index_x = floor(x(i)/2);
        index_y = floor(y(j)/2);
        % if the called points are out of range set it zero
        if (index_x>size(W_raw,5) || index_y>size(W_raw,4) || index_z>size(W_raw,3)|| index_x<=0||index_y<=0||index_z<=0)
            
            test(:,:,i,j) = 0;
        else
            temp = W_raw(:,:,index_z+1,index_y,index_x); %access the value at x,y,z
            tensor = [temp(1,1),temp(1,2);temp(2,1),temp(2,2)];% trim the last row and column for 2D
            test(:,:,i,j) = tensor; % stores Dw at the domain points
        end
        
    end
end
%%  FA calculation
FA = zeros(size(test,3), size(test,4)); % stoes FA
for i = 1:size(test,3)
    for j= 1:size(test,4)
        lambda2 = eig(test(:,:,i,j));
        if (lambda2(1)==0 && lambda2(2)==0) % if no Dw then FA is 0
            FA(i,j) = 0;
        else
            FA(i,j) = abs(lambda2(1)-lambda2(2))/sqrt((lambda2(1)^2)+ (lambda2(2)^2));
        end
    end
end
%% Plot
figure(1)
surf(x,y,FA','LineStyle','none');
view(0,90);
colorbar
colormap('jet')
shading interp
title('FA of DW')
xlabel('X');
ylabel('Y')