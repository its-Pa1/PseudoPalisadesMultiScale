% This matlab scripts read the DW data from ncdf file and save DWs for a
% given domain and plots fractional anisotropy of DW
clear all;
clc;
close all;
%%
ncdisp('data.ncdf') %it will display all the variables stored and their info
W_raw  = ncread('data.ncdf','DW'); % reading th raw data
%% memory allocations
dim = 2; % 2D or 3D
tensor = zeros(dim,dim); % temp variable store Dw
W_2D = zeros(dim,dim,size(W_raw,3)*size(W_raw,4)*size(W_raw,5)); %stores Dw at each points in the domain 
%% main loop to store data
% here we store the 2D data at each point of the domain. Unlike the other
% script, here we save the data using coulmn-major linear indexing and call
% from that setting only
for m = 1:size(W_raw,5)
    for l = 1:size(W_raw,4)
        for k = 1:size(W_raw,3)
            for j = 1:dim
                for i = 1:dim
                    index = sub2ind(size(W_raw),i,j,k,l,m);
                    tensor(i,j) = W_raw(index);
                end
            end
            mapped_index = sub2ind([size(W_raw,3),size(W_raw,4),size(W_raw,5)],k,l,m);
            
            W_2D(:,:,mapped_index) = tensor;
        end
    end
end
%%
% x = 100;
% y = 150;
% slice = 100;
% 
% index_x = floor(x/2);
% index_y = floor(y/2);
% 
% temp2 = sub2ind([size(W,3),size(W,4),size(W,5)],slice+1,index_y,index_x);
% test2 = W_2D(:,:,temp2);
%%
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
            temp = sub2ind([size(W_raw,3),size(W_raw,4),size(W_raw,5)],index_z+1,index_y,index_x);
        test(:,:,i,j) = W_2D(:,:,temp); % stores Dw at the domain points
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

