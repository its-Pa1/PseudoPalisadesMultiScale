% This matlab scripts read the DW data from ncdf file and save DWs for a
% given domain and plots fractional anisotropy of DW
clear all;
clc;
close all;
tic
%% data reading
ncdisp('data.ncdf') % it will display all the variables stored and their info
W_raw = ncread('data.ncdf','DW'); % reading th raw data
%% Transform data and memory allocations
dim = 2; % 2D or 3D
count = 1;
% As it can be seen from ncdisp results that W_raw has dimension 3x3xZxYxX
% so we transform it to XxYxZx3x3 using permute function
W_permute = permute(W_raw,[5,4,3,2,1]); 
% memory allocations
index = zeros(1,size(W_permute,1)*size(W_permute,2)*size(W_permute,3)*size(W_permute,4)*size(W_permute,5));
W_D = zeros(dim,dim,size(W_permute,1)*size(W_permute,2)*size(W_permute,3));
tensor = zeros(dim,dim); % temp variable store Dw
%% Save as single index
% as matlab store a matrix/tensor using column-major linear indexing, which
% creates truoble sometime, we convert it here to row-major indexing as C++

for i = 1:size(W_permute,1)
    for j = 1:size(W_permute,2)
        for k = 1:size(W_permute,3)
            for l = 1:size(W_permute,4)
                for m = 1:size(W_permute,5)
                    index(count) = sub2ind(size(W_permute),i,j,k,l,m);
                    count = count+1;
                end
            end
           
        end
    end   
end
% this variable W_lin stores the entire elements of Dw as a vector/array
% using 
W_lin = W_permute(index);

% now same as C++ code, we use the linear indexing(now row major) to store
% the Dw as 2D or 3D matrix in W_D
for i = 1:size(W_permute,1)
    for j = 1:size(W_permute,2)
        for k = 1:size(W_permute,3)
            for ii = 1:dim
                for jj = 1:dim
                    index2 = jj+size(W_permute,5)*(ii-1+size(W_permute,4)*(k-1+size(W_permute,3)*...
                        (j-1+size(W_permute,2)*(i-1))));
                    tensor(ii,jj) = W_lin(index2);
                end
            end
            mapped_index = k + size(W_permute,3)*((j-1)+size(W_permute,2)*(i-1));
            W_D(:,:,mapped_index)= tensor;
        end
    end
end
toc
%% To access Dw at a particular point
% x = 120;
% y = 180;
% z = 200;
% 
% index_x = floor(x/2);
% index_y = floor(y/2);
% index_z = floor(z/2);
% temp = 1+(index_z + size(W_permute,3)*((index_y)+size(W_permute,2)*(index_x)));
% % fprintf('Water Diff tensor at x = %d, y =%d, z = %d is',x,y,z);
% DW = W_D(:,:,temp)
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
        temp = 1+(index_z + size(W_permute,3)*((index_y)+size(W_permute,2)*(index_x)));
        
        if (temp>size(W_D,3))
            test(:,:,i,j) = 0;
        else
        test(:,:,i,j) = W_D(:,:,temp); % stores Dw at the domain points
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