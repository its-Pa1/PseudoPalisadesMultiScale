clear all;
clc;

radius = 20;
sigma = 5;
x=0:0.1:500;
y=0:0.1:500;
centerx = 100;
centery = 150;
dx =  (x-centerx).^2;
dy =  (y-centery).^2;
distance = sqrt(dx+dy );

f_val = exp(-(dx+dy)/((sigma^2)));
if (distance<radius)
    ff=f_val;
else
    ff=0;
end
[x,y] = meshgrid(x,y);
surf(x,y,f_val);