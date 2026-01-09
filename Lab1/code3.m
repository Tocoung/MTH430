clc;
clear;
close all;

t=[1 1.1 1.3 1.5 1.9 2.1]
y=[1.84 1.96 2.21 2.45 2.94 3.18]

t_fine=linspace(1,2.1,200);

figure;
hold on;
grid on;

plot(t,y,'ro');

for n=1:4
    p=polyfit(t,y,n);
    y_out=polyval(p,t);
    e=sum((y-y_out).^2);
    fprintf("%d\n",n);
    disp(p);
    disp(e);
    plot(t_fine,polyval(p,t_fine));
end

legend('Data','Degree1','Degree2','Degree3','Degree4','Location','northwest')