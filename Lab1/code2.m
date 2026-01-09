clc;
clear;
close all;

t=[0 0.15 0.31 0.5 0.6 0.75];
y=[1.0 1.004 1.031 1.117 1.223 1.422];

t_fine=linspace(0,0.75,200);

figure;
hold on;
grid on;

plot(t,y,'ro');

for n=1:4
    p=polyfit(t,y,n);
    y_out=polyval(p,t);
    E=sum((y-y_out).^2);
    fprintf("%d\n",n);
    disp(p);
    fprintf("%d\n",E);
    plot(t_fine,polyval(p,t_fine));
end

legend('Data','Degree 1','Degree 2','Degree 3','Degree 4','Location','northwest');