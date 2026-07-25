%% optimize x/(x^2+1) s.t. x>=0
% optimal val is 0.5
clc,clear;
close all

iter_num = 10;
x = zeros(iter_num+1,1);
x(1) = 0.1;
y = zeros(iter_num+1,1);
y(1) = 0.1;
val = zeros(iter_num+1,1);

for iter=1:iter_num
    
    y(iter+1) = x(iter)^(1/2)/(x(iter)^2+1);
    x(iter+1) = (2*y(iter+1))^(-2/3);
    val(iter+1) = x(iter+1)/(x(iter+1)^2+1);
    
    
    
end
figure
plot(val(2:iter_num+1),'b-o','linewidth',1.5);
xlabel('iteration index');
ylabel('Objective value');
grid on
ylim([min(val(2:iter_num+1)),0.501]);
