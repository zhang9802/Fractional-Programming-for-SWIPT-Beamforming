%%  sum  rate maximization for SWIPT
% sum rate maximization via quardatic transform 
% reference: Active Reconfigurable Intelligent Surface: Fully-Connected or Sub-Connected
% maximize sum rate, subject to power constraint
clc,clear;
% close all;
s = rng(0);

Nt = 6;
K = 2;
H = rand(Nt,K)+1i*rand(Nt,K);
he =  rand(Nt,1)+1i*rand(Nt,1);
iter_num = 100;
W = zeros(Nt,K+1,iter_num+1);

val = zeros(1,iter_num+1);
mu = zeros(K,iter_num+1);
v = zeros(K,iter_num+1);
Pmax = 3;
sigma = 1;
delta = 1e-4;
% initialization
W(:,:,1) = rand(Nt,K+1)+1i*rand(Nt,K+1);
mu(:,1) = rand(K,1);
v(:,1) = rand(K,1);
eta = 0.8; % energy transform efficiency
Emax = 8;
for iter = 1:iter_num

% optimize mu and v
% mu_k = zeros(K,1);
% v_k  = zeros(K,1);

for kk=1:K
   rho = real(conj(v(kk,iter))*H(:,kk)'*W(:,kk,iter));
   mu(kk,iter+1)= rho/2*(rho+sqrt(rho^2+4));
   v(kk,iter+1) = sqrt(1+mu(kk,iter+1))*H(:,kk)'*W(:,kk,iter)/(norm(H(:,kk)'*W(:,:,iter))^2+sigma);
end
   
% optimize W
cvx_begin quiet
variable W0(Nt,K+1) complex
maximize (2*sqrt(1+mu(1,iter+1))*real(v(1,iter+1)'*H(:,1)'*W0(:,1))-abs(v(1,iter+1))^2*(square_pos(norm(H(:,1)'*W0))+sigma)...
    +2*sqrt(1+mu(2,iter+1))*real(v(2,iter+1)'*H(:,2)'*W0(:,2))-abs(v(2,iter+1))^2*(square_pos(norm(H(:,2)'*W0))+sigma))

subject to
norm(W0)<=Pmax;
eta*(norm(he'*W(:,:,iter))^2+sigma+real(trace(2*W(:,:,iter)'*(he*he')*(W0-W(:,:,iter)))))>=Emax
cvx_end
W(:,:,iter+1)=W0;
val(iter+1) = log2(1+abs(H(:,1)'*W0(:,1))^2/(norm(H(:,1)'*W0(:,2:3))^2+sigma)) + log2(1+abs(H(:,2)'*W0(:,2))^2/(norm(H(:,2)'*W0(:,[1,3]))^2+sigma));
disp(iter)
tf = strcmp(cvx_status,'Solved');
if tf~=1
    break;
end
if abs(val(iter+1) - val(iter))<=delta 
    break
end
end

figure
plot(val(2:iter+1),'b-','linewidth',1.5);
xlabel('iteration index');
ylabel('Objective value');
grid on
% %% exhaustive search method
% MentCarlo_num = 1e6;
% val1 = zeros(MentCarlo_num,1);
% for ii=1:MentCarlo_num
%     W1 = randn(Nt,K)+1i*randn(Nt,K);
%     W1 = W1/norm(W1)*Pmax*rand(1);
%     val1(ii) = log2(1+abs(H(:,1)'*W1(:,1))^2/(abs(H(:,1)'*W1(:,2))^2+sigma)) + log2(1+abs(H(:,2)'*W1(:,2))^2/(abs(H(:,2)'*W1(:,1))^2+sigma)); 
% end
% disp(max(val1));
