clc
clearvars
cla
close all

%%
N = 161;
T = 100001;
L = 80;
x = linspace(-L/2,L/2,N);
t = linspace(0,L,T);
dx = L/(N-1);
dt = L/(T-1);
AA = 1;
w = 1;
k = 1;

h = dt/(dx^2);


%%
u = zeros(N,1);

%зависимость от времени
% u = zeros(N,T);
% for i = 1:T
% for j = 1:N
%    
%         u(j,i) = 0;% -(x(j)/10) ^ 2 + 16 * abs(cos(i/1000));
%   
% end
% end
% 
% for i = 1:T
% plot(u(:,i*8));
% pause(0.000001);
% end 

% кусочный потенциал
for j = 1:N
    if j >= N/2 && j <=  N
        u(j) = 5;
    else
        u(j) = 0;
    end
end


% многослойная структура
% for j = N/2-0.5:N-5
%     if mod(j,2) == 1
%         u(j) = 6;
%     else
%         u(j) = 0;
%     end
% end
% % 
% u(:) = 5.*exp( -((x-20).^6)./(4*1));
%%
A = zeros(N,N);
for i = 1:N
    for j = 1:N
        if i == j
            A(i,j) = (1 + 1i*h) - 1i * u(i) * dt;
        elseif j == i+1
            A(i,j) = (-1i / 2) * h;
        elseif j == i-1
            A(i,j) = (-1i / 2) * h;
        end


    end
end


%%

Psi = zeros(N,T, 'like', complex(0,0));

Psi(:,1) = -1.*exp(1i.*(x).*-10 ).*exp( -((x+25).^2)./(4*4));
% Psi(:,1) =  4 * exp( -((x).^2));
% Psi(:,1) =  (1.25.*exp(1i.*(x).*3 ) + 1.5.*exp(1i.*(x).*6)).*exp( -((x+25).^2)./(4*4));

%


%%
for g = 1:T-1
    Psi(:,g+1) = linsolve(A,Psi(:,g));
end


%%
 
figure(1)
mesh(abs(Psi(:,1:70000)).^2);
grid off

for i = 1:100000/16
    figure(2)
    plot(x,real(Psi(:,i*16)),x,u(:));
    pause(0.01);
end 


% под u зависящий от времени
% for i = 1:100000/16
%     figure(2)
%     plot(x,abs(Psi(:,i*16)).^2,x,u(:,i*16));
%     pause(0.01);
% end 
