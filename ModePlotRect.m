a=2;
b=a/2;
freq = 2.5E9;
c = 3E8;
m=2;
n=2;
% fc = c*100*sqrt((m/a).^2 + (n/b).^2)/2;
% fc = c*100*sqrt((m/a).^2 + (n/b).^2)/2;
epsilon = 8.8540e-12;           % Permittivity constant
epsilon_r = 1;                  % Relative Permittivity constant
mu1 = 4*pi*1E-7;
mu1_r = 1;
omega = 2*pi*freq;                 % Frequency of operation in rad/s
M=400;
k = omega/c;
g = sqrt((m*pi/a).^2 + (n*pi/b).^2);
beta = sqrt(k.^2 - g.^2);
WGlen = 0;
x = linspace(0,a,M);
y = linspace(0,b,M);
[x,y] = meshgrid(x,y);
Ex = cos(m*pi.*x./a).*sin(n*pi.*y./b).*exp(-1i*beta*WGlen);
Ey = -sin(m*pi.*x./a).*cos(n*pi.*y./b).*exp(-1i*beta*WGlen);
Hx = sin(m*pi.*x./a).*cos(n*pi.*y./b).*exp(-1i*beta*WGlen);
Hy = cos(m*pi.*x./a).*sin(n*pi.*y./b).*exp(-1i*beta*WGlen);
Hz = -cos(m*pi.*x./a).*cos(n*pi.*y./b).*exp(-1i*beta*WGlen);

% eta = 377;
% sigma = 1E6;
% Rs = sqrt((2*pi*freq*mu)./(2*sigma));
% if m ==0
%     chim = 1;
% else
%     chim = 2;
% end
% fcbyf = (fc/freq);
% fcbyfSq = fcbyf*fcbyf;
% alpha = (2*Rs*k/(eta*b*beta))*(((1+(b/a))*fcbyfSq) + (b/a)*((chim/2) - fcbyfSq)*...
% (((n*n*a*b) + (m*m*a*a))/((n*n*b*b) + (m*m*a*a))));
% gamma = alpha +1i*beta;

figure();
quiver(x,y,real(Ex),real(Ey));
title(['Plot of front view for TE_',num2str(m),'_',num2str(n),' E-Field']);
legend('E-Field');
xlabel('x-dimension 0 to a');
ylabel('y-dimension 0 to b=a/2');
figure();
quiver(x,y,real(Hx),real(Hy));
title(['Plot of front view for TE_',num2str(m),'_',num2str(n),' H-Field']);
legend('H-Field');
xlabel('x-dimension 0 to a');
ylabel('y-dimension 0 to b=a/2');
figure();
quiver(x,y,real(Ex),real(Ey));
hold on
quiver(x,y,real(Hx),real(Hy));
grid on
title(['Plot of front view for TE_',num2str(m),'_',num2str(n)]);
legend('E-Field','H-Field');
xlabel('x-dimension 0 to a');
ylabel('y-dimension 0 to b=a/2');

