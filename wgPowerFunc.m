function WGPower = wgPowerFunc(n_TE,m_TE,n_TM,m_TM,a,b,f,fc)
eta = 377;
pi=3.14;
lambda = 3*10^8/f;
k=2*pi/lambda;
% x=a/2;
% y=b/2;
WG_len = 0.5;
z = WG_len;
mu= 4*pi*10^-7;
sigma = 1.4*10^6;
epsilon = 8.85*(10^-12);
omega = 2*pi*f;
Rs = sqrt(omega*mu/(2*sigma));
Eio=10;
for ni = 1:length(n_TE)
    %unit of beta = radians/meter
    beta = sqrt((k^2) - (m_TE(ni)*pi/a)^2 - (n_TE(ni)*pi/b)^2);
    %unit of alpha = nepers/meter
    alpha_TE(ni) = ((2*Rs*k/(eta*b*beta))*(((1+(b/a))*((fc(ni)/f)^2)) +...
        (b/a)*(1 - ((fc(ni)/f)^2))*(((n_TE(ni)*n_TE(ni)*a*b) +...
        (m_TE(ni)*m_TE(ni)*a*a))/((n_TE(ni)*n_TE(ni)*b*b) + (m_TE(ni)*m_TE(ni)*a*a)))));
    att_TE(ni) = 8.686*alpha_TE(ni)*WG_len;
    hsqr = (k^2) - beta^2;
    D = Eio/((omega*mu/hsqr)*(n_TE(ni)*pi/b));
%    Hz = D * cos(m_TE(ni)*x*180/a)*cos(n_TE(ni)*y*180/b)*exp((-alpha_TE(ni)-1i*beta)*z);
%     funcEx = D * (1i*omega*mu/hsqr)*(n_TE(ni)*pi/b)* cos(m_TE(ni)*x*180/a)*...
%         sin(n_TE(ni)*y*180/b)*exp((-alpha_TE(ni)-1i*beta)*z);
%     %Ex = Eio * exp((-alpha_TE(ni)-1i*beta)*z);
%     funcEy = D * (1i*omega*mu/hsqr)*(m_TE(ni)*pi/a)* sin(m_TE(ni)*x*180/a)*...
%         cos(n_TE(ni)*y*180/b)*exp((-alpha_TE(ni)-1i*beta)*z);
%     funcHx = D * (1i*beta/hsqr)*(m_TE(ni)*pi/a)* sin(m_TE(ni)*x*180/a)*...
%         cos(n_TE(ni)*y*180/b)*exp((-alpha_TE(ni)-1i*beta)*z);
%     funcHy = D * (1i*beta/hsqr)*(n_TE(ni)*pi/b)* cos(m_TE(ni)*x*180/a)*...
%         sin(n_TE(ni)*y*180/b)*exp((-alpha_TE(ni)-1i*beta)*z);
    
    funcExHy = @(x,y) real((1i*omega*mu/hsqr).*(1i.*beta./hsqr)).*...
        ((D*(n_TE(ni).*pi./b).* cos(m_TE(ni).*x.*180./a).*...
        sin(n_TE(ni).*y.*180./b)).^2).*exp((-2.*alpha_TE(ni)).*z);
    
%     funcEyHx = @(x,y) (1i*omega*mu/hsqr)*(1i*beta/hsqr)*...
%         ((D*(m_TE(ni)*pi/a)* sin(m_TE(ni)*x*180/a)*...
%         cos(n_TE(ni)*y*180/b))^2)*exp((-2*alpha_TE(ni))*z);
    
    %fun = @(x,y) 0.5*real(funcExHy - funcEyHx);
    P_TE_end(ni) = -1*0.5*integral2(funcExHy,0,a,0,b);
    
    if n_TE(ni)==0
        chi_n=1;
    else
        chi_n=2;
    end
    if m_TE(ni)==0
        chi_m=1;
    else
        chi_m=2;
    end
    P_TE(ni) = (1/(2*eta*chi_n*chi_m))*((k*Eio)^2)*a*b*sqrt(1-((fc(ni)/f)^2))/...
        (((m_TE(ni)*pi/a)^2)+((n_TE(ni)*pi/b)^2));
    
    
%     P_TE_end(ni) = (1/(2*eta*chi_n*chi_m))*abs(Ex*conj(Ex))*...
%         a*b*sqrt(1-((fc(ni)/f)^2))/sqrt(((m_TE(ni)*pi/a)^2)+((n_TE(ni)*pi/b)^2));
end
for ni = 1:length(n_TM)
    beta = sqrt((k^2) - (m_TM(ni)*pi/a)^2 - (n_TM(ni)*pi/b)^2);
    alpha_TM(ni) = ((2*Rs*k/(eta*b*beta))*(((n_TM(ni)*n_TM(ni)*(b^3)) +...
        (m_TM(ni)*m_TM(ni)*(a^3)))/((n_TM(ni)*n_TM(ni)*b*b*a) + (m_TM(ni)*m_TM(ni)*(a^3)))));
    att_TM(ni) = 8.686*alpha_TM(ni)*WG_len;
    hsqr = (k^2) - beta^2;
    D = Eio/((beta/hsqr)*(m_TM(ni)*pi/a));
%     Ez = Eio * sin(m_TM(ni)*x*180/a)*sin(n_TM(ni)*y*180/b)*exp((-alpha_TM(ni)-1i*beta)*z);
%     Ex = Eio * (-1i*beta/hsqr)*(m_TM(ni)*pi/a)* cos(m_TM(ni)*x*180/a)* ...
%         sin(n_TM(ni)*y*180/b)*exp((-alpha-1i*beta)*z);
%     Ex = Eio * exp((-alpha_TM(ni)-1i*beta)*z);
%     Ey = Eio * (-1i*beta/hsqr)*(n_TM(ni)*pi/b)* sin(m_TM(ni)*x*180/a)* ...
%         cos(n_TM(ni)*y*180/b)*exp((-alpha_TM(ni)-1i*beta)*z);
%     Hx = Eio * (1i*omega*epsilon/hsqr)*(n_TM(ni)*pi/b)* sin(m_TM(ni)*x*180/a)* ...
%         cos(n_TM(ni)*y*180/b)*exp((-alpha_TM(ni)-1i*beta)*z);
%     Hy = Eio * (-1i*omega*epsilon/hsqr)*(m_TM(ni)*pi/a)* cos(m_TM(ni)*x*180/a)* ... 
%         sin(n_TM(ni)*y*180/b)*exp((-alpha_TM(ni)-1i*beta)*z);
    funcExHy = @(x,y) real((-1i*omega*epsilon/hsqr)*(-1i*beta/hsqr)).*...
        ((D*(m_TM(ni)*pi/a).* cos(m_TM(ni).*x.*180./a).*...
        sin(n_TM(ni).*y.*180./b)).^2).*exp((-2*alpha_TM(ni))*z);
    
    P_TM_end(ni) = -1*0.5*integral2(funcExHy,0,a,0,b);
    if n_TM(ni)==0
        chi_n=1;
    else
        chi_n=2;
    end
    if m_TM(ni)==0
        chi_m=1;
    else
        chi_m=2;
    end
    P_TM(ni) = (1/(2*eta*chi_n*chi_m))*((k*Eio)^2)*a*b*sqrt(1-((fc(ni)/f)^2))/...
        (((m_TM(ni)*pi/a)^2)+((n_TM(ni)*pi/b)^2));
%     P_TM_end(ni) = (1/(2*eta*chi_n*chi_m))*abs(Ex*conj(Hy) - Ey*conj(Hx))*...
%         sqrt(1-((fc(ni)/f)^2))/sqrt(((m_TM(ni)*pi/a)^2)+((n_TM(ni)*pi/b)^2));
%     P_TM_end(ni) = (1/(2*eta*chi_n*chi_m))*abs(Ex*conj(Ex))*...
%         a*b*sqrt(1-((fc(ni)/f)^2))/sqrt(((m_TM(ni)*pi/a)^2)+((n_TM(ni)*pi/b)^2));
end
WGPower_input = [P_TE P_TM];
WGPower_output = [P_TE_end P_TM_end];
alpha = [-1*att_TE -1*att_TM];
WGPower = [WGPower_input; WGPower_output;alpha];
end