function modePower = radResRect_SingFreq180bypi(m_TE,n_TE,m_TM,n_TM,a,b,fc_TE,fc_TM,k)
eta = 377;
% l= 0.0106;
% d=0.005;
% l= 0.0306;
% d=0.2;
l = a/8;
d = b/2;
r_TE=0;
r_TM=0;
for n = 1:length(n_TE)
    D_TE_Sq = ((cos(n_TE(n)*l*180/b) - cos(k*l*180/pi))/((1 - (((n_TE(n)*pi)/(k*b))^2))*k))^2;
    g_TE_Sq = (m_TE(n)*pi/a)^2 + (n_TE(n)*pi/b)^2;
    beta_TE = sqrt(k^2 - g_TE_Sq);
    if n_TE(n)==0
        chi_n=1;
    else
        chi_n=2;
    end
    if m_TE(n)==0
        chi_m=1;
    else
        chi_m=2;
    end
    r_TE(n) = real((eta*(pi^2)*D_TE_Sq*chi_n*chi_m*((sin(m_TE(n)*d*180/a))^2)*...
        k*(m_TE(n)^2))/(2*beta_TE*b*g_TE_Sq*((sin(k*l*180/pi))^2)*(a^3)));
    if isnan(r_TE(n)) == 1
        r_TE(n) = 0;
    end
    
end
for n = 1:length(n_TM)
    D_TM_Sq = ((cos(n_TM(n)*l*180/b) - cos(k*l*180/pi))/((1 - (((n_TM(n)*pi)/(k*b))^2))*k))^2;
    g_TM_Sq = (m_TM(n)*pi/a)^2 + (n_TM(n)*pi/b)^2;
    beta_TM = sqrt(k^2 - g_TM_Sq);
    if n_TM(n)==0
        chi_n=1;
    else
        chi_n=2;
    end
    if m_TM(n)==0
        chi_m=1;
    else
        chi_m=2;
    end
    r_TM(n) = real((eta*(pi^2)*D_TM_Sq*chi_n*chi_m*((sin(m_TM(n)*d*180/a))^2)*...
        beta_TM*((n_TM(n))^2))/(2*beta_TM*a*g_TM_Sq*((sin(k*l*180/pi))^2)*(b^3)));
    if isnan(r_TM(n)) == 1
        r_TM(n) = 0;
    end
end

r_TE = r_TE';
r_TM = r_TM';

radRes = [r_TE; r_TM];
fc = [fc_TE; fc_TM];

counting = zeros(1,length(radRes));
sumRadRes = sum(radRes);

for i = 1:length(radRes)
    p(i) = 100*radRes(i)/sumRadRes;
    if p(i) >= 5
        counting(i) = i;
    end
end
stem(fc,p');
%modePower = [sum(radRes) counting];
modePower = [p; counting];
end

