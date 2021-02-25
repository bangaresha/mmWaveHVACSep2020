function [radresTE, radresTM, gammaTE, gammaTM] = radResCyl_singletone(m_TE,n_TE,m_TM,n_TM,radius,freq,c,k,R,eta,l,coWnTE,coWnTM)

r_TE=0;
r_TM=0;

for n = 1:length(n_TE)
    gTESq(n) = (coWnTE(n))^2;
    betaTE(n) = sqrt(k^2 - gTESq(n));
    alphaTempTE(n) = (R*k)/(eta*radius*betaTE(n));
    alphaTE(n)=8.85*(alphaTempTE(n)*(((m_TE(n)^2)/...
        ((coWnTE(n)^2)- (m_TE(n)^2)))+(((c*coWnTE(n))/(2*pi*freq))^2)));
    gammaTE(n)=alphaTE(n)+1i*betaTE(n);
    chrimpTE(n) = eta*k*(m_TE(n)^2)/(betaTE(n)*...
        pi*((sin(k*l*pi/180))^2));
    radTEfunc=@(zeta) ((besselj(m_TE(n),radius.*zeta.*coWnTE(n)).*...
        sin((k.*radius.*((l/radius)-1+zeta))*180/pi))./zeta);
    radTEnum(n) = (integral(radTEfunc,1-(l/radius),1)).^2;
    radTEden(n) = ((besselj(m_TE(n),radius*coWnTE(n)))^2)*...
        ((radius*coWnTE(n))^2 - (m_TE(n))^2);
    radresTE(n)= (chrimpTE(n)*radTEnum(n))/radTEden(n);
end
for n = 1:length(n_TM)
    gTMSq = (coWnTM(n))^2;
    betaTM(n) = sqrt(k^2 - gTMSq);
    alphaTM(n) = 8.85*((R*k)/(eta*radius*betaTM(n)));
    gammaTM(n)=alphaTM(n)+1i*betaTM(n);
    chrimpTM(n) = eta*betaTM(n)/(pi*((sin(l*180*(k)/pi))^2)*...
        k*((besselj(m_TM(n),coWnTM(n)))^2));
    radTMfunc=@(zeta) (((m_TM(n)./(zeta.*coWnTM(n))).*...
        besselj(m_TM(n),zeta.*coWnTM(n))) - ...
        besselj(m_TM(n)+1,zeta.*coWnTM(n))).*m_TM(n).*...
        sin((k.*radius.*((l./radius)-1+zeta))*180/pi);
    radTMnum(n) =(integral(radTMfunc,1-(l./radius),1)).^2;
    radresTM(n)= chrimpTM(n)*radTMnum(n);
end

r_TE = radresTE';
r_TM = radresTM';
end

