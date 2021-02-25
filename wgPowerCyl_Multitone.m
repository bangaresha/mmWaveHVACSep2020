function WGPower = wgPowerCyl_Multitone(mTE,nTE,mTM,nTM,radius,freq,betaTE,alphaTE,betaTM,alphaTM,k,WG_len,mu,epsilon)
Eio=10;
for fi = 1:length(freq)
    for ni = 1:length(nTE)
        attTE(fr,ni) = alphaTE(fi,ni)*WG_len;
        hsqrTE(fr,ni) = ((k(fi))^2) - betaTE(fi,ni)^2;
        hTE(fr,ni) = sqrt(hsqrTE(fr,ni));
        om = 2*pi*freq(fi);
        %D = Eio/((omega*mu/hsqrTM)*(nTE(fi,ni)*pi/b));

        funcErHtTE = @(r,theta) Eio.*(k(fr,ni))*om*mu*((mTE(fr,ni)).^2).*...
            ((((besselj(mTE(fr,ni),r.*hTE(fr,ni)))*sin(theta.*pi.*...
            (mTE(fr,ni))./180))./((hsqrTE(fr,ni)).*r)).^2).*exp((-2.*attTE(fr,ni)).*z);
        PTEend(fr,ni) = -0.5*integral2(funcErHtTE,0,radius,0,2*pi);

%         if nTE(ni)==0
%             chi_n=1;
%         else
%             chi_n=2;
%         end
%         if mTE(ni)==0
%             chi_m=1;
%         else
%             chi_m=2;
%         end
%         P_TE(ni) = (1/(2*eta*chi_n*chi_m))*((k*Eio)^2)*a*b*sqrt(1-((fc(ni)/freq)^2))/...
%             (((mTE(ni)*pi/a)^2)+((nTE(ni)*pi/b)^2));
    end
    for ni = 1:length(nTM)
        attTM(fr,ni) = alphaTM(ni)*WG_len;
        hsqrTM(fr,ni) = (k^2) - betaTM^2;
        hTM(fr,ni) = sqqrt(hsqrTM(fr,ni));
        %D = Eio/((betaTM/hsqr)*(mTM(ni)*pi/a));
        om = 2*pi*freq(fr);
        funcErHtTM = @(r,theta) ((k(fr))*om*epsilon*((mTM(fr,ni)).^2).*...
            ((((besselj(mTM(fr,ni),r.*hTM(fr,ni)))*sin(theta.*pi.*...
            (mTM(fr,ni))./180))./((hsqrTM(fr,ni)).*r)).^2).*Eio).*...
            exp((-2.*attTM(fr,ni)).*z);

        PTMend(fr,ni) = -1*0.5*integral2(funcErHtTM,0,radius,0,2*pi);
%         if nTM(ni)==0
%             chi_n=1;
%         else
%             chi_n=2;
%         end
%         if mTM(ni)==0
%             chi_m=1;
%         else
%             chi_m=2;
%         end
%         P_TM(ni) = (1/(2*eta*chi_n*chi_m))*((k*Eio)^2)*a*b*sqrt(1-((fc(ni)/freq)^2))/...
%             (((mTM(ni)*pi/a)^2)+((nTM(ni)*pi/b)^2));
    end
%WGPower_input = [P_TE P_TM];
WGPower_output = [PTEend PTMend];
alpha = [-1*attTE -1*attTM];
WGPower = [WGPower_input; WGPower_output;alpha];
end