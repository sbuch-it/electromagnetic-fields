%% calcola i coefficienti di riflessione e trasmissione tra 2 mezzzi
clear all
close all

%% mezzo 1
epsr1=1;
mur1=1;

%% mezzo 2
epsr2=2.2;
mur2=1;

eta1=sqrt(mur1/epsr1); % impedenze normalizzate
eta2=sqrt(mur2/epsr2);

thetainc=linspace(0,pi/2,1001); % angolo di incidenza
thetatra=conj(asin(sqrt(epsr1*mur1)/sqrt(epsr2*mur2)*sin(thetainc))); % legge di Snell

Gammapar=(eta1*cos(thetainc)-eta2*cos(thetatra))./(eta1*cos(thetainc)+eta2*cos(thetatra));
Gammaprp=(eta2./cos(thetatra)-eta1./cos(thetainc))./(eta2./cos(thetatra)+eta1./cos(thetainc));

Taupar=2*eta2*cos(thetainc)./(eta1*cos(thetainc)+eta2*cos(thetatra));
Tauprp=2*eta2./cos(thetatra)./(eta2./cos(thetatra)+eta1./cos(thetainc));

%% plot
figure(1)
plot(thetainc,abs(Gammapar),thetainc,abs(Gammaprp),thetainc,abs(Taupar),thetainc,abs(Tauprp),'linewidth',2)
xlabel('\theta^{inc} (rad)'), ylabel('|\Gamma|')
legend('\Gamma_{||}','\Gamma_{\perp}','T_{||}','T_{\perp}')

figure(2)
plot(thetainc,real(Gammapar),thetainc,real(Gammaprp),thetainc,imag(Gammapar),thetainc,imag(Gammaprp),'linewidth',2)
xlabel('\theta^{inc} (rad)'), ylabel('\Gamma')
legend('\Ree\{\Gamma_{||}\}','\Ree\{\Gamma_{\perp}\}','\Imm\{\Gamma_{||}\}','\Imm\{\Gamma_{\perp}\}')

figure(3)
plot(thetainc,real(Taupar),thetainc,real(Tauprp),thetainc,imag(Taupar),thetainc,imag(Tauprp),'linewidth',2)
xlabel('\theta^{inc} (rad)'), ylabel('T')
legend('\Ree\{T_{||}\}','\Ree\{T_{\perp}\}','\Imm\{T_{||}\}','\Imm\{T_{\perp}\}')