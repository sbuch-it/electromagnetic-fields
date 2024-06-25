%% esercizio 1
clear all
close all

f=3e9; % frequenza
omega=2*pi*f; % frequenza angolare
c=299792458; % velocità della luce
mu0=4e-7*pi; % permeabilità del vuoto
eps0=1/(c^2*mu0); % permittività del vuoto
eta0=sqrt(mu0/eps0); % impedenza dello spazio libero
k0=omega/c; % numero d'onda nello spazio libero
lambda0=c/f; % lunghezza d'onda nel vuoto

thetainc=pi/3; % angolo di incidenza

%% mezzo 1
epsr1=1; mur1=1; % aria
%  epsr1=2.2; mur1=1; %dielettrico
% epsr1=4*(1-1i*4e-4); mur1=1; % dielettrico con perdite

eta1=eta0*sqrt(mur1/epsr1); % impedenza caratteristica
k1=k0*(-1i)*sqrt(-mur1*epsr1); % numero d'onda
Zpar1=eta1*cos(thetainc);
Zprp1=eta1/cos(thetainc);

uinc=[sin(thetainc); 0; cos(thetainc)];
uincpar=[cos(thetainc); 0; -sin(thetainc)];
uprp=[0;1;0];
uref=[sin(thetainc); 0; -cos(thetainc)];
urefpar=[-cos(thetainc); 0; -sin(thetainc)];

%% mezzo 2
%  epsr2=1; mur2=1; % aria
% epsr2=2.2; mur2=1; % PTFE
sigma=5.8e7; epsr2=1-1i*sigma/(omega*eps0);  mur2=1;% conduttore
% sigma=5e-3; epsr2=15-1i*sigma/(omega*eps0); mur2=1;% suolo
% sigma=5e-3; epsr2=80-1i*sigma/(omega*eps0); mur2=1;% Sea Water
%epsr2=4-0.4i; mur2=1; % dielettrico con perdite
%fp=sqrt(6)*1e9; epsr2=1-(fp/f)^2; mur2=1; % plasma

eta2=eta0*sqrt(mur2/epsr2); % impedenza caratteristica
k2=k0*(-1i)*sqrt(-mur2*epsr2); % costante di propagazione
%delta=(1-1i)/k2;
thetatra=(asin(k1*sin(thetainc)/k2)); % angolo di trasmissione
Zpar2=eta2*cos(thetatra);
Zprp2=eta2/cos(thetatra);
utra=[sin(thetatra); 0; (cos(thetatra))];
utrapar=[cos(thetatra); 0; -sin(thetatra)];

%% campo incidente
% S0=100; % W/m^2
%E0=sqrt(2*eta1*S0);
E0=1; % V/m
E0par=E0;
E0prp=0;

% E0par=E0/sqrt(2);
% E0prp=-1i*E0/sqrt(2);

Gammapar=(Zpar1-Zpar2)/(Zpar1+Zpar2);
Gammaprp=(Zprp2-Zprp1)/(Zprp2+Zprp1);

Taupar=eta2/eta1*(1+Gammapar);
Tauprp=1+Gammaprp;

[X,Z]=meshgrid(linspace(-2*lambda0,2*lambda0,1001),linspace(-2*lambda0,2*lambda0,1001));
Propinc=exp(-1i*k1*(X*uinc(1)+Z*uinc(3))); Propinc(Z>=0)=0; % fattori di propagazione
Propref=exp(-1i*k1*(X*uref(1)+Z*uref(3))); Propref(Z>=0)=0;
Proptra=exp(-1i*k2*(X*utra(1)+Z*utra(3))); Proptra(Z<0)=0;

Uzn=ones(size(Z)); Uzn(Z>=0)=0;
Uzp=ones(size(Z)); Uzp(Z<0)=0;

for n=1:3
    Einc(:,:,n)=(E0par*uincpar(n)+E0prp*uprp(n))*Propinc;
    Hinc(:,:,n)=(E0par*uprp(n)-E0prp*uincpar(n))/eta1*Propinc;
    Erif(:,:,n)=(E0par*Gammapar*urefpar(n)+E0prp*Gammaprp*uprp(n))*Propref;
    Hrif(:,:,n)=(E0par*Gammapar*uprp(n)-E0prp*Gammaprp*urefpar(n))/eta1*Propref;
    Etra(:,:,n)=(E0par*Taupar*utrapar(n)+E0prp*Tauprp*uprp(n))*Proptra;
    Htra(:,:,n)=(E0par*Taupar*uprp(n)-E0prp*Tauprp*utrapar(n))/eta2*Proptra;
end

E=Einc+Erif+Etra;
H=Hinc+Hrif+Htra;
S=1/2*cross(E,conj(H),3);
Eabs=sqrt(sum(abs(E).^2,3));
Habs=sqrt(sum(abs(H).^2,3));

Nt=36;
cax=[-1 1]*max(max(Eabs));
figure(1)
pcolor(X,Z,real(E(:,:,1))); shading flat;axis equal tight; colormap jet; caxis(cax); colorbar, title('e_x(t) (V/m)'), xlabel('x (m)'), ylabel('z (m)')
% for nt=1:Nt
%    pcolor(X,Z,real(E(:,:,1)*exp(1i*2*pi*(nt-1)/Nt))); shading flat;axis equal tight; colormap jet; caxis(cax); colorbar, title('e_x(t) (V/m)'), xlabel('x (m)'), ylabel('z (m)')
%     Mov(nt)=getframe;
% end
% movie(Mov,100,2*Nt)

figure(2)
pcolor(X,Z,real(E(:,:,2))); shading flat;axis equal tight; colormap jet; caxis(cax); colorbar, title('e_y(t) (V/m)'), xlabel('x (m)'), ylabel('z (m)')
figure(3)
pcolor(X,Z,real(E(:,:,3))); shading flat;axis equal tight; colormap jet; caxis(cax); colorbar, title('e_z(t) (V/m)'), xlabel('x (m)'), ylabel('z (m)')

cax=[-1 1]*max(max(Habs));
figure(11)
pcolor(X,Z,real(H(:,:,1))); shading flat;axis equal tight; colormap jet; caxis(cax); colorbar, title('h_x(t) (A/m)'), xlabel('x (m)'), ylabel('z (m)')
figure(12)
pcolor(X,Z,real(H(:,:,2))); shading flat;axis equal tight; colormap jet; caxis(cax); colorbar, title('h_y(t) (A/m)'), xlabel('x (m)'), ylabel('z (m)')
figure(13)
pcolor(X,Z,real(H(:,:,3))); shading flat;axis equal tight; colormap jet; caxis(cax); colorbar, title('h_z(t) (A/m)'), xlabel('x (m)'), ylabel('z (m)')

figure(20)
cax=[0 1]*max(max(Eabs));
pcolor(X,Z,Eabs); shading flat;axis equal tight; colormap jet; caxis(cax); colorbar, title('|E| (V/m)'), xlabel('x (m)'), ylabel('z (m)')
figure(21)
cax=[0 1]*max(max(Habs));
pcolor(X,Z,Habs); shading flat;axis equal tight; colormap jet; caxis(cax); colorbar, title('|H| (A/m)'), xlabel('x (m)'), ylabel('z (m)')

cax=[-1 1]*max(max(sqrt(dot(S,S,3))));
figure(30)
pcolor(X,Z,real(S(:,:,1))); shading flat;axis equal tight; colormap jet; caxis(cax); colorbar, title('\Ree S_x (W/m^2)'), xlabel('x (m)'), ylabel('z (m)')
figure(31)
pcolor(X,Z,real(S(:,:,2))); shading flat;axis equal tight; colormap jet; caxis(cax); colorbar, title('\Ree S_y (W/m^2)'), xlabel('x (m)'), ylabel('z (m)')
figure(32)
pcolor(X,Z,real(S(:,:,3))); shading flat;axis equal tight; colormap jet; caxis(cax); colorbar, title('\Ree S_z (W/m^2)'), xlabel('x (m)'), ylabel('z (m)')
figure(40)
pcolor(X,Z,imag(S(:,:,1))); shading flat;axis equal tight; colormap jet; caxis(cax); colorbar, title('\Imm S_x (VAR/m^2)'), xlabel('x (m)'), ylabel('z (m)')
figure(41)
pcolor(X,Z,imag(S(:,:,2))); shading flat;axis equal tight; colormap jet; caxis(cax); colorbar, title('\Imm S_y (VAR/m^2)'), xlabel('x (m)'), ylabel('z (m)')
figure(42)
pcolor(X,Z,imag(S(:,:,3))); shading flat;axis equal tight; colormap jet; caxis(cax); colorbar, title('\Imm S_z (VAR/m^2)'), xlabel('x (m)'), ylabel('z (m)')

we=eps0/2*real(epsr1*Uzn+epsr2*Uzp).*Eabs.^2;
wm=mu0/2*real(mur1*Uzn+mur2*Uzp).*Habs.^2;

cax=[0,max(max(max(we)),max(max(wm)))];
figure(50)
pcolor(X,Z,we); shading flat;axis equal tight; colormap jet; caxis(cax); colorbar, title('w_e=1/2\epsilon|E|^2 (J/m^3)'), xlabel('x (m)'), ylabel('z (m)')
figure(51)
pcolor(X,Z,wm); shading flat;axis equal tight; colormap jet; caxis(cax); colorbar, title('w_m=1/2\mu|H|^2 (J/m^3)'), xlabel('x (m)'), ylabel('z (m)')
