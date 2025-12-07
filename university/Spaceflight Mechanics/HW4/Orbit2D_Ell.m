function[tabella_posPer]=Orbit2D_Ell(a,e,dt,anomalia_vera_0,modo,ind_fig,color)
%   ORBIT2D_ELL
%   disegna un orbita dati a, e, il tempo di volo e ni0
%   Trova gli assi perifocali e disegna il plot dell'orbita ellittica in 2D

if a<0
    error('a deve essere un valore in [0,+inf)');
end
if e<0 || e>=1
    error('e deve essere un valore in [0,1)');
end

% selezione stella/pianeta/satellite
muVect=[10^20,132712000000,22033,324860,398600,4890,...
    42828,126686535,37931188,5793940,6836530];
modi_disp={'cerchio','sole','mercurio','venere','terra','luna',...
    'marte','giove','saturno','urano','nettuno'};
confronto=strcmp(modo,modi_disp);
if sum(confronto)==0
    error('Il character array inserito non è valido');
end
mu=muVect(confronto==1);

moto_medio=sqrt(mu/a^3);
T_orb=2*pi/moto_medio;
if confronto(1)==1
    dt=T_orb;
end
eps=10^(-7);

E_iniziale=2*atan(sqrt((1-e)/(1+e))*tan(anomalia_vera_0/2));
M_iniziale=E_iniziale-e*sin(E_iniziale);

tempo=linspace(0,dt,1000);
anomalia_vera=zeros(1,length(tempo));
for it=1:length(tempo)
    t=tempo(it);
    M=moto_medio*t+M_iniziale;
    err=eps+1;
    E0=pi;
    while err>eps
            En=E0-(E0-e*sin(E0)-M)/(1-e*cos(E0));
            err=abs(En-E0);
            E0=En;
    end
    anomalia_vera(it)=2*atan(sqrt((1+e)/(1-e))*tan(En/2));
end
p=a*(1-e^2);
r=p./(1+e*cos(anomalia_vera));
xPer=r.*cos(anomalia_vera);
yPer=r.*sin(anomalia_vera);

tabella_posPer=[xPer;yPer];
%figure;
%hold on
%plot(xPer,yPer,'k');
%plot(0,0,'ko');
%plot(xPer(1),yPer(1),'ro');
%hold off
%daspect([1 1 1]);
%grid on
%ax=gca;
%ax.XAxis.Exponent=0;
%ax.YAxis.Exponent=0;