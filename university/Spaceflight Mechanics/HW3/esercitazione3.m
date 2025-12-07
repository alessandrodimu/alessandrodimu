% ESERCITAZIONE 3
% PARTE 1

% dati sistema Terra-Luna
mu_terra=398600;
mu_luna=4890;
raggio_terra=6371;
raggio_luna=1738;
D_terra_luna=384400;
T_luna=2*pi*sqrt(D_terra_luna^3/mu_terra);
n_luna=2*pi/T_luna;
v_luna=n_luna*D_terra_luna;
raggio_soi_luna=66300;

% parametri dell'esercitazione parte 1
r_ini=6700;

% caso hohmann (caso limite)
a_hohmann=(r_ini+D_terra_luna)/2;
e_hohmann=(D_terra_luna-r_ini)/(D_terra_luna+r_ini);
v_hohmann=sqrt(mu_terra/a_hohmann*(1+e_hohmann)/(1-e_hohmann));

% r_0=6700, v_0=varia, phi_0=0
v_ini=linspace(v_hohmann,12,1000);
a_t=mu_terra./(2*mu_terra/r_ini-v_ini.^2);           % tramite energia
e_t=sqrt(1-(r_ini*v_ini).^2./(mu_terra*a_t));        % tramite momento ang
p_t=a_t.*(1-e_t.^2);
ni_fin=real(acos((p_t-D_terra_luna)./(D_terra_luna*e_t)));
ni_ini=0;
deltaNi=ni_fin-ni_ini;
deltaT=zeros(1,length(v_ini));

for i=1:length(v_ini)
    if e_t(i)<1
        Ef=2*atan(sqrt((1-e_t(i))/(1+e_t(i)))*tan(deltaNi(i)/2));
        deltaT(i)=(Ef-e_t(i)*sin(Ef))*sqrt(a_t(i)^3/mu_terra);
    elseif e_t(i)>1
        Hf=2*atanh(sqrt((e_t(i)-1)/(e_t(i)+1))*tan(deltaNi(i)/2));
        deltaT(i)=(e_t(i)*sinh(Hf)-Hf)*sqrt(-a_t(i)^3/mu_terra);
    end
end
rapportoNiT=diff(deltaNi)./diff(deltaT);
rapportoNiT(end+1)=rapportoNiT(end);
theta_luna_ini=deltaNi-n_luna*deltaT;

% disegno i grafici
figure(1)
plot(v_ini,e_t);
xlabel('v_0 (km/s)');
ylabel('e');
title('Eccentricita al variare della velocita iniziale');
grid on

figure(2)
plot(v_ini,a_t);
xlabel('v_0 (km/s)');
ylabel('a (km)');
title('Semiasse maggiore al variare della velocita iniziale');
grid on

figure(3)
deltaNi=deltaNi*180/pi;
plot(v_ini,deltaNi);
xlabel('v_0 (km/s)');
ylabel('angolo percorso dalla sonda (rad)');
title('Angolo percorso dalla sonda al variare della velocita iniziale');
grid on

figure(4)
plot(v_ini,deltaT/3600);
xlabel('v_0 (km/s)');
ylabel('T_(volo) (ore)');
title('Tempo di volo al variare della velocita iniziale');
grid on

% controllare BENE questa figura
figure(5)
plot(v_ini,rapportoNiT);
xlabel('v_0 (km/s)');
ylabel('deltaNi/deltaT (rad/s)');
title('deltaNi/deltaT al variare della velocita iniziale');
grid on

figure(6)
plot(v_ini,theta_luna_ini*180/pi);
xlabel('v_0 (km/s)');
ylabel('theta luna_0 (gradi)');
title('Angolo sonda/Luna alla partenza al variare della velocita iniziale');
grid on


% ESERCITAZIONE 3
% PARTE 2

% parametri dell'esercitazione parte 2
r_periselenio=2000;

% usiamo il metodo delle coniche raccordate,
% vediamo gli istanti del volo
%   0: partenza dall'orbita di parcheggio, sistema geocentrico
%   1: entrata nella soi lunare, sistema geocentrico
%   2: entrata nella soi lunare, sistema selenocentrico
%   3: periselenio, sistema selenocentrico
%   4: uscita dalla soi lunare, sistema selenocentrico
%   5: uscita dalla soi lunare, sistema geocentrico

% posizione, velocit� e angolo di volo nel punto 0
r_0=r_ini;
velocita_0=10.87;
phi_0=0;
energia_01=velocita_0^2/2-mu_terra/r_0;
momento_01=velocita_0*r_0*cos(phi_0);
a_01=-mu_terra/(2*energia_01);
p_01=momento_01^2/mu_terra;
e_01=sqrt(1-p_01/a_01);
ni_0=0;
E_0=0;
M_0=0;

% faccio ipotesi su lambda_1
lambda_1=linspace(0,pi,10000);

% posizione, velocit� e angolo di volo nel punto 1
r_1=sqrt(D_terra_luna^2+raggio_soi_luna^2+...
    2*D_terra_luna.*raggio_soi_luna.*cos(lambda_1));
velocita_1=sqrt(2*energia_01+(2*mu_terra)./r_1);
phi_1=acos(momento_01./(r_1.*velocita_1));
ni_1=acos((p_01-r_1)./(r_1*e_01));
E_1=2*atan(sqrt((1-e_01)/(1+e_01)).*tan(ni_1/2));
M_1=E_1-e_01.*sin(E_1);
gamma_1=-asin(raggio_soi_luna./r_1.*sin(lambda_1));
tempo_volo_01=(M_1-M_0).*sqrt(a_01^3/mu_terra);

% posizione, velocit� e angolo di volo nel punto 2
r_2=raggio_soi_luna;
velocita_2=sqrt(velocita_1.^2+v_luna^2-...
    (2*velocita_1*v_luna).*cos(phi_1+gamma_1));
supporto_2=asin(velocita_1./velocita_2.*cos(lambda_1+gamma_1+phi_1)-...
    v_luna./velocita_2.*cos(lambda_1));
energia_234=velocita_2.^2/2-mu_luna/raggio_soi_luna;
momento_234=raggio_soi_luna*velocita_2.*sin(supporto_2);
p_234=momento_234.^2/mu_luna;
e_234=sqrt(1+2*energia_234.*momento_234.^2/mu_luna^2);
a_234=p_234./(1-e_234.^2);
ni_2=-acos((p_234-raggio_soi_luna)./(raggio_soi_luna*e_234));
F_2=2*atanh(sqrt((e_234-1)./(e_234+1)).*tan(ni_2/2));
M_2=e_234.*sinh(F_2)-F_2;

% posizione, velocit� e angolo di volo nel punto 3
r_3=p_234./(1+e_234);
velocita_3=sqrt(2*(energia_234+mu_luna./r_3));
phi_3=0;
ni_3=0;

% calcolo il deltaV totale per cattura in un orbita lunare
deltaV1=velocita_0-sqrt(mu_terra/r_0);
a_mo=(r_3+40000)/2;
velocita_3_mo=sqrt(2*mu_luna./r_3-mu_luna./a_mo);
deltaV2=velocita_3-velocita_3_mo;
deltaVtot=deltaV1+deltaV2;

% interpolo gli indici di r_3 corrispondenti o vicini a r_periselenio scelto
ind_min=find(r_3==min(r_3));
r_3a=r_3(1:ind_min);
r_3b=r_3(ind_min:length(r_3));
ind1=interp1(r_3a,1:length(r_3a),r_periselenio,'nearest');
ind2=interp1(r_3b,1:length(r_3b),r_periselenio,'nearest');
indici_trovati=[ind1 ind_min+ind2-1];

% posizione, velocit� e angolo di volo nel punto 4
r_4=raggio_soi_luna;
velocita_4=velocita_2;
supporto_4=pi-supporto_2;
ni_4=-ni_2;
F_4=2*atanh(sqrt((e_234-1)./(e_234+1)).*tan(ni_4/2));
M_4=e_234.*sinh(F_4)-F_4;
tempo_volo_24=(M_4-M_2).*sqrt(-a_234.^3/mu_luna);

% posizione, velocit� e angolo di volo nel punto 5
lambda_5=lambda_1+2*[-ni_2(1:ind_min) ni_2(ind_min+1:end)];
r_5=sqrt(raggio_soi_luna^2+D_terra_luna^2+...
    2*D_terra_luna*raggio_soi_luna*cos(lambda_5));
velocita_5=sqrt((velocita_4.*cos(supporto_4+lambda_5)).^2+...
    (velocita_4.*sin(supporto_4+lambda_5)+v_luna).^2);
energia_5=velocita_5.^2/2-mu_terra./r_5;
a_5=-mu_terra./(2*energia_5);

chi_1=atan2((velocita_4.*sin(supporto_4+lambda_5)+v_luna),...
    (velocita_4.*cos(supporto_4+lambda_5)));
chi_2=asin(raggio_soi_luna./r_5.*sin(lambda_5));
chi_3=chi_1-chi_2;

momento_5=r_5.*velocita_5.*sin(chi_3);
e_5=sqrt(1+2*energia_5.*(momento_5/mu_terra).^2);
p_5=a_5.*(1-e_5.^2);
ni_5=acos((p_5-r_5)./(r_5.*e_5));            % non so piu se con - o +
phi_5=acos(momento_5./(r_5.*velocita_5));    % non so piu se con - o +


% traccio figure per aiutare la comprensione
figure(7)
plot(lambda_1*180/pi,momento_234);
hold on
plot(lambda_1(ind_min)*180/pi*ones(100),linspace(-6*10^4,6*10^4,100),'k--');
hold off
title('momento tratto luna(lambda1)');
xlabel('lambda1');
ylabel('momento');
text(lambda_1(ind_min)*180/pi-25,0,'orbita retrograda','FontWeight','bold');
text(lambda_1(ind_min)*180/pi+10,0,'orbita diretta','FontWeight','bold');
xlim([lambda_1(ind_min)*180/pi-30,lambda_1(ind_min)*180/pi+30]);
ax=gca;
ax.YAxis.Exponent=0;
grid on

figure(8)
plot(lambda_1*180/pi,r_3);
hold on
plot(lambda_1(ind_min)*180/pi*ones(100),linspace(0,35000,100),'k--');
plot(linspace(lambda_1(ind_min)*180/pi-30,...
    lambda_1(ind_min)*180/pi+30,100),r_periselenio*ones(100),'k--');
hold off
title('rp(lambda1)');
xlabel('lambda1');
ylabel('rp');
text(lambda_1(ind_min)*180/pi-25,5000,'orbita retrograda','FontWeight','bold');
text(lambda_1(ind_min)*180/pi+10,5000,'orbita diretta','FontWeight','bold');
xlim([lambda_1(ind_min)*180/pi-30,lambda_1(ind_min)*180/pi+30]);
ax=gca;
ax.YAxis.Exponent=0;
grid on

% trovo i valori voluti con le due scelte di r_3
lambda_1=lambda_1(indici_trovati);
r_1=r_1(indici_trovati);
velocita_1=velocita_1(indici_trovati);
phi_1=phi_1(indici_trovati);
gamma_1=gamma_1(indici_trovati);
ni_1=ni_1(indici_trovati);
tempo_volo_01=tempo_volo_01(indici_trovati);

velocita_2=velocita_2(indici_trovati);
energia_234=energia_234(indici_trovati);
momento_234=momento_234(indici_trovati);
ni_2=ni_2(indici_trovati);
tempo_volo_24=tempo_volo_24(indici_trovati);

velocita_3=velocita_3(indici_trovati);
r_3=r_3(indici_trovati);
deltaVtot=deltaVtot(indici_trovati);

velocita_4=velocita_4(indici_trovati);
ni_4=ni_4(indici_trovati);

lambda_5=lambda_5(indici_trovati);
r_5=r_5(indici_trovati);
velocita_5=velocita_5(indici_trovati);
phi_5=phi_5(indici_trovati);
ni_5=ni_5(indici_trovati);
energia_5=energia_5(indici_trovati);
momento_5=momento_5(indici_trovati);
