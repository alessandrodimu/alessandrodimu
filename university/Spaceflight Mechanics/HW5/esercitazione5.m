% ESERCITAZIONE 5
% PARTE 1

m_u3=2000;
dv=9300;
g0=9.81;
Isp1=268;
Isp2=285;
Isp3=320;
eps1=9;
eps2=9.7;
eps3=7;


k=exp(dv/(g0*(Isp1+Isp2+Isp3)));
lambda3=(eps3/k-1)/(eps3-1);
m_i3=m_u3/lambda3;
m_f3=m_i3/k;
m_p3=m_i3-m_f3;
m_s3=m_i3-m_u3-m_p3;
m_tot3=m_p3+m_s3;
tot3_s3=m_tot3/m_s3;
DV_3=g0*Isp3*log(k);

m_u2=m_i3;
lambda2=(eps2/k-1)/(eps2-1);
m_i2=m_u2/lambda2;
m_f2=m_i2/k;
m_p2=m_i2-m_f2;
m_s2=m_i2-m_u2-m_p2;
m_tot2=m_p2+m_s2;
tot2_s2=m_tot2/m_s2;
DV_2=g0*Isp2*log(k);

m_u1=m_i2;
lambda1=(eps1/k-1)/(eps1-1);
m_i1=m_u1/lambda1;
m_f1=m_i1/k;
m_p1=m_i1-m_f1;
m_s1=m_i1-m_u1-m_p1;
m_tot1=m_p1+m_s1;
tot1_s1=m_tot1/m_s1;
DV_1=g0*Isp1*log(k);