function dw = dyn_eq(t,w)
% DYN_EQ contains the alternate forms of the dynamics equations for a
% spacecraft with three wheels and one damper (homework 1)

eta=w(1:10);
xi=w(11);
E=w(12);

v_P=w(1:3);
omega=w(4:6);
xi_dot=w(7);
omega_S=w(8:10);

dw=zeros(12,1);

global m_b m_w m_d cd kd kw b_w I_s I_t in

b1=[1;0;0];
b2=[0;1;0];
b3=[0;0;1];
n=b1;
r_d=xi*n;

vp_til=skew(v_P);
omega_til=skew(omega);
n_til=skew(n);
r_d_til=skew(r_d);


% estimation of matrix MT
M=m_b+3*m_w+m_d;

S_BP=0;
S_WP=m_w*(b_w*b1)+m_w*(b_w*b2)+m_w*(b_w*b3);
S_DP=m_d*xi*b1;
S_P=S_BP+S_WP+S_DP;

J_BP=diag([350,300,400]);
J_WP=diag([I_s,I_t,I_t])+m_w*(norm(b_w*b1)^2*eye(3,3)-(b_w*b1)*(b_w*b1)')+...
     diag([I_t,I_s,I_t])+m_w*(norm(b_w*b2)^2*eye(3,3)-(b_w*b2)*(b_w*b2)')+...
     diag([I_t,I_t,I_s])+m_w*(norm(b_w*b3)^2*eye(3,3)-(b_w*b3)*(b_w*b3)');
J_DP=m_d*(norm(xi*b1)^2*eye(3,3)-(xi*b1)*(xi*b1'));
J_P=J_BP+J_WP+J_DP;

S_Ptil=skew(S_P);

MT=[M*eye(3,3) -S_Ptil m_d*b1 zeros(3,1) zeros(3,1) zeros(3,1);...
    S_Ptil J_P zeros(3,1) I_s*b1 I_s*b2 I_s*b3;...
    m_d*b1' zeros(1,3) m_d 0 0 0;...
    zeros(1,3) I_s*b1' 0 I_s 0 0;...
    zeros(1,3) I_s*b2' 0 0 I_s 0;...
    zeros(1,3) I_s*b3' 0 0 0 I_s];


% estimation of matrix MT_dot
S_P_dot=m_d*xi_dot*b1;
J_P_dot=m_d*[0 0 0;...
             0 2*xi*xi_dot 0;...
             0 0 2*xi*xi_dot];
S_P_dottil=skew(S_P_dot);

MT_dot=[zeros(3,3) -S_P_dottil zeros(3,4);...
        S_P_dottil J_P_dot zeros(3,4);...
        zeros(4,10)];


% estimation of matrix A
A=[-omega_til zeros(3,3) zeros(3,4);...
   -vp_til -omega_til zeros(3,4);...
   zeros(1,10);...
   zeros(3,10)];
    

% alternate forms of the dynamics equations
switch in
    case 1      % c does not depend on omega_dot
        ga1=0;
        ga2=0;
        ga3=0;
        
        c=[zeros(3,1);...
           zeros(3,1);...
           m_d*omega'*n_til*(v_P-r_d_til*omega)-cd*xi_dot-kd*xi;...
           ga1;...
           ga2;...
           ga3];
        
        eta_dot=MT\((A*MT-MT_dot)*eta+c);
       
        
    case 2      % c depends on omega_dot in a certain way
        B=[zeros(3,10);...
           zeros(3,10);...
           zeros(1,10);...
           zeros(1,3) I_s*b1' zeros(1,4);...
           zeros(1,3) I_s*b2' zeros(1,4);...
           zeros(1,3) I_s*b3' zeros(1,4)];
       
        c0=[zeros(3,1);...
            zeros(3,1);...
            m_d*omega'*n_til*(v_P-r_d_til*omega)-cd*xi_dot-kd*xi;...
            zeros(3,1)];
        
        eta_dot=(MT-B)\((A*MT-MT_dot)*eta+c0);
        
        omega_dot=eta_dot(4:6);
        ga1=I_s*b1'*omega_dot;
        ga2=I_s*b2'*omega_dot;
        ga3=I_s*b3'*omega_dot;

        
    case 3      % c depends on omega_dot in a certain way
        ga2=-kw*omega_S(2);
        ga3=-kw*omega_S(3);
        
        B=[zeros(3,10);...
           zeros(3,10);...
           zeros(1,10);...
           zeros(1,3) I_s*b1' zeros(1,4);...
           zeros(1,10);...
           zeros(1,10)];
       
        c0=[zeros(3,1);...
            zeros(3,1);...
            m_d*omega'*n_til*(v_P-r_d_til*omega)-cd*xi_dot-kd*xi;...
            0;...
            ga2;
            ga3];
        
        eta_dot=(MT-B)\((A*MT-MT_dot)*eta+c0);
        
        omega_dot=eta_dot(4:6);
        ga1=I_s*b1'*omega_dot;
        

    case 4      % c does not depend on omega_dot
        ga1=-kw*omega_S(1);
        ga2=-kw*omega_S(2);
        ga3=-kw*omega_S(3);
        
        c=[zeros(3,1);...
           zeros(3,1);...
           m_d*omega'*n_til*(v_P-r_d_til*omega)-cd*xi_dot-kd*xi;...
           ga1;...
           ga2;...
           ga3];
        
        eta_dot=MT\((A*MT-MT_dot)*eta+c);
end


% mechanical energy rate
E_dot=-cd*xi_dot^2+ga1*omega_S(1)+ga2*omega_S(2)+ga3*omega_S(3);


% integration of the parameters
dw(1:10)=eta_dot;
dw(11)=xi_dot;
dw(12)=E_dot;

end

