function dw=state_eqs(t,w,u)

global mu L2

x=w(1);
y=w(2);
z=w(3);
xdot=w(4);
ydot=w(5);
zdot=w(6);

dw=zeros(6,1);

% state equations (non linear system)
dw(1)=xdot;
dw(2)=ydot;
dw(3)=zdot;
dw(4)=2*ydot+(L2+x)-(1-mu)*((L2+x)+mu)/(((mu+(L2+x))^2+y^2+z^2)^(3/2))...
      -mu*((L2+x)+mu-1)/(((L2+x)+mu-1)^2+y^2+z^2)^(3/2)+u(1);
dw(5)=-2*xdot+y-(1-mu)*y/(((mu+(L2+x))^2+y^2+z^2)^(3/2))...
      -mu*y/(((L2+x)+mu-1)^2+y^2+z^2)^(3/2)+u(2);
dw(6)=-(1-mu)*z/(((mu+(L2+x))^2+y^2+z^2)^(3/2))...
      -mu*z/(((L2+x)+mu-1)^2+y^2+z^2)^(3/2)+u(3);

end

