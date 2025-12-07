function [thnorm] = anorm (theta,icampo)
       %RESTITUISCE L'ANGOLO THETA NORMALIZZATO TRA:
       %  0  E PI2  SE ICAMPO = 1
       % -PI E PI   SE ICAMPO = 0
    pi2=2*pi;
	thnorm=theta-floor(theta/pi2)*pi2;
    if(icampo == 0)
        if(thnorm > pi)
            thnorm=thnorm-pi2;
        end;
    end;
