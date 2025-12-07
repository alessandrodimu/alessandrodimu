function[qopt,Ptt]=quest(OBS,meass,indexspan)
% QUEST algorithm
%
% OBS         retains information over W, V and sigma for the assigned sensor
% meass       rows to be considered for further calculation in the algorithm
% indexspan   index interval over which to estimate qopt

% intialization
qopt=zeros(length(indexspan),4);
Ptt=zeros(3,3,length(indexspan));

for i=indexspan
    sigma=OBS(meass,7,i);
    sigmatot2=sum(1./(sigma.^2));
    sigmatot=sqrt(1/sigmatot2);
    
    a=(sigmatot./sigma).^2;
    a=a/sum(a);
    
    W=OBS(meass,1:3,i);
    V=OBS(meass,4:6,i);
    
    B=zeros(3,3);
    for ii=1:length(meass)
        B=B+a(ii)*(W(ii,:)'*V(ii,:));
    end
      
    Z=zeros(3,1);
    for ii=1:length(meass)
        Z=Z+a(ii)*cross(W(ii,:),V(ii,:))';
    end
    S=B'+B;
    sig=trace(B);
    
    K=[S-sig*eye(3) Z ; Z' sig];
    
%     [V,~]=eig(K);
%     qopt(i,:)=V(:,4)';

    adjS=det(S)*inv(S);
    k=trace(adjS);
    lambda0=1;
    for ii=1:3
        lambda=lambda0;
        alpha=lambda^2-sig^2+k;
        beta=lambda-sig;
        gamma=(lambda+sig)*alpha-det(S);
        X=(alpha*eye(3,3)+beta*S+S^2)*Z;
        
        % note2self: |gamma| cannot be allowed to become too small, otherwise
        % the method of sequential rotations need to be invoked.
        phi=gamma*(lambda-sig)-Z'*X;
        phi_dot=(gamma^2+norm(X)^2)/gamma;
        lambda0=lambda-phi/phi_dot;
    end
    qopt(i,:)=1/sqrt(gamma^2+norm(X)^2)*[X ; gamma];
    
    supp1=zeros(3,3);
    for ii=1:length(meass)
        supp1=supp1+a(ii)*(W(ii,:)'*W(ii,:));
    end
    Pqq=0.25*sigmatot^2*inv(eye(3)-supp1);
    Ptt(:,:,i)=4*Pqq;
end
qopt=qopt(indexspan,:);
Ptt=Ptt(:,:,indexspan);

