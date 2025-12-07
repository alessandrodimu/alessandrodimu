function[qdet,Ptt]=triad(OBS,meass,indexspan)
% TRIAD algortihm
%
% OBS           retains information over W, V and sigma for the assigned sensor
% meass         rows to be considered for further calculation in the algorithm
%               only 2 rows will be considered here
% indexspan     index interval over which to estimate qopt

% check for errors in the model
if length(meass)~=2
    error('Only 2 measurements must be taken into account');
end
    
% initialization
qdet=zeros(length(indexspan),4);
Ptt=zeros(3,3,length(indexspan));

for i=indexspan
    % get values from OBS matrix
    V1=OBS(meass(1),4:6,i)';
    V2=OBS(meass(2),4:6,i)';
    
    W1=OBS(meass(1),1:3,i)';
    W2=OBS(meass(2),1:3,i)';
    
    % create two triads of orthonormal reference
    r1=V1;
    r2=cross(V1,V2)/norm(cross(V1,V2));
    r3=cross(V1,cross(V1,V2))/norm(cross(V1,V2));
    
    s1=W1;
    s2=cross(W1,W2)/norm(cross(W1,W2));
    s3=cross(W1,cross(W1,W2))/norm(cross(W1,W2));
    
    Mref=[r1 r2 r3];
    Mobs=[s1 s2 s3];
    
    A=Mobs*Mref';
    
    q4=0.5*sqrt(1+A(1,1)+A(2,2)+A(3,3));
    qdet(i,:)=[0.25/q4*(A(2,3)-A(3,2)) ...
               0.25/q4*(A(3,1)-A(1,3)) ...
               0.25/q4*(A(1,2)-A(2,1)) ...
               q4];
    
    sig1=OBS(meass(1),7,i);
    sig2=OBS(meass(2),7,i);
    Ptt(:,:,i)=sig1^2*eye(3,3)+((sig2^2-sig1^2)*(W1*W1')+...
                sig1^2*dot(W1,W2)*(W1*W2'+W2*W1'))/norm(cross(W1,W2))^2;
end
qdet=qdet(indexspan,:);
Ptt=Ptt(:,:,indexspan);

