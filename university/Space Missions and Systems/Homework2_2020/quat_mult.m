function [delta_q]=quat_mult(qtrue,qopt)
% qtrue   quaternione vero
% qopt    quaternione trovato dall'algoritmo

[n1,m1]=size(qtrue);
if n1~=1 || m1~=4
    error('qtrue must be a 1x4 vector');
end

[n2,m2]=size(qopt);
if n2~=1 || m2~=4
    error('qopt must be a 1x4 vector');
end

qopt(1:3)=-qopt(1:3);
delta_q=[qtrue(4)  qtrue(3) -qtrue(2) qtrue(1) ;...
        -qtrue(3)  qtrue(4)  qtrue(1) qtrue(2) ;...
         qtrue(2) -qtrue(1)  qtrue(4) qtrue(3) ;...
        -qtrue(1) -qtrue(2) -qtrue(3) qtrue(4)]*qopt';
end

