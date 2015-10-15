%
% Compute pull-in voltage vs beam length 
% for lengths from 20 to 100 microns
%
len1=20;
len2=100;
L_beam=linspace(len1, len2);
num=length(L_beam);
VV=zeros(num,1);
for iter=num:-1:1
    %
    % change the length of the beam:
    %
    pgap = gap_setup('c', 80, L_beam(iter));
    [u,V]=gap_pullin(pgap);
    VV(iter)=V;
end

plot(L_beam, VV);
xlabel('Beam Length [micron]');
ylabel('Pull-in voltage [V]');