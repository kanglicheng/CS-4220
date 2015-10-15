%
% Newton point continuation (controlling tip displacement)
%

pgap = gap_setup('c', 80, 70);
num=100;
nelt = pgap.nelt;
Ce = pgap.Ce;
M = pgap.M;
K = pgap.K;
N = pgap.N;
wg = pgap.wg;
Itip = pgap.Itip;
ndof = pgap.ndof;

grange=0.8*pgap.g;
tipDisp=linspace(0.3, grange, num);
tipDisp=-tipDisp;
V=0;
u=zeros(ndof,1);
Volt=zeros(num,1);
PDstate=zeros(num,1);
% === Continuation loop
for iter=1:num
    d=tipDisp(iter);
    [u, V, rnorms] = gap_solve_Nd(pgap, u, V, d);
    fprintf('\n***d=%g,   V=%g\n', d, V);
    Volt(iter)=V;
    %%%%%% adding a section of test code for PD of J:
    [f, Ju, JV]=gap_force(pgap,u, V);
    JF=K-Ju;
    [R,p]=chol(JF);
    PDstate(iter)=p;    
end



% === Plot V vs tip displacement
plot(Volt, -tipDisp);
grid on;
ylabel('-Tip displacement [microns]');
xlabel('Voltage [V]');

figure(2);
plot(Volt,PDstate);
xlabel('Voltage [V]');
