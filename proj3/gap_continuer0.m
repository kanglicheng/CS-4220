%
% Fixed point continuation (controlling voltage)
%

% === Load gap finite element code
% argument in gap_setup:  'c': for cantilever (this problem)
%                         80: number of elements in the beam
pgap = gap_setup('c', 80);

nelt = pgap.nelt;
Ce = pgap.Ce;
M = pgap.M;
K = pgap.K;
N = pgap.N;
wg = pgap.wg;
Itip = pgap.Itip;
ndof = pgap.ndof; %ndof is the "active" degrees of freedom

% === Continuation loop
% V: voltage varying between 0 and 10
%
numV=200;
V=linspace(0,6.975,numV);

utip=zeros(numV,1);
u=zeros(ndof,1); % initial guess for u
for kk=1:numV
    [u, rnorms]=gap_solve_fpV(pgap, u, V(kk));
    utip(kk)=u(pgap.Itip);
end
% === Plot V vs tip displacement
plot(V, utip);
xlabel('Voltage');
ylabel('Tip dispalcement');