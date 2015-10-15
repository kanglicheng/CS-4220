%
% Newton continuation (controlling voltage)
%

% === Load gap finite element code

pgap = gap_setup('c', 80, 68.2);

nelt = pgap.nelt;
Ce = pgap.Ce;
M = pgap.M;
K = pgap.K;
N = pgap.N;
wg = pgap.wg;
Itip = pgap.Itip;
ndof = pgap.ndof;

% === Continuation loop

numV=200;
V=linspace(6,10,numV);

utip=zeros(numV,1);
u=zeros(ndof,1); % initial guess for u
for iter=1:numV
    [u, rnorms]=gap_solve_NV(pgap, u, V(iter));
    utip(iter)=u(pgap.Itip);
end
% === Plot V vs tip displacement
plot(V, utip);
xlabel('Voltage');
ylabel('Tip dispalcement');