%Stability Conditions for the Burgers equation with reaction term
% u_t =  u_xx -u*u_x + beta*u^2
%
%Considers functionals of the form
% /
% | p(x)u^2 dx
% /D
%
tic

clear all;

n = 1;%number of dependent variables
ns = 1;%number of the spatial variables (ni + 1)
ord = 2;%this is the order of the derivatives appearing on your monomials

dMp = 6;%degrees of the LF polynomial in the independent variable x
degmax = 8;%degree of the polynomials in the integration by parts
dRegindep = 8;%degrees of the multipliers for the local checks

Du = sym('Du',[n,ord+1]); Du = sym(Du,'real');
DULIM = sym('DULIM',[n,ord]); DULIM = sym(DULIM,'real');
DLLIM = sym('DLLIM',[n,ord]); DLLIM = sym(DLLIM,'real');
x = sym('x',[ns,1]);x = sym(x,'real');

u = Du(1,1);u_x = Du(1,2);u_xx = Du(1,3);
%this is the time-derivative
beta = .4*(pi/sqrt(2));
u_t = u_xx -u*u_x + beta*u^2;

%B MATRIX SETS THE BOUNDARY CONDITIONS
%this is Dirichlet
B = [1 0 0 0; 0 0 1 0];

%each multiplies ub = [u(1)  u_x(1)  u(0)  u_x(0) ] to get B*vec_UL = 0
vec_UL = [DULIM(1:n*ord)  DLLIM(1:n*ord)].';

%SETS THE SOS program
dvartable = [Du(1:n*(ord+1)) vec_UL']; ivartable = [x];
prog = sosprogram([dvartable ivartable]);

%the weighting polynomial p(x)
degp = 0:dMp;%the function must be positive at least in the interval
degh = 0:degmax; 
degreg = 0:dRegindep;
 

%--------------------------------------------------------------------------
%DEFINES THE LYAPUNOV FUNCTION
%--------------------------------------------------------------------------
nV = n;
nP = nV*(nV+1)/2;
p = sym('p',[ns,1]);p = sym(p,'real');
P = p;
%P = mpvar('P',[nV,nV]);%the lines below define a symmetric matrix
%col = vec(repmat(1:nV,[nV 1]));row = vec(repmat(1:nV,[1 nV])');
%P(row>=col) = p; P = (P+P'-diag(P(row==col)));%just builds the symmetric matrix

%this is the derivative of the Lyapunov function
Vdot = 2*u*P*u_t;%considering the nonlinear terms

[F,G,H,Z,IPMON] = SOSINTCONSTR_MAT(-Vdot,Du);

neg = size(G.g,1);%number of elements in matrix G
neh = size(H.h,1);%number of elements in matrix H

[nrep,~] = size(F);%number of elements in matrix F

%--------------------------------------------------------------------------
%CREATES THE VARIABLES AND REPLACES THEM IN THE ENTRIES
%BUILDING THE POLYNOMIAL MATRIX REPRESENTATION
%--------------------------------------------------------------------------

%the weighting polynomial p(x)
[prog,px] = sospolymatrixvar(prog,monomials(x,degp), [nP,1]);%it is a quadratic function on the dependent variables
[prog,gx] = sospolymatrixvar(prog,monomials(x,degp), [neg,1]);%it is a quadratic function on the dependent variables
[prog,hx] = sospolymatrixvar(prog,monomials(x,degh), [neh,1]);%it is a quadratic function on the dependent variables
hxdiff = diff(hx,x);



Fv = subs(F,p,px);
if ~isempty(G.g)
    Gv = subs(G.G,G.g,gx);
else
    Gv = zeros(size(H.H));
end
Hv = subs(H.H,H.h,hx);
Hdiffv = subs(H.Hdiff,H.h,hxdiff);

Mx = Fv+Gv+(Hv+Hdiffv);


%--------------------------------------------------------------------------
%SETS THE CONDITIONS IN THE DOMAIN
%--------------------------------------------------------------------------
Llim = 0; Ulim = 1;
pD = -(x(1)-Llim)*(x(1)-Ulim);%pD is positive in the set of interest


eps = 0.001;
Px = subs(P,p,px);
%prog = sosmatrixineq(prog,Px-eps*eye(nV),'quadraticMineq');


[prog,NPx] = sospolymatrixvar(prog,monomials(x,degreg),[nV,nV],'symmetric');
prog = sosmatrixineq(prog,NPx,'quadraticMineq');
LocPx = (Px-eps*eye(nV))  - pD*NPx;
prog = sosmatrixineq(prog,LocPx,'quadraticMineq');


[prog,Nx] = sospolymatrixvar(prog,monomials(x,degreg),[nrep,nrep],'symmetric');
prog = sosmatrixineq(prog,Nx,'quadraticMineq');
Locx = Mx  - pD*Nx;
prog = sosmatrixineq(prog,Locx,'quadraticMineq');

 
%--------------------------------------------------------------------------
%AND THE VALUES FOR H AT THE BOUNDARIES
%--------------------------------------------------------------------------
bcvar = B;%for this example we take it as a matrix
[prog,hb] =  BCCONSTRAINTS(prog,x,[DULIM; DLLIM],hx,IPMON.Lf,IPMON.Lof,bcvar,'linearBC');
 

toc
prog = sossolve(prog);
toc

%--------------------------------------------------------------------------
%EXTRACTS THE SOLUTION
%--------------------------------------------------------------------------

SOL_coeffsV = sosgetsol(prog,px);
Px = subs(P,p,SOL_coeffsV);
FSOL = subs(F,p,SOL_coeffsV);
if ~isempty(G.g)
    SOL_coeffsG = sosgetsol(prog,gx);
    GSOL = subs(G.G,G.g,SOL_coeffsG);
else
    GSOL = zeros(size(H.H));
end
SOL_coeffsH = sosgetsol(prog,hx);
HSOL = subs(H.H,H.h,SOL_coeffsH);
SOL_coeffsHdiff  = diff(SOL_coeffsH,x);
HdiffSOL = subs(H.Hdiff,H.h,SOL_coeffsHdiff);
HbarSOL = HSOL+HdiffSOL;
 
Mx = FSOL+GSOL+HbarSOL;


%--------------------------------------------------------------------------
%PLOTS THE RESULT
%--------------------------------------------------------------------------

Llim = 0; Ulim = 1;
Nfact = subs(Px,x,0);%the value at zero
len = Llim:0.01:Ulim; Gnum = double(Nfact)^(-1)*subs(Px,x,len);
figure;plot(len,Gnum,'LineWidth',1.5)
xlabel('x')
ylabel('p(x)')
 