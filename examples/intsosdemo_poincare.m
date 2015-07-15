%Conditions for the Poincare inequality
%
%Solves the inequality of the form
% /
% | Cu_x^2 - u^2 dx
% /D
%
tic

%clear all;

n = 1;%number of dependent variables
ns = 1;%number of the spatial variables (ni + 1)
ord = 1;%this is the order of the derivatives appearing on your monomials

dMp = 10;%degrees of the LF polynomial in the independent variable x
degmax = 12;%degree of the polynomials in the integration by parts
dRegindep = 12;%degrees of the multipliers for the local checks

Du = sym('Du',[n,ord+1]); Du = sym(Du,'real');
DULIM = sym('DULIM',[n,ord]); DULIM = sym(DULIM,'real');
DLLIM = sym('DLLIM',[n,ord]); DLLIM = sym(DLLIM,'real');
x = sym('x',[ns,1]);x = sym(x,'real');

u = Du(1,1);u_x = Du(1,2);
%this is the parameter to optimize
C = 1.5*pi^(-2);

%B MATRIX SETS THE BOUNDARY CONDITIONS
%this is Dirichlet
B = [1 0; 0 1];

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
%OBTAINS THE DERIVATIVE AND THE FTC TERMS
%--------------------------------------------------------------------------

%this is the derivative of the Lyapunov function
exprf = C*u_x^2 - u^2;%considering the nonlinear terms

[exprh,Zd,IPMON] = SOSINTCONSTR(exprf,Du);

neh = size(IPMON.Lf,1);%number of elements in matrix H
 
%--------------------------------------------------------------------------
%CREATES THE VARIABLES AND REPLACES THEM IN THE ENTRIES
%BUILDING THE POLYNOMIAL MATRIX REPRESENTATION
%--------------------------------------------------------------------------

%the weighting polynomial p(x)
[prog,hx] = sospolymatrixvar(prog,monomials(x,degh), [neh,1]);%it is a quadratic function on the dependent variables
hxdiff = diff(hx,x);

exprfv = exprf;
exprhv = subs(exprh.h,exprh.hvar,hx);
exprhdiffv = subs(exprh.hdiff,exprh.hvar,hxdiff);

exprx = exprfv+(exprhv+exprhdiffv);


%--------------------------------------------------------------------------
%SETS THE CONDITIONS IN THE DOMAIN
%--------------------------------------------------------------------------
Llim = 0; Ulim = 1;
pD = -(x(1)-Llim)*(x(1)-Ulim);%pD is positive in the set of interest
 

nZd = size(Zd);
[prog,Nx] = sospolymatrixvar(prog,monomials(x,degreg),[1,nZd]);
NVdot = Nx*Zd;
prog = sosineq(prog,NVdot,'sparsemultipartite',{x,vec(Du)});
Locx = exprx-pD*(NVdot);
prog = sosineq(prog,Locx,'sparsemultipartite',{x,vec(Du)});

 
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



FSOL = subs(exprf);
SOL_coeffsH = sosgetsol(prog,hx);
HSOL = subs(exprh.h,exprh.hvar,SOL_coeffsH);
SOL_coeffsHdiff  = diff(SOL_coeffsH,x);
HdiffSOL = subs(exprh.hdiff,exprh.hvar,SOL_coeffsHdiff);
HbarSOL = HSOL+HdiffSOL;
 
Mx = FSOL+HbarSOL;

