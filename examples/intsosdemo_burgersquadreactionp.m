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

dMp = 5;%degrees of the LF polynomial in the independent variable x
degmax = 7;%degree of the polynomials in the integration by parts
dRegindep = 7;%degrees of the multipliers for the local checks

Du = mpvar('Du',[n,ord+1]);
DULIM = mpvar('DULIM',[n,ord]);
DLLIM = mpvar('DLLIM',[n,ord]);
x = mpvar('x',[ns,1]);

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
Zp = monomials(Du(:,1),2);nZp = size(Zp,1);
p = mpvar('p',[1 nZp]);
Vp = p*Zp;


%--------------------------------------------------------------------------
%OBTAINS THE DERIVATIVE AND THE FTC TERMS
%--------------------------------------------------------------------------

%this is the derivative of the Lyapunov function
exprf = -(jacobian(Vp,Du(:,1))*u_t);%considering the nonlinear terms

[exprh,Zd,IPMON] = SOSINTCONSTR(exprf,Du);

neh = size(IPMON.Lf,1);%number of elements in matrix H
 
%--------------------------------------------------------------------------
%CREATES THE VARIABLES AND REPLACES THEM IN THE ENTRIES
%BUILDING THE POLYNOMIAL MATRIX REPRESENTATION
%--------------------------------------------------------------------------

%the weighting polynomial p(x)
[prog,px] = sospolymatrixvar(prog,monomials(x,degp), [nZp,1]);%it is a quadratic function on the dependent variables
[prog,hx] = sospolymatrixvar(prog,monomials(x,degh), [neh,1]);%it is a quadratic function on the dependent variables
hxdiff = diff(hx,x);

exprfv = subs(exprf,p,px);
exprhv = subs(exprh.h,exprh.hvar,hx);
exprhdiffv = subs(exprh.hdiff,exprh.hvar,hxdiff);

exprx = exprfv+(exprhv+exprhdiffv);


%--------------------------------------------------------------------------
%SETS THE CONDITIONS IN THE DOMAIN
%--------------------------------------------------------------------------
Llim = 0; Ulim = 1;
pD = -(x(1)-Llim)*(x(1)-Ulim);%pD is positive in the set of interest

eps = 0.001;
epsx = eps*Du(:,1).'*Du(:,1);
Vpx = subs(Vp,p,px);
[prog,NPx] = sospolymatrixvar(prog,monomials(x,degreg),[1,nZp]);
NV = NPx*Zp;
prog = sosineq(prog,NV,'sparsemultipartite',{x,vec(Du)});
LocVx = (Vpx-epsx)  - pD*NV;
prog = sosineq(prog,LocVx,'sparsemultipartite',{x,vec(Du)});

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

SOL_coeffsV = sosgetsol(prog,px);
Px = subs(Vp,p,SOL_coeffsV);
FSOL = subs(exprf,p,SOL_coeffsV);
SOL_coeffsH = sosgetsol(prog,hx);
HSOL = subs(exprh.h,exprh.hvar,SOL_coeffsH);
SOL_coeffsHdiff  = diff(SOL_coeffsH,x);
HdiffSOL = subs(exprh.hdiff,exprh.hvar,SOL_coeffsHdiff);
HbarSOL = HSOL+HdiffSOL;
 
Mx = FSOL+HbarSOL;


%--------------------------------------------------------------------------
%PLOTS THE RESULT
%--------------------------------------------------------------------------

Llim = 0; Ulim = 1;
Px = subs(Px,Du(1),1);%the value at zero
Nfact = subs(Px,x,0);
len = Llim:0.01:Ulim; Gnum = double(Nfact)^(-1)*subs(Px,x,len);
figure;plot(len,Gnum,'LineWidth',1.5)
xlabel('x')
ylabel('p(x)')
 