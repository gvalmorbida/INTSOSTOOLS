%Stability Conditions for a Reaction--Diffusion system
% u_t =   u_xx  + u
%
%Considers functionals of the form
% /
% | p(x)u^2 dx
% /D
%
tic
clear all;

n = 4;%number of dependent variables
ns = 1;%number of the spatial variables (ni + 1)
ord = 2;%this is the order of the derivatives appearing on your monomials

dMp = 2;%degrees of the LF polynomial in the independent variable x
degmax = 4;%degree of the polynomials in the integration by parts
dRegindep = 4;%degrees of the multipliers for the local checks

Du = sym('Du',[n,ord+1]); Du = sym(Du,'real');
DULIM = sym('DULIM',[n,ord]); DULIM = sym(DULIM,'real');
DLLIM = sym('DLLIM',[n,ord]); DLLIM = sym(DLLIM,'real');
x = sym('x',[ns,1]);x = sym(x,'real');

%sets the parameters
Re =0.1;
Ro = 0;
alpha = pi/2;

%these are the terms for the time-derivative
u_dif = (1/(4*Re))*[Du(:,3)];%the diffusion terms

u_t = u_dif;%the linear terms on the time-derivative

%B MATRIX SETS THE BOUNDARY CONDITIONS
bc = [DLLIM(:,1); DULIM(:,1)];
%each multiplies ub = [u(1) v(1) u(0) v(0) u_x(1) v_x(1) u_x(0) v_x(0)] to get B*ub = 0
B = double(jacobian(bc,[ vec(DULIM); vec(DLLIM)]));%B MATRIX SETS THE BOUNDARY CONDITIONS


%each multiplies ub = [u(1)  u_x(1)  u(0)  u_x(0) ] to get B*vec_UL = 0
vec_UL = [DULIM(1:n*ord)  DLLIM(1:n*ord)].';

%SETS THE SOS program
dvartable = [Du(1:n*(ord+1)) vec_UL']; ivartable = [x];
prog = sosprogram([dvartable ivartable]);


%Define the degree of the Lyapunov funciton and the multipliers on the
%spatial variable
degp = 0:dMp;%the function must be positive at least in the interval
degh = 0:degmax; 
degreg = 0:dRegindep;
 

%--------------------------------------------------------------------------
%DEFINES THE LYAPUNOV FUNCTION
%--------------------------------------------------------------------------
%a quadratic for
degp = 2;
Zp = monomials(Du(:,1),degp);nZp = size(Zp,1);
p = sym('p',[1 nZp]);p = sym(p,'real');
Vp = p*Zp;


%--------------------------------------------------------------------------
%OBTAINS THE DERIVATIVE AND THE FTC TERMS
%--------------------------------------------------------------------------

%this is the derivative of the Lyapunov function
exprf = -(jacobian(Vp,Du(:,1))*u_t);%considering the nonlinear terms

[exprh,Zd,IPMON] = SOSINTCONSTR(exprf,Du);
 
%--------------------------------------------------------------------------
%CREATES THE VARIABLES AND REPLACES THEM IN THE ENTRIES,
%BUILDING THE POLYNOMIAL MATRIX REPRESENTATION
%--------------------------------------------------------------------------

%the weighting polynomial p(x)
degp = 0:dMp;%the function must be positive at least in the interval
[prog,px] = sospolymatrixvar(prog,monomials(x,degp), [nZp,1]);%it is a quadratic function on the dependent variables
neh = size(IPMON.Lf,1);%number of elements in matrix H
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
 

%--------------------------------------------------------------------------
%SOLVE THE SOSP
%--------------------------------------------------------------------------
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