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

dMp = 16;%degrees of the LF polynomial in the independent variable x
degmax = 18;%degree of the polynomials in the integration by parts
dRegindep = 18;%degrees of the multipliers for the local checks

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
expr  = C*u_x^2 - u^2;%considering the nonlinear terms


[F,G,H,Z,IPMON] = SOSINTCONSTR_MAT(expr,Du);

neg = size(G.g,1);%number of elements in matrix G
neh = size(H.h,1);%number of elements in matrix H

[nrep,~] = size(F);%number of elements in matrix F

%--------------------------------------------------------------------------
%CREATES THE VARIABLES AND REPLACES THEM IN THE ENTRIES
%BUILDING THE POLYNOMIAL MATRIX REPRESENTATION
%--------------------------------------------------------------------------


[prog,gx] = sospolymatrixvar(prog,monomials(x,degp), [neg,1]);%it is a quadratic function on the dependent variables
[prog,hx] = sospolymatrixvar(prog,monomials(x,degh), [neh,1]);%it is a quadratic function on the dependent variables
hxdiff = diff(hx,x);



Fv = F;
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

 
FSOL = F;
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