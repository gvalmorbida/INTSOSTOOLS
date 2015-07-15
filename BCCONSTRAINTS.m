function  [prog,hb] = BCCONSTRAINTS(prog,x,BCVAR,h,L,Lo,bg,option)
%--------------------------------------------------------------------------
%BCCONSTRAINTS.m
%
%This function takes a vector of functions (r), the corresponding list of
%monomials in the dependent variables (L,idx,Lf), the matrix describing a
%set of (linear) boundary conditions B, the boundaries, the order o and
%number of dependent variables n.
%
%THIS IS THE CODE FOR ONE SPATIAL VARIABLE
%
%inputs -           prog:           the SOS program
%                   x:              the vector with the independent variables
%                   BCVAR:          Variables associated to the boundary
%                                   conditions. The format is [BCVAR(1); BCVAR(0)]
%                                   boundary constraints
%                   h:              either a vector of polynomials
%                                   associated to the monomials of the polynomials
%                                   or a matrix with the structure for the linear case
%                   L,Lo:           the set of monomials related to the entries in h
%                   bg:             either a matrix in the linear BC case or
%                                   an expression establishing the boundary conditions for
%                                   the general case. For the matrix case,
%                                   B must multiply the vector [Dalpha_u(1);Dalpha_u(0)], Dalpha = [u; u_x; u_xx; ...]
%                   option:         option to simplify the boundary condition in the case of a quadratic inequality and
%                                   linear BC
%
%outputs -          prog:         the SOS program with the inequality constraint involving the boundary
%                   hb:           the set of constraints on the boundaries
%
%--------------------------------------------------------------------------
%
% This file is part of INTSOSTOOLS - INTSOSTOOLS ver 1.00 a plug-in to SOSOTOOLS for integral inequalities.
%
% Copyright (C)2015  G. Valmorbida, A. Papachristodoulou

% Department of Engineering Science, University of Oxford, Oxford, U.K.
%
% Send bug reports and feedback to: giorgio.valmorbida@eng.ox.ac.uk
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------

% 06/07/15 - GV


%PREVIOUS VERSION - ALLOWING FOR THE BOUNDARY TO BE VARIABLE.
%lb = bound(1);
%ub = bound(2);

lb = 0;
ub = 1;

lb = -1;
ub = 1;


%1)manipulates the symbolic variables given by the values at the
%boundary
%the order to relate with B will be e.g. [u(1); u_x(1); ... ; u(0); u_x(0); ... ]

%first extracts the order from the vector of variables at the
%boundaries
 

if nargin==8
    if option=='linearBC'
        %in this case the boundary conditions are determined by the matrix bg, 
        
        [n1,n2] = size(h);
        [n,ord] = size(BCVAR);n = n/2;
        
        if n1==n2&&n1==n*ord%THE EXPRESSION IS QUADRATIC ON THE BOUNDARY VARIABLES
            
            %1)create the block-diagonal matrix substituting the variables at the
            %boundaries
            
            Hbar = blkdiag(subs(h,x,ub), -subs(h,x,lb));%in this case h must be a matrix
            
            %2)computes the nullspace of B
            Bperp = null(bg);
            
            %3)Sets the matrix at the boundaries
            hb = Bperp'*Hbar*Bperp;
            
            %4)and sets the inequality
            if isa(hb,'polynomial')
                if hb.coefficient~=0
                    prog = sosmatrixineq(prog,-hb,'quadraticMineq');
                end
            else
                if hb~=0
                    prog = sosmatrixineq(prog,-hb,'quadraticMineq');
                end
            end
 
        else
            
            %THE BOUNDARY CONDITION IS DEFINED AS bg*BCVARvec=0
            [n,ord] = size(BCVAR);n = n/2;%the number of variables and the order of the expression
            BCVARvec = [vec((BCVAR(1:n,:))); vec((BCVAR(n+1:2*n,:)))];
            
            
            %THE MATRIX MUST HAVE FULL ROW RANK!!
            [m,~] = size(bg);
            %substitutes the variables defined by the boundary conditions
            [Q,R,e] = qr(bg,0);
            if isempty(BCVARvec(e(m+1:end)))
                BCVARvec_subs  = zeros(size(BCVARvec));
            else
                BCVARvec_subs  = subs(BCVARvec,BCVARvec(e(1:m)),-inv(R(:,1:m))*R(:,m+1:end)*BCVARvec(e(m+1:end)));
            end
            %2)build the vectors that receive the monomials evaluated at the boundaries
            %the input to MONLISTDER must be a matrix 
            ud_u = MONLISTDER(reshape(BCVARvec_subs(1:n*ord,:),n,ord),L,Lo);
            ud_l = MONLISTDER(reshape(BCVARvec_subs((n*ord)+1:2*(n*ord),:),n,ord),L,Lo);
            
            %3)build the polynomial h evaluated at the boundaries
            boundary_values = [subs(h,x,ub).' -subs(h,x,lb).']*[ud_u;ud_l];%this is the nonlinear expression in the boundary values, which will be used to define a constraint on 
            
            if isa(boundary_values,'polynomial')
                if boundary_values.coefficient~=0
                    prog = sosineq(prog,-boundary_values);
                end
            else
                if boundary_values~=0
                    prog = sosineq(prog,-boundary_values);
                end
            end
            hb = -boundary_values;
   
        end
        
    end
else
    
    [n,ord] = size(BCVAR);n = n/2;%the number of variables and the order of the expression
    BCVARvec = [vec((BCVAR(1:n,:)).'); vec((BCVAR(n+1:2*n,:)).')];
    
    %2)build the vectors that receive the monomials evaluated at the boundaries
    ud_u = MONLISTDER(BCVAR(1:n,:),L,Lo);
    ud_l = MONLISTDER(BCVAR(n+1:2*n,:),L,Lo);
    
    %3)build the polynomial h evaluated at the boundaries
    boundary_values = subs(h,x,ub).'*ud_u-subs(h,x,lb).'*ud_l;
    

    
    %4)complete the expression by adding terms which multiply the BC
    %expressions
    
    %4a)obtain the degree of the boundary expressions vector on the
    %variables BCVAR and the degree of the integrations by part expressions
    
    %gets the degree of the expression in h
    degzm = min(sum(L,2));
    degzM = max(sum(L,2));
    
    %gets the degree of the boundary conditions expressions, in most of the
    %cases it will be a linear expression but this code allows for very
    %general expressions
    
    if isfield(prog,'symvartable')
        
        %gets the degree of the integration by parts expressions
        etaK = sym(zeros(size(bg,1),1));%create the symbolic variable to receive the multipliers
        charvartable = converttochar([vec(BCVAR)']);
        for i = 1:size(bg,1)
            degcheck = evalin(symengine,char(bg(i)));
            
            degcheck = feval(symengine,'expand',degcheck);
            degcheck = feval(symengine,'collect',degcheck,charvartable);
            degcheckmon = feval(symengine,'poly2list',degcheck,charvartable);
            mondeg = zeros(length(degcheckmon),1);
            for j = 1:length(degcheckmon)
                test = degcheckmon(j);
                mondeg(j) = sum(double(test(2)));
            end
            degBCM = degzM- max(mondeg);%writes the maximum degre to be added to each BC constraint
            degBCm = degzm- min(mondeg);%writes the minimum degre to be added to each BC constraint
            
            %4b)compute the degree of the multiplier and builds it
            
            eta = monomials(BCVARvec,degBCm:degBCM);%for each element in bg it corresponds a different degree on the expression. It is performed not to exceed the maximum degree in h
 
            [prog,etaKi] = sospolyvar(prog,eta);
            etaK(i) = etaKi;
        end
        
        
    else
        
        etaK = mpvar('etaK',[size(bg,1),1]);
        for i = 1:size(bg,1)
            expri = bg(i);
            % Collect expri(x,c) = g0(c)+g(c)*ht(x) where x are the
            % vars in vartable and c are the decision variables.
            [g0,g,ht] = collect(expri,setdiff(expri.varname,prog.vartable));
            if ~isequal(g0,0)
                ht = [1; ht];
            end
            % Check the max/min deg
            mondegmax = ht.maxdeg;
            mondegmin = ht.mindeg;

            degBCM = degzM- mondegmax;%writes the maximum degre to be added to each BC constraint
            degBCm = degzm- mondegmin;%writes the minimum degre to be added to each BC constraint
  
            eta = monomials(BCVARvec,degBCm:degBCM);%for each element in bg it corresponds a different degree on the expression. It is performed not to exceed the maximum degree in h
 
            [prog,etaKi] = sospolyvar(prog,eta);
            etaK(i) = etaKi;
        
        end
    
    end
    
    %4c)completes the expression
    BV_expr = boundary_values + etaK.'*bg;
    
    %4d)and sets the inequality

    if isa(BV_expr,'polynomial')
        if BV_expr.coefficient~=0
            prog = sosineq(prog,-BV_expr);
        end
    else
        if BV_expr~=0
            prog = sosineq(prog,-BV_expr);
        end
    end
    hb = -BV_expr;
end