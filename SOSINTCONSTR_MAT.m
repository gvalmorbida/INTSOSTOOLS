function  [F,Gst,Hst,Z,IPMON] = SOSINTCONSTR_MAT(expr,Dvar)

%--------------------------------------------------------------------------
%SOSINTCONSTR_MAT.m
%
%This function builds the matrix in the quadratic-like form in the
%integrand, including integration by parts and the gram matrix
%
%inputs - expr: the integrand, it must be a function of the dependent
%               variables DepVar and the independent variables 
%
%outputs - F,G,H: the matrices defining the complete-form integrand
%          Z:     the vector of monomials defining the quadratic form
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


[n,~] = size(Dvar);%obtains the number of variables and the order

%==========================================================================
%STEP 1)
%EXTRACT THE LIST OF MONOMIALS ON THE DEPENDENT VARIABLE IN THE EXPRESSION
%==========================================================================
%WE HAVE THE LIST AS THE POWERS OF ENTRIES OF VECTOR
%vec(Dvar)' = [u1' u1_x' u1_xx' ... u1_ord' ... un' un_x' un_xx' ... un_ord'];
[coeffs,Zf,L,Lo] = GETINTEQUATION(expr,Dvar);%Lt and Lo have the same number of rows
%the rows are ordered
%monlistexpr = mysympower(vec(Dvar)',unique(Zf,'rows'));%WARNING!!! MYSYMPOWER DOES NOT RETURN A VECTOR OF MONOMIALS IN THE SAME ORDER AS THE VECTOR Zf!!

%the following lines extract the 'deepest' row in Lo corresponding to each
%of the present monomials in Lt.
%it relies on the fact that the rows are ordered as we are using the option 'stable' as an input
[Lmax,ind] = unique(L,'rows','stable');%unique returns the index of the last element repeated on the array
Lomax=sparse(length(ind),size(Lo,2));%allocates space for Lomax
for i = 1:length(ind)-1
    [~,indLomax] = max(sum(Lo(ind(i):ind(i+1)-1,1:sum(Lmax(i,:),2)),2));%the sum over Lu(i,:) defines the degree, therefore the number of elements in Lo to be in the sum
    Lomax(i,:) = Lo((ind(i)-1)+indLomax,:); %needs an initial index since extracts from the original Lo
end
[~,indLomax] = max(sum(Lo(ind(length(ind)):size(Lo,1),1:sum(Lmax(length(ind),:),2)),2));%obtains the last element (Lo(1: ... )
if length(ind)==1
    Lomax(end,:) = Lo(indLomax,:);
else
    Lomax(end,:) = Lo((ind(length(ind))-1)+indLomax,:);
end
% Lt and Lo - contain all monomials
% Lmax and Lomax - contain only the deepest monomials of each degree
[Lmax,Lomax] = INCLUDEEVEN(Lmax,Lomax);%includes the terms of even degree, higher than the odd degree, to be used in the computation of the monomials for integration by parts
[Lf,Lof,LfTREE,LofTREE] = IPMONOMIALS_V2(Lmax,Lomax);%this is the command to include the monomials from integration by parts.
IPMON.Lf = Lf; IPMON.Lof = Lof;
%==========================================================================
%STEP 2)
%SETS THE MONOMIAL VECTOR FORMING THE QUADRATIC EXPRESSION
%==========================================================================
LtQ = unique([LfTREE LofTREE;Lmax Lomax],'rows','stable');
LQ = LtQ(:,1:n);LoQ = LtQ(:,n+1:end);%the matrix with the degree and the order without repeated elements
monforQuad = vec(MONLISTDER(Dvar,LQ,LoQ));%uses the function MONLISTDER (originally written for computing monomials to enter the derivative)
[monlistQuad] = MONFORQUADEXPR_V2(monforQuad,Dvar);

%EXTRACTS THE LIST OF MONOMIALS
Z = [];
if isa(monlistQuad,'polynomial')
    cvartable =  char(Dvar.varname);
    [~,nvar] = size(monlistQuad.degmat);
    for i = 1:size(monlistQuad,1)
        [g0,g,h] = collect(monlistQuad(i),setdiff(monlistQuad.varname,cvartable));
        g = g(:);
        if ~isequal(g0,0)
            g = [g0;g];
            h = [1; h];
        end
        [nmon,~] = size(h.degmat);
        
        % Reorder the monomials matrix with variables in the order
        % listed in cvartable and sorted monomials
        Zt = zeros(nmon,nvar);
        [~,idx]=ismember(h.varname,cvartable);
        Zt(:,idx) = full(h.degmat);
        Z(i,:) = double(Zt);
    end
    
else
    
    charvartable = converttochar(vec(Dvar)');%obtains the string of characters to be passed to the symbolic engine
    coefmon = feval(symengine,'poly2list',sum(monlistQuad),charvartable);%calls the symbolic engine to obtain the list of monomials on the expression
    coefmon = coefmon.';%this is the list of coefficients and exponents
    for i = 1:size(coefmon,1)
        dummyvar = reshape(coefmon(i,:),2,1);
        Z(i,:) = double(dummyvar(2));
    end
    
end
Z = sparse(Z);


%THE FOLLOWING IS NOT REQUIRED AND CAN BE REMOVED. ONLY HELPS VISUALIZING THE RESULTS
[~,ind] = sort(sum(Z,2));
Z = Z(ind,:);

%OBTAINS THE MATRIX A, TO BE USED IN THE CONSTRUCTION OF THE QUADRATIC EXPRESSION
[A,ZZ] = getconstraint(Z);%ZZ contains the monomials of Z kron Z, and A is the matrix combining the elements of vec(Q) as A*vec(Q)=F' and we have Z'QZ = F*ZZ
[rnz] = find(sum(A,2));
A = A(rnz,:);

%==========================================================================
%STEP 3)
%COMPUTE A QUADRATIC REPRESENTATION OF THE INPUT EXPRESSION
%==========================================================================
[~,i1,i2] = intersect(Zf,ZZ,'rows');%finds the common rows in Zf in ZZ to establish the relation between the two and set the representation of F in the extended vector
M_ZZtoZf = sparse(i1,i2,ones(size(Zf,1),1),size(Zf,1),size(ZZ,1));%this matrix transforms the original expression into
fb = (coeffs.'*M_ZZtoZf).';%computes the vector of coefficients multiplying the vector of monomials ZZ

iA = pinv(full(A));
Fmat = mat(iA*fb);
if isa(expr,'polynomial')
    %THERE IS NO LEFTMULT FUNCTION FOR fb
    [U,S,V] = svd(full(A));%computes the matrix yielding the quadratic expression
    Fmat = mat(V*(S\(U'))*fb);
else
    [U,S,V] = svd(full(A));%computes the matrix yielding the quadratic expression
    Fmat = mat(V*(S\(U'))*fb);
    %Fmat = mat(A\fb);%computes the matrix yielding the quadratic expression
end
F = 1/2*(Fmat+Fmat.');%F has the dimensions of the monlistQuad vector.


%==========================================================================
%STEP 4)
%COMPUTES THE QUADRATIC REPRESENTATION FOR THE GRAM MATRIX
%==========================================================================
dimG = size(Z,1)^2;
[INDR,INDC] = ind2sub(sqrt(dimG),[1:dimG]);%extracts the indices
REMOVE = find(INDR>INDC);%select the indices of symmetric elements to remove
RR = INDR(REMOVE);%the indices of the rows to remove
CR = INDC(REMOVE);%the indices of the columns to keep
indADDsource = sub2ind([sqrt(dimG) sqrt(dimG)],RR,CR);%elements to be removed
indADDdest = sub2ind([sqrt(dimG) sqrt(dimG)],CR,RR);%elements to be kept
Aman = A;
Aman(:,indADDdest) =  Aman(:,indADDdest) + Aman(:,indADDsource);%adds the equivalent terms
[indKEEP,~] = setdiff([1:dimG],indADDsource,'stable');%this is the index of the elements we kept
Aman = Aman(:,indKEEP);%updates A only with the elements to keep
An = null(full(Aman));%NOT EFFICIENT! BUT NEEDED TO AVOID THE EQUALITY CONSTRAINTS AND HELPS TO OBTAIN THE MINIMAL NUMBER OF VARIABLES
An(abs(An)<1e-12)=0;
[~,IND2] = size(An);

if isa(expr,'polynomial')
    g = mpvar('g',[IND2,1]); 
else
    g = sym('g',[IND2,1]); g = sym(g,'real');%generates the gram matrix variables
end
Ag = An*g;%multiplies by the variables

    if isa(expr,'polynomial')
        Adg = polynomial(zeros(1,dimG));
    else
        Adg = sym(zeros(1,dimG));
    end


Adg(indKEEP) = Ag;
Adg(:,indADDsource) = Adg(:,indADDdest);
G = mat(Adg);

Gst.g = g;
Gst.G = G;

%==========================================================================
%STEP 5)
%COMPUTE THE QUADRATIC REPRESENTATION FOR THE INTEGRATION BY PARTS
%==========================================================================
if ~isempty(Lf)%if Lf is empty, then the expression contains only terms of order zero
    
    monIP = MONLISTDER(Dvar,Lf,Lof);%obtains the the monomials to be included in the integration by parts
    monIPdiff = MONLISTDERDUDX(Dvar,Lf,Lof);%obtains the derivative of monIP
    %it is important to have the each entry on monIPdiff correspond to the same entry of monIP
    %==========================================================================
    %STEP 5.1)
    %CREATE A VECTOR OF VARIABLES
    %==========================================================================
    n_polIP = size(monIP,1);%this is the number of functions to be introduced in the integration by parts
    
    if isa(expr,'polynomial')
        h = mpvar('h',[n_polIP,1]); 
    else
        h = sym('h',[n_polIP,1]); h = sym(h,'real');%define the variables to be introduced by integration by parts
    end
    
    
    %==========================================================================
    %STEP 5.2)
    %COMPUTES THE EXPRESSIONS IN THE KERNEL COMING FROM THE REPRESENTATION OF THE CURL
    %==========================================================================
    %EXTRACTS THE LIST OF MONOMIALS IN monIP
    Zh = []; %the loop below guarantees that the monomials are in the order corresponding to the vector h
    M_ZZtoZh =sparse(size(monIP,1),size(ZZ,1));
    
    for j = 1:length(monIP)%for each element in monIP
        
        if  isa(expr,'polynomial')
            nvar = size(ZZ,2);
            [g0,g,ht] = collect(monIP(j),setdiff(monIP(j).varname,cvartable));
            g = g(:);
            if ~isequal(g0,0)
                g = [g0;g];
                ht = [1; ht];
            end
            [nmon,~] = size(ht.degmat);
            Zh = zeros(nmon,nvar);
            [~,idx]=ismember(ht.varname,cvartable);
            Zh(:,idx) = full(ht.degmat);
            Zh = sparse(Zh);
        else
            coefmon = feval(symengine,'poly2list',monIP(j),charvartable).';%calls the symbolic engine to obtain the list of monomials on the expression
            dummyvar = reshape(coefmon,2,1);
            Zh = double(dummyvar(2));
        end
        
        [~,~,i2] = intersect(Zh,ZZ,'rows');%finds the common rows in Zf in ZZ to establish the relation between the two and set the representation of F in the extended vector
        M_ZZtoZh(j,i2) =  1;%this matrix transforms the original expression into the h representation
                            %it is indexed on j since in this case it is performe row-wise
    end
    hb = (h.'*M_ZZtoZh).';%computes the vector of coefficients multiplying the vector of monomials ZZ
    
    Hmat = mat(iA*hb);
    %Hmat = mat(V*(S\(U'*hb)));
    
    if isa(expr,'polynomial')
        %THERE IS NO LEFTMULT FUNCTION FOR fb
        %[U,S,V] = svd(full(A));%computes the matrix yielding the quadratic expression
        Hmat = mat(V*(S\(U'))*hb);
    else
        Hmat = mat(V*(S\(U'))*hb);
        %Hmat = mat(A\hb);%computes the matrix yielding the quadratic expression (converts to a square matrix the vector of n^2 entries)
    end
    
    
    H = 1/2*(Hmat+Hmat.');%F has the dimensions of the monlistQuad vector.
    
    %EXTRACTS THE LIST OF MONOMIALS IN monIPdiff
    Zhdiff = [];
    M_ZZtoZhdiff =sparse(size(monIPdiff,1),size(ZZ,1));
    for j = 1:length(monIP)%since the derivative may have more than one monomial, needs to do for each element in monIP
        
        if  isa(expr,'polynomial')
            nvar = size(ZZ,2);
            [g0,g,ht] = collect(monIPdiff(j),setdiff(monIPdiff(j).varname,cvartable));
            g = g(:);
            if ~isequal(g0,0)
                g = [g0;g];
                ht = [1; ht];
            end
            [nmon,~] = size(ht.degmat);
            Zhdiff = zeros(nmon,nvar);
            [~,idx]=ismember(ht.varname,cvartable);
            Zhdiff(:,idx) = full(ht.degmat);
            Zhdiff = sparse(Zhdiff);
            [~,~,i2] = intersect(Zhdiff,ZZ,'rows');%finds the common rows in Zf in ZZ to establish the relation between the two and set the representation of F in the extended vector
            M_ZZtoZhdiff(j,i2) =  double(g);%this matrix transforms the original expression into the diff_h representation
        else
            coefmon = feval(symengine,'poly2list',monIPdiff(j),charvartable).';%calls the symbolic engine to obtain the list of monomials on the expression
            for i = 1:size(coefmon,1)%for each monomial appearing in the derivative of monIP(j)
                dummyvar = reshape(coefmon(i,:),2,1);
                Zhdiff(i,:) = double(dummyvar(2));%KEEPS THE PREVIOUS ROWS OF Zhdiff EVEN IF NOT NECESSARY
                [~,~,i2] = intersect(Zhdiff(i,:),ZZ,'rows');%finds the common rows in Zf in ZZ to establish the relation between the two and set the representation of F in the extended vector
                M_ZZtoZhdiff(j,i2) =  double(dummyvar(1));%this matrix transforms the original expression into the diff_h representation
                %it is indexed on j since in this case it is performe row-wise
            end
        end
        
        
        
    end
    
    
    hbdiff = (h.'*M_ZZtoZhdiff).';%computes the vector of coefficients multiplying the vector of monomials ZZ
    
    Hmatdiff = mat(iA*hbdiff);
    %Hmatdiff = mat(V*(S\(U'*hbdiff)));
    
    if isa(expr,'polynomial')
        Hmatdiff = mat(V*(S\(U'))*hbdiff);
    else
        Hmatdiff = mat(V*(S\(U'))*hbdiff);
        %Hmatdiff = mat(A\hbdiff);%computes the matrix yielding the quadratic expression
    end
    
    
    Hdiff = 1/2*(Hmatdiff+Hmatdiff.');%F has the dimensions of the monlistQuad vector.
    
    %==========================================================================
    %STEP 5.3)
    %SETS THE INTEGRATION BY PARTS MATRICES
    %==========================================================================
    
    Hst.h = h;%defines the list of variables
    %Hst.monIP = Zh;%defines the list of dependent variables which are related to each entry of h
    Hst.H =  Hdiff;%WARNING!!! THIS IS THE CORRECT WAY, Hst.H is defined by matrix Hdiff %defines the matrix with functions h
    Hst.Hdiff =  H;%WARNING!!! THIS IS THE CORRECT WAY, Hst.Hdiff is defined by matrix H %defines the matrix with derivatives of h
    
else%there is no spatial derivative in the expression, therefore no H matrices
    Hst.h = [];
    Hst.monIP = [];
    Hst.H = zeros(size(monlistQuad),size(monlistQuad));
    Hst.Hdiff = zeros(size(monlistQuad),size(monlistQuad));
end