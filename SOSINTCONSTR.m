function  [exprh,Zd,IPMON] = SOSINTCONSTR(expr,Dvar)

%--------------------------------------------------------------------------
% SOSINTCONSTR_V2.m
% 
% This function extracts the list of monomials on variables Dvar of the
% integrand expr. It returns the list of monomials to be introduced in the
% integration by parts and a structure containing the expressions of the
% monomials to be introduced by the FTC.
% 
% THIS IS THE CODE FOR ONE SPATIAL VARIABLE AND THE [0,1] DOMAIN
% 
% inputs - expr: the integrand, it must be a function of the dependent 
%               variables DepVar and the independent variables IndepVar
%               
%         Dvar: the matrix containing the dependent variables
%
% outputs - exprh:   structure containing the fields
%                       .hvar - the list of functions
%                       .h -    a polynomial on the dependent variable Dvar which
%                               coefficients are h
%                       .hdiff -	a polynomial on the dependent variable Dvar which
%                                   coefficients are h which must receive the derivative of
%                                   hvar
%           Zd:     the list of the monomials in expr
%           IPMON:  structure describing the polynomials on Dvar in
%                   expr.hiff, with fields
%                   .L -	matrix containing the degree on each variable
%                   .Lo -	matrix containing the order of each element
%                   each row of IPMON corresponds to an element in exprh.h 
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




[n,ord] = size(Dvar);
ord = ord-1;
%==========================================================================
%STEP 1)
%EXTRACT THE LIST OF MONOMIALS ON THE DEPENDENT VARIABLE IN THE EXPRESSION
%==========================================================================
 

%WE HAVE THE LIST AS THE POWERS OF ENTRIES OF VECTOR
%vec(u_D)' = [u1' u1_x' u1_xx' ... u1_ord' ... un' un_x' un_xx' ... un_ord'];
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

%==========================================================================
%STEP 2)
%EXTRACT THE LIST OF IPMONOMIALS
%==========================================================================
[Lf,Lof,LfTREE,LofTREE] = IPMONOMIALS_V2(Lmax,Lomax);%this is the command to include the monomials from integration by parts.
IPMON.Lf = Lf; IPMON.Lof = Lof;

%==========================================================================
%STEP 3)
%SETS THE MONOMIAL VECTOR 
%==========================================================================
LtQ = unique([LfTREE LofTREE;Lmax Lomax],'rows','stable');
LQ = LtQ(:,1:n);LoQ = LtQ(:,n+1:end);%the matrix with the degree and the order without repeated elements
Zd = vec(MONLISTDER(Dvar,LQ,LoQ));%uses the function MONLISTDER (originally written for computing monomials to enter the derivative)
 

if ~isempty(Lf)%if L is empty, then the expression contains only terms of order zero
    
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
    
    
    exprh.h = monIPdiff.'*h;
    exprh.hdiff = monIP.'*h;
    exprh.hvar = h;
    
else
    exprh.h = 0;
    exprh.hdiff = 0;
    exprh.hvar = 0;
end
 