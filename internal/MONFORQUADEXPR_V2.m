function [monlist] = MONFORQUADEXPR_V2(expr,vartable)

%--------------------------------------------------------------------------
%MONFORQUADEXPR.m
%
%This function takes monomials from a list and obtaining the minimal vector
%of monomials to be written as a quadratic expression
%
%THIS IS THE CODE FOR ONE SPATIAL VARIABLE
%
%inputs -           expr:       a polynomial expression in the dependent variables
%                               (LATER VERSIONS MUST ALLOW ALSO TO BE ON THE the 
%                               independent variables and the decision variables)
%                   vartable:   
%
%
%outputs -          monlist:    the list of monomials to build the quadratic
%                               expression
%
%--------------------------------------------------------------------------

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

% 06/07/15 - GV





%1)get the full list of monomials in the expression expr

if isa(vartable,'polynomial')%if it is a pvar
    
    
    cvartable =  char(vartable.varname);
    Zfull = [];
   % [~,nvar] = size(expr.degmat);%gets the number of dependent variables
    [nrvar,ncvar] = size(vartable);%gets the number of dependent variables
    nvar = nrvar*ncvar;
    for i = 1:length(expr)

        [g0,g,h] = collect(expr(i),setdiff(expr(i).varname,cvartable));
        g = g(:);
        if ~isequal(g0,0)
            g = [g0;g];
            h = [1; h];
        end
        [nmon,~] = size(h.degmat);
        
        % Reorder the monomials matrix with variables in the order
        % listed in cvartable and sorted monomials
        Z = zeros(nmon,nvar);
        [~,idx]=ismember(h.varname,cvartable);
        Z(:,idx) = full(h.degmat);
        Z = sparse(Z);
        Zfull = sortNoRepeat(Zfull,Z);
    end
    
else
    charvartable = converttochar([vec(vartable')']);
    Zfull = [];
    for i = 1:length(expr)
        Z0 = [];
        coefmon = feval(symengine,'poly2list',expr(i),charvartable);
        for j = 1:length(coefmon)
            dvar = coefmon(j);
            Z0 = [Z0; double(dvar(2))];
        end
        Zfull = sortNoRepeat(Zfull,Z0);
    end
    
end

numstates = size(Zfull,2);
degZfull = sum(Zfull,2);

Zn = [];
for i = min(degZfull):max(degZfull)
    
    setind = find(degZfull==i);
    
    Zt = Zfull(setind,:);
    
    % Creating extra variable
    maxdeg = full(max(sum(Zt,2)));
    mindeg = full(min(sum(Zt,2)));
    Z = monomials(numstates,[floor(mindeg/2):ceil(maxdeg/2)]);
    
    % Discarding unnecessary monomials
    maxdegree = ceil(max(Zt,[],1)/2);
    mindegree = floor(min(Zt,[],1)/2);
    j = 1;
    while (j <= size(Z,1))
        Zdummy1 = maxdegree-Z(j,:);
        Zdummy2 = Z(j,:)-mindegree;
        idx = find([Zdummy1, Zdummy2]<0);
        if ~isempty(idx)
            Z = [Z(1:j-1,:); Z(j+1:end,:)];
        else
            j = j+1;
        end;
    end;
    Zn = [Zn; Z];
end

if isa(vartable,'sym')
    monlist  = mysympower(vec(vartable')',unique(Zn,'rows'));%Then need to set up the monomials vector for the final expression. 
elseif isa(vartable,'polynomial')
  Zn = unique(Zn,'rows');
  coefftemp = speye(size(Zn,1));
  monlist = polynomial(coefftemp,Zn,vartable.varname,[size(Zn,1),1]);
end;


