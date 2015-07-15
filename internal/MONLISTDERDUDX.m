function  [monlist] = MONLISTDERDUDX(DepVar,L,Lo)
%--------------------------------------------------------------------------
%MONLISTDERDUDX.m
%
%This function takes a list of monomials on a variable and its (spatial)
%derivatives and returns a vector containing the DERIVATIVES OF THE monomials 
%as symbolic variables. These are taken from the list of monomials Depbar
%
%THIS IS THE CODE FOR ONE SPATIAL VARIABLE
%
%inputs - DepVar: DepVar contains all the dependent variables and its derivatives 
%         L:      matrix containing the list of monomials that appear on the
%                 base monomials, and the indexes for the related monomials
%                 with the spatial derivatives in matrix Lo. 
%                 L in R^(nmb X (n+2))
%         Lf:     matrix containing the list of spatial derivatives on each
%                 monomial.
%                 Lo in R^(nm X dg)
%         idx:    Index list relating rows of L to blocks of rows in Lf
%
%outputs - monlist: The list containing the derivatives of the monomials
%                   described by L, Lf, idx
%
%example:
%           DepVar = [u1 u1x u1xx; u2 u2x u2xx];
%           L = [1 1; 0 2]; u1u2, u2^2
%           idx = [1 2; 3 3];    the first two are related to u1u2, the third is
%                               related to u2^2
%           Lf = [0 2; 1 1; 1 1];the derivative for u1u2xx, u1xu2x, u2x^2
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






degL = sum(L,2);
degmax = max(degL);

if isa(DepVar,'polynomial')
    monlist = polynomial(zeros([size(Lo,1),1]));
else
    monlist = sym(zeros([size(Lo,1),1]));
end



for i = 1:degmax
    
    if isa(DepVar,'polynomial')
        monlistbase = polynomial(zeros([size(Lo,1),1]));
    else
        monlistbase = sym(zeros([size(Lo,1),1]));
    end
    

    Lob = Lo;
    Lob(:,i) = Lob(:,i)+1;%adds one to the ith column of Lfb, this is differentiating
    
    setind = find(gt(degL,i-1));%only need to differentiate the ith entry if the degree (strictly) greater than i-1
    %Lin = L(setind,:);%just uses the monomials of degree (strictly) higher than i-1
    %idxin = idx(setind,:);
    
    %monlistbase(min(idxin):max(idxin)) = MONLISTDER_MOD(DepVar,Lin,Lfb,idxin);%builds the list of monomials using MONLISTDER
    monlistbase(setind)  = MONLISTDER(DepVar,L(setind,:),Lob(setind,:));%builds the list of monomials using MONLISTDER
    
    monlist = monlist + monlistbase;%updates the derivative expression
end