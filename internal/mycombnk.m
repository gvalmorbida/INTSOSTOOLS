function M  = mycombnk(ord,ncol)

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



%the funciton avoids the call to combnk (matlab statistics toolbox function).
%it gives an ordered set of rows. 

M0 = (0:ord)';
M = M0;
for i = 2:ncol
    M = [repmat(M,[ord+1 1])  vec(repmat(M0', [size(M,1) 1]))];
    IND = find(M(:,i-1)<=M(:,i));
    M = M(IND,:);
end
%mat = sparse(mat);


