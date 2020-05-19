%  FILE:
%      write_bathy.m
%   
%  AUTHORS:
%      Sasan Tavakkol <tavakkol@usc.edu> or <sasantavakkol@yahoo.com>
%   
%    LAST UPDATE:
%      26-Aug-2016
% 
%  COPYRIGHT:
%      Copyright (C) 2016, Sasan Tavakkol.
%   
%      While this file is not a part of Celeris software, it is included in
%      the project to assist using the software. This file is also covered
%      under the same license as Celeris.
%   
%      you can redistribute this file
%      and/or modify it under the terms of the GNU General Public
%      License as published by the Free Software Foundation, either
%      version 3 of the License, or (at your option) any later version.
%   
%      Celeris is distributed in the hope that it will be
%      useful, but WITHOUT ANY WARRANTY; without even the implied
%      warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%      See the GNU General Public License for more details.
%   
%      You should have received a copy of the GNU General Public License
%      along with the Shallow Water Demo. If not, see
%      <http://www.gnu.org/licenses/>.
   


function write_bathy (name,nx,ny,B)
    fileID = fopen(name,'w');

    fprintf(fileID,'[nx] %d \n', nx);
    fprintf(fileID,'[ny] %d \n', ny);
    fprintf(fileID,'\n\n\n====================================\n');
    for i=1:ny
        for j=1:nx

            fprintf(fileID,'%0.8f ', B(j,i));

        end
                fprintf(fileID,'\n');
    end

    fclose(fileID);
end