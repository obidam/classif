function [l p] = plot_GMMellipse(GMM,idref,id,k)
% plot_GMMellipse H1LINE
%
% [ax el] = plot_GMMellipse(mix,idref,id,k)
%
% Inputs:
%
% Outputs:
%
% Eg:
%
% See Also: 
%
% Copyright (c) 2014, Guillaume Maze (Ifremer, Laboratoire de Physique des Oceans).
% For more information, see the http://codes.guillaumemaze.org
% Created: 2014-12-10 (G. Maze)

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 	* Redistributions of source code must retain the above copyright notice, this list of 
% 	conditions and the following disclaimer.
% 	* Redistributions in binary form must reproduce the above copyright notice, this list 
% 	of conditions and the following disclaimer in the documentation and/or other materials 
% 	provided with the distribution.
% 	* Neither the name of the Ifremer, Laboratoire de Physique des Oceans nor the names of its contributors may be used 
%	to endorse or promote products derived from this software without specific prior 
%	written permission.
%
% THIS SOFTWARE IS PROVIDED BY Guillaume Maze ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, 
% INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL Guillaume Maze BE LIABLE FOR ANY 
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
% LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR 
% BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
% STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

switch GMM.covar_type
	case 'full',      [v,d] = eig(GMM.covars([idref id],[idref id],k));
	case {'diag','mix'}, [v,d] = eig(diag(GMM.covars(k,[idref id])));
	case 'spherical', [v,d] = eig(diag([1 1]*GMM.covars(k)));
end% switch

ALPH = 'ABCDEFG';
for idir = 1 : 2	
	% Ensure that eigenvector has unit length
	v(:,idir) = v(:,idir)/norm(v(:,idir));
	start = GMM.centres(k,[idref id])-sqrt(d(idir,idir))*(v(:,idir)');
	endpt = GMM.centres(k,[idref id])+sqrt(d(idir,idir))*(v(:,idir)');
	if 1
		col = 'kk';
		linex = [start(1) endpt(1)];
		liney = [start(2) endpt(2)];
		l(idir) = line(linex, liney, 'Color', col(idir), 'LineWidth', 2);
	else
		col = 'rb';
		linex = [start(1) endpt(1)];
		liney = [start(2) endpt(2)];
		l(idir) = line(linex, liney, 'Color', col(idir), 'LineWidth', 2);
		text(GMM.centres(k,idref),GMM.centres(k,id),'C','Color', 'k','fontsize',14);
		text(start(1),start(2),[ALPH(idir) '-1'],'Color', col(idir),'fontsize',14);
		text(endpt(1),endpt(2),[ALPH(idir) '-2'],'Color', col(idir),'fontsize',14);
	end% if 
end

% Plot ellipses of one standard deviation
theta = 0:0.02:2*pi;
x = sqrt(d(1,1))*cos(theta);
y = sqrt(d(2,2))*sin(theta);
% Rotate ellipse axes
ellipse = (v*([x; y]))';
ellipse2 = (v*([x; y]))';
% Adjust centre
ellipse  = ellipse + ones(length(theta), 1)*GMM.centres(k,[idref id]);
ellipse2 = 2*ellipse2 + ones(length(theta), 1)*GMM.centres(k,[idref id]);
p(1) = plot(ellipse(:,1), ellipse(:,2), 'k-','LineWidth', 2);
p(2) = plot(ellipse2(:,1), ellipse2(:,2), 'k-','LineWidth', 1);

%in1s = inpolygon(OBS(idref,labels==k),OBS(id,labels==k),ellipse(:,1), ellipse(:,2));
%in2s = inpolygon(OBS(idref,labels==k),OBS(id,labels==k),ellipse2(:,1), ellipse2(:,2));
%leg(1) = {sprintf('1 std (%0.1f%%)',length(find(in1s==1))*100/Klen(ik))};
%leg(2) = {sprintf('2 std (%0.1f%%)',length(find(in2s==1))*100/Klen(ik))};
%legend(p,leg,'location','best');

end %functionplot_GMMellipse