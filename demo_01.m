%DEF GMM simple demo
%REQ Netlab
%
% This demo will create a 2-dimensional training set from 2 random Gaussian distributions and
% - plot the training set PDF
% - train a GMM with ncentres classes and different covariance form (3)
% - plot posteriors (proba for each data to belong to each of the class)
% - plot model PDF
% - plot all the above using a new dataset (a gridded version of the plan)
%
% You may want tro try:
% - to change the number of classes (ncentres)
% - to change the number of points in the training dataset (ndata)
% - to change the demo training set (from 2 to 4)
% 
% 
% 
% Copyright (c) 2016, Guillaume Maze (Ifremer, Laboratoire de Physique des Oceans).
% For more information, see the http://codes.guillaumemaze.org
% Created: 2016-03-14 (G. Maze)

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
clear
addpath('~/matlab/Netlab3.3');

%- Define the number of class to use in the GMM
ncentres = 2; % this is usually called K

%- Create a 2 dimensional dataset:
randn('state', 42);
rand('state', 42);

% Total nb of point in the training set:
ndata = 2000; 

% Dimension of the dataset
input_dim = 2; 

% Select the demo training set:
iset = 1; % 2 equi-probable classes with nearly equal covariance
% iset = 2; % 2 clearly distinct classes with anisotrope covariances
% iset = 3; % 2 overlapping classes with anisotrope covariances
% iset = 4; % another example of 2 clearly distinct classes with anisotrope covariances
% iset = 5; % 3 distinct classes, one having a much smaller covariance than the other 2

% Generate dataset:
[data, datacentres, datapriors, datastandev] = dem2ddat(ndata,iset); 

% Create a Netlab structure with theoretical values:
%mix = gmm(input_dim, size(datacentres,1), 'full'); % Init the classification model with Netlab

% Create a regular 2-dimensional space for "continuous" mapping (ie a new dataset for prediction using trainged GMM)
axl = fix(abs(xtrm(data(:)))+1); % for axis limits
x   = linspace(-axl,axl,100); y = x; % axis
[Xm, Ym] = meshgrid(x,y); % plan
X = Xm(:); Y = Ym(:);

%- Plot the dataset and its main variance axis (EOF)
figure; hold on
plot(data(:, 1), data(:, 2), 'k.','markersize',10,'color',[1 1 1]/2)
grid on, box on, axis square,axis([-1 1 -1 1]*axl),set(gca,'xtick',-axl:axl,'ytick',-axl:axl);
% Plot eigen vectors, main variance axis of the dataset:
[v,d] = eig(data' * data);
for idir = 1 : 2
	% Ensure that eigenvector has unit length
	v(:,idir) = v(:,idir)/norm(v(:,idir));
	start = mean(data)-sqrt(d(idir,idir))*(v(:,idir)');
	endpt = mean(data)+sqrt(d(idir,idir))*(v(:,idir)');
	linex = [start(1) endpt(1)];
	liney = [start(2) endpt(2)];
	l(idir) = line(linex, liney, 'Color', [1 1 1]/2, 'LineWidth', 1,'linestyle','--');
end
title(sprintf('The collection of 2D data\nwith covariance main axis (eofs)'),'fontweight','bold','fontsize',14)

%- Plot observed PDF:
histo = bin2mat(data(:, 1), data(:, 2), ones([ndata 1]), Xm, Ym, @sum); % Histogram
histo(isnan(histo)) = 0;
dspdf = histo./(sum(histo(:)*diff(x(1:2))^2)); % Dataset PDF (2D integral goes to 1)
disp('This should be 1:')
	sum(sum(dspdf*diff(x(1:2)),1)*diff(x(1:2)),2)
dspdf(dspdf==0) = NaN;

figure
pcolor(x,y,dspdf);
shading flat;colorbar
grid on, box on, axis square,axis([-1 1 -1 1]*axl),set(gca,'xtick',-axl:axl,'ytick',-axl:axl);
title(sprintf('Observed PDF'),'fontweight','bold','fontsize',14)

%- Now we classify data using different covariance matrix form
ctype = {'spherical','diag','full'};
ctype = {'full'};
for ii = 1 : length(ctype)
	covar_type = ctype{ii};
	%-- Set up and train the mixture model
	mix = gmm(input_dim, ncentres, covar_type); % Init the classification model with Netlab

	% Set options to init the GMM with a Kmean:
	options = foptions;
	options(14) = 5;	% Just use 5 iterations of k-means in initialisation
	% Initialise the model parameters from the data
	mix = gmminit(mix, data, options);

	% Set options to train the GMM:
	options = foptions;
	options(1)  = 0;		% Prints out error values.
	options(14) = 50;		% Max. Number of iterations.
	[mix, options, errlog] = gmmem(mix, data, options); % E-M algorithm
	% Read the covariance matrix from the GMM output Netlab structure in the form [input_dim,input_dim,ncentres]
	clear covar
	switch covar_type
		case 'full',      
			covar = mix.covars; % Already [input_dim,input_dim,ncentres]
			covarlab = 'Plain';
		case 'diag',      
			for k = 1 : ncentres
				covar(:,:,k) = diag(mix.covars(k,1:input_dim));
			end% for k
			covarlab = 'Diagonal';
		case 'spherical', 
			for k = 1 : ncentres
				covar(:,:,k) = diag([1 1]*mix.covars(k));
			end% for k
			covarlab = 'Isotrope';			
	end% switch
	covar
	
	%-- Plot classification results (GMM 1 and 2 std contours for each class)
	figure; hold on
	plot(data(:, 1), data(:, 2), 'k.','markersize',10,'color',[1 1 1]/2)
	grid on, box on, axis square,axis([-1 1 -1 1]*axl),
	set(gca,'xtick',-axl:axl,'ytick',-axl:axl);
	col = hsv(ncentres);
	for k = 1 : ncentres
		[l p] = plot_GMMellipse(mix,1,2,k);
		set([l p],'color',col(k,:))
	end% for
	axis square,axis([-1 1 -1 1]*axl)
	title(sprintf('GMM class variance contours\n%s covariance matrix and K=%i',covarlab,ncentres),'fontweight','bold','fontsize',14)
	
end% for ii


%- Compare GMM class with theoretical values:
if 0
	figure; hold on
	for ki = 1 : ncentres
		[l p] = plot_GMMellipse(mix,1,2,ki);
		set([l p],'color',[1 1 1]*0,'linewidth',2)
	end% for
	grid on, box on, axis square,axis([-1 1 -1 1]*axl),set(gca,'xtick',-axl:axl,'ytick',-axl:axl);
end% if 

%- Compute GMM probabilities and hard classify data (labels)
post = gmmpost(mix,data); % Proba for one data to belong to one class
[~,labels] = max(post,[],2); % hard classification
prob = gmmprob(mix,data); % Dataset GMM model pdf
acti = gmmactiv(mix,data); % Gaussians for each component for each data

%- Plot model PDF with training set:
figure; hold on
s=scatter(data(:, 1), data(:, 2), 20, prob);
set(s,'markerfacecolor','flat','marker','+');
colorbar
grid on, box on, axis square,axis([-1 1 -1 1]*axl),set(gca,'xtick',-axl:axl,'ytick',-axl:axl);
title(sprintf('Model PDF = Mixture P(X)\nTraining dataset'),'fontweight','bold','fontsize',14)

%- Plot model PDF with new data:
% This is the model PDF values on the regular 2-dimensional space
figure
mpdf = gmmprob(mix, [X Y]);
mpdf = reshape(mpdf, length(x), length(y));
disp('This should be 1:')
	sum(mpdf(:)*diff(x(1:2))*diff(x(1:2)))
pcolor(x,y,mpdf);
shading flat;colorbar
grid on, box on, axis square,axis([-1 1 -1 1]*axl),set(gca,'xtick',-axl:axl,'ytick',-axl:axl);
title(sprintf('Model PDF = Mixture P(X)\nPrediction on new data (gridded plan)'),'fontweight','bold','fontsize',14)

%- Compare model with observed PDF:
figure; hold on
pcolor(x,y,mpdf./dspdf);
shading flat;colorbar
grid on, box on, axis square,axis([-1 1 -1 1]*axl),set(gca,'xtick',-axl:axl,'ytick',-axl:axl);
title(sprintf('Model vs observed PDF (ratio)\nIncrease size of training set to improve this'),'fontweight','bold','fontsize',14)

%- Superimpose GMM components ellipses on posteriors (component probabilities)
for k = 1 : ncentres	
	figure; hold on
	s = scatter(data(:, 1), data(:, 2), 20, post(:,k));
	set(s,'markerfacecolor','flat','marker','+');
	colorbar
	for ki = 1 : ncentres
		[l p] = plot_GMMellipse(mix,1,2,ki);
		if ki == k
			set([l p],'color',[1 1 1]/2,'linewidth',2)
		else
			set([l p],'color',[1 1 1]/2,'linestyle','--')
		end% if 
	end% for
	grid on, box on, axis square,axis([-1 1 -1 1]*axl),set(gca,'xtick',-axl:axl,'ytick',-axl:axl);
	title(sprintf('Component posteriors P(K_%i|X)\nTraining dataset',k),'fontweight','bold','fontsize',14)	
end% for k

%- Plot labels (hard classification)
figure;hold on
colormap(hsv(ncentres));
s=scatter(data(:, 1), data(:, 2), 20, labels);
set(s,'markerfacecolor','flat','marker','+');
cl = colorbar;
clabelcmap(cl,[1 ncentres],ncentres,1,' Class %0.0f');
grid on, box on, axis square,axis([-1 1 -1 1]*axl),set(gca,'xtick',-axl:axl,'ytick',-axl:axl);
title(sprintf('Data classification\nTraining dataset'),'fontweight','bold','fontsize',14)

%- Plot GMM probabilities and hard classification (labels) on regular 2-dimensional space 
Z = gmmpost(mix, [X Y]);
Z = reshape(Z, length(x), length(y), ncentres);
for k = 1 : ncentres	
	figure; hold on
	pcolor(x,y,squeeze(Z(:,:,k)));
	shading flat;colorbar
	%plot(data(:, 1), data(:, 2), 'k.','markersize',6,'color',[1 1 1]/2)
	for ki = 1 : ncentres
		[l p] = plot_GMMellipse(mix,1,2,ki);
		if ki == k
			set([l p],'color',[1 1 1],'linewidth',2)
		else
			set([l p],'color',[1 1 1],'linestyle','--')
		end% if 
	end% for
	grid on, box on, axis square,axis([-1 1 -1 1]*axl),set(gca,'xtick',-axl:axl,'ytick',-axl:axl);
	title(sprintf('Component posteriors P(K_%i|X)\nPrediction on new data (gridded plan)',k),'fontweight','bold','fontsize',14)	
end% for k

figure; hold on
Z = gmmpost(mix, [X Y]);
Z = reshape(Z, length(x), length(y), ncentres);
[~,Z] = max(Z,[],3);
pcolor(x,y,Z);
shading flat;colorbar
grid on, box on, axis square,axis([-1 1 -1 1]*axl),set(gca,'xtick',-axl:axl,'ytick',-axl:axl);
title(sprintf('Data labels\nPrediction on new data (gridded plan)'),'fontweight','bold','fontsize',14)

%- Plot label metric for training set (see Maze et al, POC, 2016)
Plist = [0 0.33 0.66 0.9 .99 1];
rowl0 = {'Unlikely','As likely as not','Likely','Very Likely','Virtually certain'};
cmap2 = hsv(length(Plist)-1); cmap2 = cmap2([1 5 4 2 3],:);

post = gmmpost(mix,data); % Proba for one data to belong to each class
maxP = max(post,[],2); % Max posterior
R    = (maxP-1/ncentres)*ncentres/(ncentres-1); % Scale it

% Create colorscale for each item
c  = R;
cc = [];
for ir = 1 : length(Plist)-1
	ic = find(c>=Plist(ir)&c<=Plist(ir+1));
	cc(ic,:) = repmat(cmap2(ir,:),[length(ic) 1]);
end% for ip

figure; hold on
colormap(cmap2);
s=scatter(data(:, 1), data(:, 2), 20, cc);
set(s,'markerfacecolor','flat','marker','+');
cl = colorbar; 
caxis([0 1]);
clabelcmap(cl,[0 1],size(cmap2,1),1,'%0.2f');		
set(cl,'yticklabel',rowl0,'fontsize',12,'fontweight','bold');
grid on, box on, axis square,axis([-1 1 -1 1]*axl),set(gca,'xtick',-axl:axl,'ytick',-axl:axl);
for ki = 1 : ncentres
	[l p] = plot_GMMellipse(mix,1,2,ki);
	set([l p],'color',[1 1 1]*0,'linewidth',2)
end% for
title(sprintf('Label metric\nTraining dataset'),'fontweight','bold','fontsize',14)


%- Plot label metric for gridded plan
POST = gmmpost(mix, [X Y]);
maxP = max(POST,[],2);
R    = (maxP-1/ncentres)*ncentres/(ncentres-1);
% Create colorscale for each item
c = R;
cc = [];
for ir = 1 : length(Plist)-1
	ic = find(c>=Plist(ir)&c<=Plist(ir+1));
	cc(ic,:) = repmat(cmap2(ir,:),[length(ic) 1]);
end% for ip
R  = reshape(R, length(x), length(y));
cc = reshape(cc, length(x), length(y), 3);

figure; hold on
colormap(cmap2);
surf(x,y,R,cc);
view(2);shading flat;
cl = colorbar; 
caxis([0 1]);
clabelcmap(cl,[0 1],size(cmap2,1),1,'%0.2f');		
set(cl,'yticklabel',rowl0,'fontsize',12,'fontweight','bold');
grid on, box on, axis square,axis([-1 1 -1 1]*axl),set(gca,'xtick',-axl:axl,'ytick',-axl:axl);
for ki = 1 : ncentres
	[l p] = plot_GMMellipse(mix,1,2,ki);
	set([l p],'color',[1 1 1]*0,'linewidth',2)
	for i=1:2;set(l(i),'zdata',ones(size(get(l(i),'ydata'))));end
	for i=1:2;set(p(i),'zdata',ones(size(get(p(i),'ydata'))));end
end% for
title(sprintf('Label metric\nPrediction on new data (gridded plan)'),'fontweight','bold','fontsize',14)








