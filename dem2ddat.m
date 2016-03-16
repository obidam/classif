function [data, c, prior, sd] = dem2ddat(ndata,config)
%DEM2DDAT Generates two dimensional data for demos (modified version from Netlab dem2ddat).
%
%	[DATA, CENTRES, PRIORS, STANDEV] = DEM2DDAT(NDATA,ISET) also returns a matrix containing the
%	centres of the Gaussian distributions.
%
% iset = 1; % 2 equi-probable classes with equal covariance
% iset = 2; % 2 clearly distinct classes with anisotrope covariances
% iset = 3; % 2 overlapping classes with anisotrope covariances
% iset = 4; % another example of 2 clearly distinct classes with anisotrope covariances
% iset = 5; % 4 distinct classes, one having a much smaller covariance than the other 3
%
% Revised: 2016-03-14 (G. Maze) Added more demo datasets

%	Copyright (c) Ian T Nabney (1996-2001)

input_dim = 2;

% Fix seed for reproducible results
randn('state', 42);

% Generate white data in input_dim dimensional space
data = randn(ndata, input_dim);

switch config
	case 0
		ncentres = 2;
	
		% Priors for the three clusters
		prior(1) = 0.3;
		prior(2) = 0.5;
		prior(3) = 0.2;

		% Cluster centres
		c = [2.0, 3.5; 0.0, 0.0; 0.0, 2.0];

		% Cluster standard deviations
		sd  = [0.2 0.2 1.0];

		% Put first cluster at (2, 3.5)
		data(1:prior(1)*ndata, 1) = data(1:prior(1)*ndata, 1) * 0.2 + c(1,1);
		data(1:prior(1)*ndata, 2) = data(1:prior(1)*ndata, 2) * 0.2 + c(1,2);

		% Leave second cluster at (0,0)
		data((prior(1)*ndata + 1):(prior(2)+prior(1))*ndata, :) = ...
			data((prior(1)*ndata + 1):(prior(2)+prior(1))*ndata, :) * 0.5;

		% Put third cluster at (0,2)
		data((prior(1)+prior(2))*ndata +1:ndata, 2) = ...
			data((prior(1)+prior(2))*ndata+1:ndata, 2) + c(3, 2);
			
	case 1
		ncentres = 2;

		% Priors for the 2 clusters
		prior(1) = 0.5;
		prior(2) = 0.5;

		% Cluster centres
		c = [-2.0, -2.0; 2.0, 2.0];

		% Cluster standard deviations
		sd  = [0.8 0.8];
		covar = [sd(1) 0; 0 sd(2)];

		% Put first cluster at (2, 3.5)
		data(1:prior(1)*ndata, 1) = data(1:prior(1)*ndata, 1) * sd(1) + c(1,1);
		data(1:prior(1)*ndata, 2) = data(1:prior(1)*ndata, 2) * sd(1) + c(1,2);

		% Put second cluster at (0,0)
		data((prior(1)*ndata + 1):(prior(2)+prior(1))*ndata, 1) = ...
			data((prior(1)*ndata + 1):(prior(2)+prior(1))*ndata, 1) * sd(2) + c(2,1);
		data((prior(1)*ndata + 1):(prior(2)+prior(1))*ndata, 2) = ...
			data((prior(1)*ndata + 1):(prior(2)+prior(1))*ndata, 2) * sd(2) + c(2,2);

		
	case 2
		ncentres = 2;

		% Priors for the 2 clusters
		prior(1) = 0.8;
		prior(2) = 0.2;

		% Cluster centres
		c = [0.0, 0.0; 2.0, 3.0];

		% Cluster standard deviations
		sd    = [0.5 0.5];
		covar = [1 .3; .3 1];

		% Put first cluster 
		
		% 
		x = data(1:prior(1)*ndata, 1) * sd(1) + c(1,1);
		y = data(1:prior(1)*ndata, 2) * sd(1) * 2 + c(1,2);
		[th,r] = cart2pol(x,y);
		[x,y]  = pol2cart(th-pi/3,r);
		data(1:prior(1)*ndata, 1) = x;
		data(1:prior(1)*ndata, 2) = y;

		% Leave second cluster at (0,0)
		x = data((prior(1)*ndata + 1):(prior(2)+prior(1))*ndata, 1) * sd(2) + c(2,1);
		y = data((prior(1)*ndata + 1):(prior(2)+prior(1))*ndata, 2) * sd(2) + c(2,2);
		data((prior(1)*ndata + 1):(prior(2)+prior(1))*ndata, 1) = x;	
		data((prior(1)*ndata + 1):(prior(2)+prior(1))*ndata, 2) = y;


	case 20
		ncentres = 2;

		% Priors for the 2 clusters
		prior(1) = 0.8;
		prior(2) = 0.2;

		% Cluster centres
		c = [0.0, 0.0; 2.0, 3.0];

		% Cluster standard deviations
		covar = [1 .3; .3 1];

		ind{1} = 1:prior(1)*ndata;
		for ic = 2 : ncentres
			ip = 1:prior(ic)*ndata;
			ind{ic} = ip+max(ind{ic-1});
		end% for ic
		for ic = 1 : ncentres
			ip = ind{ic};
			
		end% for ic
		

		
	case 3
		ncentres = 2;

		% Priors for the clusters
		prior(1) = 0.5;
		prior(2) = 0.5;

		% Cluster centres
		c = [0.0, 0.0; 2.0, 2.0];

		% Cluster standard deviations
		sd  = [1 1];

		% Put first cluster 

		x = data(1:prior(1)*ndata, 1) * sd(1) + c(1,1);
		y = data(1:prior(1)*ndata, 2) * sd(1) * 2 + c(1,2);
		[th,r] = cart2pol(x,y);
		[x,y]  = pol2cart(th-pi/4,r);
		data(1:prior(1)*ndata, 1) = x;
		data(1:prior(1)*ndata, 2) = y;

		% Leave second cluster at (0,0)
		x = data((prior(1)*ndata + 1):(prior(2)+prior(1))*ndata, 1) * sd(2) + c(2,1);
		y = data((prior(1)*ndata + 1):(prior(2)+prior(1))*ndata, 2) * sd(2) + c(2,2);
		data((prior(1)*ndata + 1):(prior(2)+prior(1))*ndata, 1) = x;	
		data((prior(1)*ndata + 1):(prior(2)+prior(1))*ndata, 2) = y;

	case 4
		ncentres = 2;

		% Priors for the three clusters
		prior(1) = 0.8;
		prior(2) = 0.2;

		% Cluster centres
		c = [0.0, 0.0; 2.0, 3.0];


		% Put first cluster 
		sd  = [0.5 0.5];
		x = data(1:prior(1)*ndata, 1) * sd(1) + c(1,1);
		y = data(1:prior(1)*ndata, 2) * sd(2) * 2 + c(1,2);
		[th,r] = cart2pol(x,y);
		[x,y]  = pol2cart(th-pi/3,r);
		data(1:prior(1)*ndata, 1) = x;
		data(1:prior(1)*ndata, 2) = y;

		% Leave second cluster at (0,0)
		sd  = [.3 .7];
		x = data((prior(1)*ndata + 1):(prior(2)+prior(1))*ndata, 1) * sd(1) + c(2,1);
		y = data((prior(1)*ndata + 1):(prior(2)+prior(1))*ndata, 2) * sd(2) + c(2,2);
		data((prior(1)*ndata + 1):(prior(2)+prior(1))*ndata, 1) = x;	
		data((prior(1)*ndata + 1):(prior(2)+prior(1))*ndata, 2) = y;


	case 5
		ncentres = 4;

		% Priors for the three clusters
		prior(1) = 0.1;
		prior(2) = 0.3;
		prior(3) = 0.5;
		prior(4) = 0.1;

		% Cluster centres
		c = [0.0, 0.0; 2.0, 3.0; 1.0 0.0; -1.0 -2.0];
%		c = [0.0, 0.0; 2.0, 3.0; 1.0 0.0;  1.0 .0];

		% Cluster standard deviations
		sd  = [.4 .5 .8 .1];

		% Put first cluster at (2, 3.5)
		ind{1} = 1:prior(1)*ndata;
		for ic = 2 : ncentres
			ip = 1:prior(ic)*ndata;
			ind{ic} = ip+max(ind{ic-1});
		end% for ic
		
		for ic = 1 : ncentres
			data(ind{ic},1) = data(ind{ic},1)*sd(ic) + c(ic,1);
			data(ind{ic},2) = data(ind{ic},2)*sd(ic) + c(ic,2);
		end% for ic
		
	case 6
		ncentres = 4;

		% Priors for the three clusters
		prior(1) = 0.1;
		prior(2) = 0.3;
		prior(3) = 0.5;
		prior(4) = 0.1;

		% Cluster centres
		c = [0.0, 0.0; 2.0, 3.0; 1.0 0.0; -1.0 -2.0];
		c = [0.0, 0.0; 2.0, 3.0; 1.0 0.0;  1.0 .0];

		% Cluster standard deviations
		sd  = [.3 .5 .8 .1];

		% Put first cluster at (2, 3.5)
		ind{1} = 1:prior(1)*ndata;
		for ic = 2 : ncentres
			ip = 1:prior(ic)*ndata;
			ind{ic} = ip+max(ind{ic-1});
		end% for ic

		for ic = 1 : ncentres
			data(ind{ic},1) = data(ind{ic},1)*sd(ic) + c(ic,1);
			data(ind{ic},2) = data(ind{ic},2)*sd(ic) + c(ic,2);
		end% for ic

end% switch 
