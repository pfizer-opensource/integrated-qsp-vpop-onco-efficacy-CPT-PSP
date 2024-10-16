function [pd,mu,sigma,r] = get_distribution_from_data(data,data_to_fit,ln_transform,gm_comp,do_plot)

% This function fits a Gaussian mixture model to the clinical endpoints
% that will be used for scoring the plausible and virtual population. This
% function includes the option to log inputs, allowing for a mixture of
% normal and log-normal distributed data. This function also includes code
% to optionally plot the Gaussian mixture fit alongside the clinical data
% to assess the fit quality.

%Input
% data -- clinical data to be fit 
% data_to_fit -- columns names of outputs in the clinical data to fit
% ln_transform -- boolean indicating which of above columns to log transform
% gm_comp -- number of Guassian components to use
% do_plot -- boolean indicating if goodness of fit plot should be generated

%Output
% pd -- probability distribution function for Guassian fit
% mu -- means of Guassian mixture
% sigma -- covariances of gaussian mixture
% r -- data columns used for fitting, w/ log transformed if indicated

% Gather the variables of interest
len_vars = length(data_to_fit);
r = data{:,data_to_fit};
ln_ind = find(ln_transform);
r(:,ln_ind) = log(r(:,ln_ind));

% r(:,1)=log(r(:,1));
% r(:,2)=log((r(:,2)/100+1)+0.1);
% r(:,3)=log(r(:,3)./(501-r(:,3)));

r = r(~any(isnan(r'),1)',:); % Clear any rows with NaN, not absolutely necessary, but avoids a warning message
c = corrcoef(r);
% Fit a multivariate Gaussian model to the data:
pd = fitgmdist(r,gm_comp,'Replicates',10);
mu = pd.mu;
sigma = pd.Sigma;
alpha = pd.ComponentProportion;

if(do_plot)
    total_plots = len_vars*(len_vars-1)/2;
    ncols = ceil(total_plots/2);
    [x,y] = meshgrid(1:len_vars, 1:len_vars);
    m = [x(:), y(:)]; n = m(m(:,2) < m(:,1),:);
    figure(111)
    %     figure(1)
    for i = 1:total_plots
        subplot(2,ncols,i)
        scatter(r(:,n(i,1)),r(:,n(i,2)));
        hold on;
        [prob, X1Grid, X2Grid] = get_marginal(r,mu,sigma,alpha,n(i,1),n(i,2));
        contour(X1Grid, X2Grid, prob);
%         surf([X1Grid, X1Grid], [X2Grid, X2Grid], [prob,prob]);
        box on;
        xlabel(data_to_fit{n(i,1)}); ylabel(data_to_fit{n(i,2)});
        drawnow;
    end    
    
    %%%% 1-D histograms
    figure(112)
    ncols = ceil(len_vars/2);
    for i = 1:len_vars
    subplot(2,ncols,i)
    h = histogram(r(:,i), 8, 'Normalization','pdf');
    hold on;
    [prob, X1Grid, ~] = get_marginal(r,mu,sigma,alpha,i,i);
    plot(X1Grid,prob);
    box on;
   xlabel(data_to_fit{i});
    drawnow;
    end
end

end
function [prob_marginal, X1Grid, X2Grid] = get_marginal(r,mu,sigma,alpha,i,j)
if (i ~= j)
    [X1Grid, X2Grid] = meshgrid(sortrows(r(:,i)), sortrows(r(:,j)));
    pd = gmdistribution(mu(:,[i,j]),sigma([i,j],[i,j],:),alpha);
    prob_marginal = reshape(arrayfun(@(x,y) pdf(pd, [x,y]), X1Grid(:), X2Grid(:)), size(X1Grid,1), size(X1Grid,2));
else
    X1Grid = sortrows(r(:,i));
    pd = gmdistribution(mu(:,i),sigma(i,i,:),alpha);
    prob_marginal = arrayfun(@(x) pdf(pd, x), X1Grid(:));
    X2Grid = [];
end
end
