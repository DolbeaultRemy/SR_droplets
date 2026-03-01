function cumulants = computeCumulants(x, maxOrder)
%% computeCumulants
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Compute cumulants up to specified order from data vector x.
%
% Inputs:
%   x        - 1D numeric vector (may contain NaNs)
%   maxOrder - maximum order of cumulants to compute (default: 6)
%
% Output:
%   cumulants - vector [kappa_1, ..., kappa_maxOrder]
%
% Notes:
%   Syntax- cumulants = computeCumulants(x, maxOrder).

    if nargin < 2
        maxOrder = 6;
    end
    
    x = x(:);
    x = x(~isnan(x)); % Remove NaNs
    
    if isempty(x)
        cumulants = NaN(1, maxOrder);
        return;
    end
    
    mu1 = mean(x, 'omitnan');
    x_centered = x - mu1;
    
    cumulants = zeros(1, maxOrder);
    cumulants(1) = mu1;
    
    mu = zeros(1, maxOrder);
    for k = 2:maxOrder
        mu(k) = mean(x_centered.^k, 'omitnan');
    end
    
    if maxOrder >= 2
        cumulants(2) = mu(2);
    end
    if maxOrder >= 3
        cumulants(3) = mu(3);
    end
    if maxOrder >= 4
        cumulants(4) = mu(4) - 3 * mu(2)^2;
    end
    if maxOrder >= 5
        cumulants(5) = mu(5) - 10 * mu(3) * mu(2);
    end
    if maxOrder >= 6
        cumulants(6) = mu(6) - 15 * mu(4) * mu(2) - 10 * mu(3)^2 + 30 * mu(2)^3;
    end

end
