function [x, S] = OPCA_TwoFacies(fala, randVector, xm, gamma)
%  Authors: H. X. Vo and L. J. Durlorfky
%% This function is not written for efficiency but rather than for easy to understand
a = fala*randVector + xm; % PCA solution
x = zeros(length(a), 1);
mu = zeros(length(a), 1);  % Lagrange multiplier
eta = zeros(length(a), 1); % Lagrange multiplier
LagrangianGrad = zeros(1, length(a));
S = zeros(length(a), length(randVector)); % sensitivity matrix: S = dx/dxi
for i=1:length(a)
    if gamma > 1 % this scenario is not used as we will get a concave combined objective function 
        f0 = (1-gamma) * 0 - 2 * (a(i) - gamma/2) * 0;
        f1 = (1-gamma) * 1 - 2 * (a(i) - gamma/2) * 1;
        if f1 < f0
            x(i) = 1;
            mu(i) = 0;
            eta(i) = 2 * a(i) + gamma - 2;
        elseif f0 <= f1
            x(i) = 0;
            mu(i) = gamma - 2 * a(i);
            eta(i) = 0;
        end        
    elseif gamma == 1 % this scenario is not used as we will get a linear combined objective function 
        f0 = -2 * (a(i) - gamma/2) * 0;
        f1 = -2 * (a(i) - gamma/2) * 1;
        if f1 < f0
            x(i) = 1;
            mu(i) = 0;
            eta(i) = 2 * a(i) + gamma - 2;
        elseif f0 <= f1
            x(i) = 0;
            mu(i) = gamma - 2 * a(i);
            eta(i) = 0;
        end
    else % we use this scenario because we want a convex combined objective function
        xS = (a(i) - gamma/2) / (1 - gamma);
        if (xS < 1) && (xS > 0)
            x(i) = xS;
            mu(i) = 0;
            eta(i) = 0;
        elseif (xS >= 1)
            x(i) = 1;
            mu(i) = 0;
            eta(i) = 2 * a(i) + gamma - 2;
        else
            x(i) = 0;
            mu(i) = gamma - 2 * a(i);
            eta(i) = 0;
        end
        for j=1:length(randVector)
%             da_idxi_j = Sig(j)*U(i, j);
            da_idxi_j = fala(i, j);
            S(i, j) =  (2*da_idxi_j * x(i) * (1 - x(i))) /...
                (mu(i) + 2 * x(i) + eta(i) * x(i) - 2 * gamma * x(i) - mu(i) * x(i) + 2 * gamma * x(i)^2 - 2 * x(i)^2); % sensitivity matrix dx/dxi
        end
    end
    LagrangianGrad(i) = 2 * (1 - gamma) * x(i) - 2 * (a(i) - gamma/2) - mu(i) + eta(i);
end

end