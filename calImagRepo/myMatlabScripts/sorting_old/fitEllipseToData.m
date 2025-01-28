function [r_ellipse,X0,Y0] = fitEllipseToData(X,confidence)
%Remove nans from data
    data = X(~isnan(X(:,1)) & ~isnan(X(:,2)),:);

    % Calculate the eigenvectors and eigenvalues
    [coeff,~,latent] = pca(data);

    % Calculate the angle between the x-axis and the largest eigenvector
    angle = atan2(coeff(2,1),coeff(1,1));

    % This angle is between -pi and pi.
    % Let's shift it such that the angle is between 0 and 2pi
    if(angle < 0)
        angle = angle + 2*pi;
    end

    % Get the coordinates of the data mean
    avg = mean(data);

    % Get the 95% confidence interval error ellipse
    chiVal = chi2inv(confidence,2);
    theta_grid = linspace(0,2*pi);
    phi = angle;
    X0=avg(1);
    Y0=avg(2);
    a=sqrt(chiVal*latent(1));
    b=sqrt(chiVal*latent(end));

    % the ellipse in x and y coordinates 
    ellipse_x_r  = a*cos( theta_grid );
    ellipse_y_r  = b*sin( theta_grid );

    %Define a rotation matrix
    R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];

    %let's rotate the ellipse to some angle phi
    r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;    
end

