function plot2Dpca(X)

if nargin < 1
    X = mvnrnd([0,0],[10,9.9;9.9,10],1000);
    figure;
end

coeff = pca(X);
meanM = mean(X);
T = (X-meanM)*coeff + meanM;
Tl = (X-meanM)*coeff(:,1);
Xrec = Tl*coeff(:,1)'  + meanM;
x1 = linspace(min(X(:,1)),max(X(:,1)),1000);
x2 = linspace(min(X(:,2)),max(X(:,2)),1000);
x1L = coeff(:,1)*(x1 - mean(x1)) + meanM';
x2L = coeff(:,2)*(x2 - mean(x2)) + meanM';
plot(X(:,1),X(:,2),'x');hold on;
plot(T(:,1),T(:,2),'x')
plot(Xrec(:,1),Xrec(:,2),'x');
plot(x1L(1,:),x1L(2,:),'k');
plot(x2L(1,:),x2L(2,:),'--k');
legend('Original','PCspace','Reconstructed','PC1','PC2')
axis equal;
end

