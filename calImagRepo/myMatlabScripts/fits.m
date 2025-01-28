%%
%close all;
    figure;
    dataTypes = {'raw','shuffle'};
for k = 1:length(dataTypes)
    dataType = dataTypes{k};
    animalData = data.(['m2024']).decCells;
    conf = 1;
    %%
    x = animalData(conf).vec(sessIdx).cellVect;
    switch dataType
        case 'raw'
            y = mean(mean(animalData(conf).vec(sessIdx).rawData,3));
        case 'shuffle'
            y = mean(mean(animalData(conf).vec(sessIdx).shuffledData,3));
    end
    x = x(12:end);
    y = y(12:end);
    %%
    %Exponential fit
    expFn = @(b,x) b(4).^(-b(1)*x + b(2)) + b(3);
    b0Exp = [1/90,2,1,1];

    wExp = ones([1,length(x)]);%[1:length(x)].^2;
    [betaExp,~,~,~,mseExp] = nlinfit(x,y,expFn,b0Exp','Weights',wExp);
    mseExp = mseExp./std(y);
    %%
    %Lineal fit
    xLog = x%log10(x);
    yLog = log10(y);

    linFn = @(b,x)b(1)*x + b(2);
    b0Lin = [1,2];

    wLin = [10^-14,ones([1,length(x)-1])];%[1:length(x)].^2;
    [betaLin,~,~,~,mseLin] = nlinfit(xLog,yLog,linFn,b0Lin','Weights',wLin);
    mseLin = mseLin./std(yLog);
    %%
    %Polynomial fit 2 order
    pol2Fn = @(b,x) b(1)*(x.^-b(2)) + b(3);
    b0Pol2 = [1,1,0];

    [betaPol2,~,~,~,msePower] = nlinfit(x,y,pol2Fn,b0Pol2');
    msePower = msePower./std(y);
    %%
    nPoints = 1000;


    annotation('textbox',[0,0,1,k/2-.4],'String',[dataTypes{k}],'FitBoxToText','off','FontSize',20,'FontWeight','bold','LineStyle','none');

    
    subplot_er(2,3,1+(3*(k-1)))
    plot(x,y,'.-','LineWidth',2,'MarkerSize',12);
    hold on;
    xTest = linspace(x(1),x(end)*2,nPoints);
    plot(xTest,expFn(betaExp,xTest),'LineWidth',2)
    legend('Original','Fit')
    grid on;
    title(['Exponential fit, MSE/var(y): ',num2str(mseExp)]);
    xlabel('x')
    ylabel('y')
    axis square

    subplot_er(2,3,2+(3*(k-1)))
    plot(xLog,yLog,'.-','LineWidth',2,'MarkerSize',12);
    hold on;
    xTest = linspace(xLog(1),xLog(end)*2,nPoints);
    plot(xTest,linFn(betaLin,xTest),'LineWidth',2)
    grid on;
    title(['Linear, MSE/var(log(y)): ', num2str(mseLin)]);
    xlabel('Log(x)')
    ylabel('Log(y)')
    axis square

    subplot_er(2,3,3+(3*(k-1)))
    plot(x,y,'.-','LineWidth',2,'MarkerSize',12);
    hold on;
    xTest = linspace(x(1),x(end)*2,nPoints);
    plot(xTest,pol2Fn(betaPol2,xTest),'LineWidth',2)
    grid on;
    title(['Power Law, MSE/var(y): ', num2str(msePower)])
    xlabel('x')
    ylabel('y')
    axis square
end
