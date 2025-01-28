function [y] = TBcostFn(x,allCellMaps,allCellMapsMax)

    ioptions.RegisType = 2;
    ioptions.meanSubtract = 0;
    ioptions.complementMatrix = 0;
    ioptions.normalizeType = 'divideByLowpass';
    ioptions.RegisType = 2;
    ioptions.SmoothX = x(1);
    ioptions.SmoothY = x(1);
    ioptions.minGain = 0;
    ioptions.Levels = 6;
    ioptions.Lastlevels = 1;
    ioptions.Epsilon = 1.1921e-07;
    ioptions.zapMean = 0;
    ioptions.parallel = 1;
    ioptions.cropCoords = [];
    ioptions.closeMatlabPool = 0;
    ioptions.removeEdges = 0;
    ioptions.registrationFxn = 'transfturboreg';
    ioptions.turboregRotation = 1;

    MprFullTB = cell([1,length(allCellMapsMax)]);
    iterations = x(2);
    templateIdx = round(length(allCellMapsMax)./2);
    for j = 1:length(allCellMapsMax)
        ioptions.altMovieRegister = allCellMaps{j};
        MprFullTB{j} = allCellMapsMax{j};
        for iter = 1:iterations
            [~,coords] = turboregMovie(cat(3,allCellMapsMax{templateIdx},MprFullTB{j}),'options',ioptions);
            ioptions.altMovieRegister = turboregMovie(allCellMaps{j},'precomputedRegistrationCooords',coords);
            MprFullTB{j} = turboregMovie(allCellMapsMax{j},'precomputedRegistrationCooords',coords);
        end
    end
    corrMatTB = zeros([1,length(allCellMapsMax)]);
    for j = 1:length(allCellMapsMax)
        corrMatTB(j) = corr2(allCellMapsMax{templateIdx},MprFullTB{j});
    end
    y = -mean(abs(corrMatTB(setdiff(1:length(allCellMapsMax),templateIdx))));
end

