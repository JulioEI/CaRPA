function [dataPerAnimalConcat] = concatDataPerAnimal(dataPerAnimal,varargin)
    doMean = false;
    if ~isempty(varargin)
        if strcmp(varargin{1},'doMean')
            doMean = varargin{2};
        end
    end
    animalFields = fields(dataPerAnimal);
    dataPerAnimalConcat = [];
    clear dayFeat
    for animal = animalFields'
        dataPerDay = dataPerAnimal.(animal{1});
        dayFields = fields(dataPerDay);
        for day = dayFields'
            dataPerFeat = dataPerDay.(day{1});
            featFields = fields(dataPerFeat);
            for feat = featFields'
                cellData = dataPerFeat.(feat{1});
                if iscell(cellData)
                    %concatEvent = zeros([sum(cellfun(@length,cellData)),1]);
                    concatEvent = [];
                    for eventData = cellData'
                        if doMean
                            concatEvent = [concatEvent;nanmean(eventData{1})'];
                        else
                            concatEvent = [concatEvent;eventData{1}'];
                        end
                    end
                    cellData = concatEvent;
                end
                try
                    dayFeat.(feat{1}) = [dayFeat.(feat{1});cellData];
                catch
                    dayFeat.(feat{1}) = cellData;
                end
            end
        end
        dataPerAnimalConcat.(animal{1}) = dayFeat;
        clear dayFeat
    end
end

