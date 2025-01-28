function [realTraces,realParams,imgs,eventTimes] = getSimData(simDataStore,ind1,ind2,ind3,imgsStorePath)

realTraces = simDataStore.traces{ind1,ind2,ind3};
realParams = simDataStore.params{ind1,ind2,ind3};
eventTimes = simDataStore.eventTimes{ind1,ind2,ind3};

imgs = load([imgsStorePath filesep 'imgs_' num2str(ind1) '_' num2str(ind2) '_' num2str(ind3)],'imgs');
imgs=imgs.imgs;
