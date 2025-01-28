function saveSimResults(simDataStore,simStorePath) %#ok<INUSL>

if exist([simStorePath filesep 'simResults' '.mat'], 'file')
    saveInd=1;
    while exist([simStorePath filesep 'simResults' '_' num2str(saveInd) '.mat'], 'file')
        saveInd=saveInd+1;
    end
    save([simStorePath filesep 'simResults' '_' num2str(saveInd) '.mat'],'simDataStore', '-v7.3')
else
    save([simStorePath filesep 'simResults'],'simDataStore','-v7.3')
end