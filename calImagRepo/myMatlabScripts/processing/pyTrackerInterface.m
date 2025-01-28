function [position,velocity,score] = pyTrackerInterface(movie)
% [position,velocity,score] = pyTrackerInterface(Mouse3014A_17_11_1113_20_14)
%Put the pyTracker file in the python path
trackerPath = fileparts(which('pyTrackerAdaptative.py'));
if count(py.sys.path,trackerPath) == 0
    insert(py.sys.path,int32(0),trackerPath);
end
py.importlib.import_module('pyTrackerAdaptative');

%Call the python tracker
position = py.pyTrackerAdaptative.tracker(movie)

%Check for read errors
if strcmp(char(position),'ERROR')
    %Fix video
    mov = VideoReader(movie);
    [pathstr,name,ext] = fileparts(movie);
    outputVideoName = [pathstr,filesep,name,'_fix',ext];
    
    v = VideoWriter(outputVideoName);
    open(v)
    for k = 1:mov.NumberOfFrames
        frame = mov.read(k);
        writeVideo(v,frame);
    end
    close(v)
    %Run tracker
    position = py.pyTrackerAdaptative.tracker(outputVideoName);
    %Delete fixed movie
    delete(outputVideoName)
end

%Convert the numpy array to a matlab matrix
positionList = cell(position.tolist);
positionCell = cellfun(@(x) cell(x), positionList, 'UniformOutput', 0);
positionMat = [cellfun(@(x) x{1}, positionCell)',cellfun(@(x) x{2}, positionCell)'];

%Postprocess traces
[position,velocity,score] = postProcessTracesTracker(positionMat);

%Show results
% plotTracesOverMovie(position,velocity,movie)



