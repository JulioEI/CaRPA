function appendPosToTracesEventsFile(fullPath,behavior)
    load(fullPath);
    if ~iscell(behavior)
        behavior = {behavior};
    end
    tracesEvents.position = [];
    tracesEvents.velocity = [];
    for behaviorFile = behavior
        [position,velocity] = getMouseTrajectory(behaviorFile{1});
        tracesEvents.position = [tracesEvents.position,position];
        tracesEvents.velocity = [tracesEvents.velocity,velocity];
    end
    save(fullPath,'tracesEvents')
end

