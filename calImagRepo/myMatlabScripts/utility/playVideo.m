function [] = playVideo(path)
    figure;
    v = VideoReader(path);
    currAxes = axes;
    while hasFrame(v)
        vidFrame = readFrame(v);
        image(vidFrame, 'Parent', currAxes);
        currAxes.Visible = 'off';
        pause(1/v.FrameRate);
    end
end

