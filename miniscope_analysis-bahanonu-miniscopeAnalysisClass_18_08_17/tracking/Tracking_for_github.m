%% code for tracking a black mouse on a white background
%INPUTS: movie_filename: the absolute or relative path of the movie
%tracking_movie_filename: the absolute or relative path of the
%tracking move
%write_tracking_movie: 1 or 0 on whether or not to save a separate
%tracking moving

%% VARIABLE PARAMETERS
% num_pixels: this is the number of pixels to use for the tracking

%OUTPUTS: (implicit) Saves output files of a tracked movie
% centroid: 2XNFrames matrix of animal centroids
%orientation: 1xN frames of animal orientaiton (doesn't really work)

%% Example
% [centroidUser,orientation] = Tracking_for_github('path/to/file.avi',...
% 'path/to/file_tracking.avi',0,'parallel',1,'num_pixels',1000);
% or to only analyze a specific subset of frames
% [centroidUser,orientation] = Tracking_for_github('path/to/file.avi',...
% 'path/to/file_tracking.avi',0,'parallel',1,'num_pixels',1000,'nFramesAnalyze',1000);

%% JESSE MARSHALL, 10/14/2015
% algorithm is well tested, but this implimentation has not been

% changelog
    % 2015.10.30 [biafra] - adding overwrite_tracking_movie, changed
    % options input to standard format for repo, added option for parallel
    % processing

function [centroid,orientation] = Tracking_for_github(movie_filename,tracking_movie_filename,write_tracking_movie,varargin)

    %% options
    % overwrite tracking movie if it already exists, 0 = no, 1 = yes
    options.overwrite_tracking_movie = 1;
    % run a parallel version of the script, 0 = no, 1 = yes
    options.parallel = 1;
    % this is the number of pixels to use for the tracking
    options.num_pixels = 1500;
    % number of frames in the movie to analyze, blank = all frames
    options.nFramesAnalyze = [];
    % do the image complement of the movie for white mice on dark background
    options.imcomplementMovie = 0;

    % get options
    options = getOptions(options,varargin);

    % display(options)
    % % unpack options into current workspace
    % fn=fieldnames(options);
    % for i=1:length(fn)
    %     eval([fn{i} '=options.' fn{i} ';']);
    % end

    overwrite_tracking_movie = options.overwrite_tracking_movie;

    num_pixels = options.num_pixels;

    %% movie information
    movie_obj = VideoReader(movie_filename);

    if isempty(options.nFramesAnalyze)
        nFrames = movie_obj.NumberOfFrames;
    else
        nFrames = options.nFramesAnalyze;
    end
    vidHeight = movie_obj.Height;
    vidWidth = movie_obj.Width;

    %% do centroid, thresholding and line plot, total distance
    centroid = zeros(2,nFrames);
    orientation = zeros(1,nFrames);
    %region_index = cell(1,nFrames);


    %% get a mask of the important area
    title('draw the mask')
    %%temp, until real reference image
    refimage = read(movie_obj, 20);
    h = figure(10)
    plot(1);hold on
    imagesc(refimage);
    figure(10)
    title('draw the mask')
    drawnow
    BW = roipoly(refimage);
    title('draw an area around the darkest portions of the mouse')
    BW_mouse = roipoly(refimage);
    BW_mouse = BW_mouse~=0;
    num_pixels = sum(BW_mouse(:));
    num_pixels
    title('draw length of container')
    [x,y] = getline(h);
    lineval = sqrt( (x(2)-x(1)).^2 + (y(2)-y(1)).^2);
    fprintf('chamber width: \t %f \n',lineval);
    BW2 = [];
    same_flag = 1;
    hold off



    %% save a copy showing the tracking quality
    if ( write_tracking_movie)
        if(~exist(tracking_movie_filename,'file') || (overwrite_tracking_movie == 1))

            writerObj2 = VideoWriter(tracking_movie_filename,'Motion JPEG AVI');
            writerObj2.Quality = 99;
            open(writerObj2);
        end
    else
        writerObj2 = [];
    end

    %% check maximum number of cores available
    if options.parallel==1
        display('opening worker')
        maxCores = feature('numCores')*2-1;
        if maxCores>7
            maxCores = 7;
        end
        % check that local matlabpool configuration is correct
        myCluster = parcluster('local');
        if myCluster.NumWorkers<maxCores
            myCluster.NumWorkers = maxCores; % 'Modified' property now TRUE
            saveProfile(myCluster);   % 'local' profile now updated
        end
        % open works = max core #, probably should do maxCores-1 for stability...
        % check whether matlabpool is already open
        if matlabpool('size') | ~options.parallel
        else
            matlabpool('open',maxCores);
        end
    end

    %% do the tracking
    tic
    parfor k = 1 : nFrames
        if (mod(k,100) == 0)
            fprintf('frame %d \n',k)
        end

        movie_obj = VideoReader(movie_filename);

        nFrames = movie_obj.NumberOfFrames;
        vidHeight = movie_obj.Height;
        vidWidth = movie_obj.Width;

        frame = (read(movie_obj, k));
        if (size(frame,3))
            frame = mean(frame,3);
        end
        movie_frame =  frame;
        movie_frame_full = movie_frame;

        movie_frame_full(BW == 0) = 255;

        %% A Separate way of doing movement tracking, by identifying connected components
        %                         movie_frame = uint8(BW_map{day_ind,mouse_number}).*(mean_movie- movie_frame);
        %
        %                         if (k==1)
        %                             threshold = max(max(max(movie_frame)))./3;
        %                         end
        %
        %                                                 %% old way of computing the centroid using bwlabel and a threshold
        %                         movie_frame(movie_frame < threshold) = 0;
        %                         comp = bwconncomp(movie_frame);
        %                         temp = zeros(1,comp.NumObjects);
        %
        %                         stats = regionprops(comp,movie_frame,...
        %                             'Area','Centroid','MeanIntensity','Orientation', 'MajorAxisLength',...
        %                             'Eccentricity','MinorAxisLength');
        %
        %                         for kk=1:comp.NumObjects
        %                             temp(kk) = stats(kk).MeanIntensity*stats(kk).Area;
        %                         end
        %                         [~,mouse_index] = max(temp);
        %
        %                         %if there are any objects, depending on the frame
        %                         %number
        %                         if (comp.NumObjects)
        %                             orientation(k) = stats(mouse_index).Orientation;
        %                             centroid(:,k) = stats(mouse_index).Centroid(1:2);
        %                         elseif (~comp.NumObjects && k>1)
        %                             orientation(k) = orientation(k-1);
        %                             centroid(:,k) = centroid(k-1);
        %                         else
        %                               orientation(k) = 100;
        %                             centroid(:,k) = 100;
        %                         end

        %% now just take center of bottom pixels
        %   movie_frame = uint8(BW).*(255-(uint8(movie_frame)-mean_movie));
        [temp,top_pixels] = sort(reshape(255-movie_frame_full,1,[]),'descend');
        top_pixels = top_pixels(1:num_pixels);
        [ pixel_x,pixel_y] = ind2sub([size(movie_frame_full,1),size(movie_frame_full,2)],top_pixels);
        cent_x = mean(pixel_x);
        cent_y = mean(pixel_y);

        centroid(:,k) = [cent_x,cent_y];
        points = [pixel_x' pixel_y'];
        points = bsxfun(@minus, points,mean(points,1));
        [~,~,V] = svd(points);
        orientation(k) = atan(V(1,1)./V(2,1))*180/pi;

        frame_temp = 255-movie_frame_full;

        %% draw a line along
        arm_length = 25;
        line_pts = arm_length;
        line_start =  [cent_x,cent_y];
        line_end = [(cent_x+arm_length.*sind(orientation(k))) (cent_y+arm_length.*cosd(orientation(k)))];
        line_x = line_start(1):(line_end(1)-line_start(1))./line_pts:line_end(1);
        line_y = line_start(2):(line_end(2)-line_start(2))./line_pts:line_end(2);

        if(~numel(line_x))
            line_x = ones(1,arm_length);
        elseif(~numel(line_y))
            line_y = ones(1,arm_length);
        end
        line_y = line_y(1:numel(line_x));

        line_x(line_x > (size(frame_temp,1)-1)) = size(frame_temp,1)-1;
        line_y(line_y > (size(frame_temp,2)-1)) = size(frame_temp,2)-1;
        line_x(line_x < 1) = 1;
        line_y(line_y < 1) = 1;

        line_inds = floor(cat(1,line_x,line_y));


        %      if (k>1)
        %      test = zeros(vidHeight,vidWidth);
        %      test(region_index{k}
        %      end
        %
        %  frame_temp(top_pixels) = 0;
        frame_temp(sub2ind(size(frame_temp),line_inds(1,:),line_inds(2,:))) = 0;
        frame_temp(sub2ind(size(frame_temp),line_inds(1,:),line_inds(2,:)+1)) = 0;
        frame_temp = uint8(frame_temp);

        %% visualize the frames
        % frame_temp = zeros(size(movie_frame_full));
        %  frame_temp(top_pixels) = 1;

        %imagesc(frame_temp)
        last_movie_frame = movie_frame;

        if write_tracking_movie==1&options.parallel==0
            if(~exist(tracking_movie_filename,'file') || (overwrite_tracking_movie == 1))
                writeVideo(writerObj2,squeeze(frame_temp));
            end
        end
    end
    toc
    %% overwrite movie
    if ( write_tracking_movie)
        if(~exist(tracking_movie_filename,'file') || (overwrite_tracking_movie == 1))
            close(writerObj2);
        end
    end
end