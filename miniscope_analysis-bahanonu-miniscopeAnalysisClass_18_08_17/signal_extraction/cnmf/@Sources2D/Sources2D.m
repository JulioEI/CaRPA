classdef Sources2D < handle
    
    % This class is a wrapper of Constrained NMF for standard 2D data. 
    % Author: Pengcheng Zhou, zhoupc1988@gmail.com with modifications from
    % Eftychios Pnevmatikakis
    
    %% properties 
    properties
        A;          % spatial components of neurons 
        C;          % temporal components of neurons 
        b;          % spatial components of backgrounds
        f;          % temporal components of backgrounds
        S;          % spike counts 
        Coor;       % neuron contours 
        options;    % options for model fitting 
        P;          % some estimated parameters 
    end
    
    %% methods
    methods
        %% constructor and options setting
        function obj = Sources2D(varargin)
            obj.options = CNMFSetParms(); 
            obj.P = struct('p', 2); 
            if nargin>0
                obj.options = CNMFSetParms(obj.options, varargin{:}); 
            end
        end 
        
        %% update parameters
        function updateParams(obj, varargin)
                obj.options = CNMFSetParms(obj.options, varargin{:});          
        end
        
        %% data preprocessing
        function Y = preprocess(obj,Y,p)
            [obj.P,Y] = preprocess_data(Y,p,obj.options);
        end
        
        %% fast initialization
        function [center] = initComponents(obj, Y, K, tau)
            [obj.A, obj.C, obj.b, obj.f, center] = initialize_components(Y, K, tau, obj.options);
        end
        
        %% update spatial components
        function updateSpatial(obj, Y)
            [obj.A, obj.b, obj.C] = update_spatial_components(Y, ...
                obj.C, obj.f, obj.A, obj.P, obj.options);
        end
        
        %% update temporal components
        function Y_res = updateTemporal(obj, Y)
            [obj.C, obj.f, Y_res, obj.P, obj.S] = update_temporal_components(...
                Y, obj.A, obj.b, obj.C, obj.f, obj.P, obj.options);
        end
        
        %% update temporal components in parallel
        function Y_res = updateTemporalParallel(obj, Y)
            [obj.C, obj.f, Y_res, obj.P, obj.S] = update_temporal_components_parallel(...
                Y, obj.A, obj.b, obj.C, obj.f, obj.P, obj.options);
        end
                
        %% merge found components
        function [nr, merged_ROIs] = merge(obj, Y_res)
            [obj.A, obj.C, nr, merged_ROIs, obj.P, obj.S] = merge_components(...
                Y_res,obj.A, obj.b, obj.C, obj.f, obj.P,obj.S, obj.options);
        end
        
        %% compute the residual
        function [Y_res] = residual(obj, Yr)
            Y_res = Yr - obj.A*obj.C - obj.b*obj.f;
        end
        
        %% take the snapshot of current results
        function [A, C, b, f, P] = snapshot(obj)
            A = obj.A;
            C = obj.C;
            b = obj.b;
            f = obj.f;
            P = obj.P;
        end
        
        %% extract DF/F signal after performing NMF
        function [C_df, Df, S_df] = extractDF_F(obj, Y, i)
            if ~exist('i', 'var')
                i = size(obj.A, 2) + 1;
            end
            
            [C_df, Df, S_df] = extract_DF_F(Y, [obj.A, obj.b],...
                [obj.C; obj.f], obj.S, i);
        end
        
        %% order_ROIs 
        function [srt] = orderROIs(obj)
            [obj.A, obj.C, obj.S, obj.P, srt] = order_ROIs(obj.A, obj.C,...
                obj.S, obj.P); 
        end 
        
        %% view contours 
        function [json_file] = viewContours(obj, Cn, contour_threshold, display)
            if or(isempty(Cn), ~exist('Cn', 'var') )
                Cn = reshape(obj.P.sn, obj.options.d1, obj.options.d2); 
            end
            [obj.Coor, json_file] = plot_contours(obj.A, Cn, ...
                contour_threshold, display); 
        end 
        
        %% plot components 
        function plotComponents(obj, Y, Cn)
            if ~exist('Cn', 'var')
                Cn = []; 
            end
            view_components(Y, obj.A, obj.C, obj.b, obj.f, Cn, obj.options); 
        end 
        
        %% make movie 
        function makePatchVideo(obj, Y) 
            make_patch_video(obj.A, obj.C, obj.b, obj.f, Y, obj.Coor,...
                obj.options); 
        end 
    end
    
end