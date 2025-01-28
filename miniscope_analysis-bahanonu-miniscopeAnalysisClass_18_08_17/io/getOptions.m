function [options] = getOptions(options,inputArgs,varargin)
    % gets default options for a function, replaces with inputArgs inputs if they are present
    % biafra ahanonu
    % started: 2013.11.04
    %
    % inputs
    %   options - structure with options
    %   inputArgs - an even numbered cell array, with {'option','value'} as the ordering. Normally pass varargin.
    %
    % NOTE
    %   use the 'options' name-value pair to input an options structure that will overwrite default options in a function, example below.
    %   options.Stargazer = 1;
    %   options.SHH = 0;
    %   getMutations(mutationList,'options',options);
    %
    %   This is in contrast to using name-value pairs, both will produce the same result.
    %   getMutations(mutationList,'Stargazer',1,'SHH',0);
    %
    % USAGE
    %   %========================
    %   options.movieType = 'tiff';
    %   % get options
    %   options = getOptions(options,varargin);
    %   % unpack options into current workspace, comment out if want to just call options structure
    %   fn=fieldnames(options);
    %   for i=1:length(fn)
    %       eval([fn{i} '=options.' fn{i} ';']);
    %   end
    %   %========================

    % changelog
    % 2014.02.12 [11:56:00] - added feature to allow input of an options structure that contains the options instead of having to input multiple name-value pairs. - Biafra
    % 2014.07.10 [05:19:00] - added displayed warning if an option is input that was not present (this usually indicates typo) - Lacey (merged)
    % 2014.12.10 [19:32:54] - now gets calling function and uses that to get default options - Biafra
    % 2015.08.24 [23:31:36] - updated comments. - Biafra
    % 2015.12.03 [13:52:15] - Added recursive aspect to mirrorRightStruct and added support for handling struct name-value inputs. mirrorRightStruct checks that struct options input by the user are struct in the input options. - Biafra

    % TODO
    % allow input of an option structure - DONE!
    % call settings function to have defaults for all functions in a single place - DONE!
    % allow recursive overwriting of options structure - DONE!
    % Type checking of all field names input by the user?

    %========================
    % whether getOptions should use recursive structures
    goptions.recursiveStructs = 1;
    % whether to show warning if option input that is not in original structure
    goptions.showWarnings = 1;
    % filter through options
    for i = 1:2:length(varargin)
        inputField = varargin{i};
        if isfield(goptions, inputField)
            inputValue = varargin{i+1};
            goptions.(inputField) = inputValue;
        end
    end
    % don't do this! recursion with no base case waiting to happen...
    % goptions = getOptions(goptions,varargin);
    %========================

    % Get default options for a function
    [ST,I] = dbstack;
    % fieldnames(ST)
    parentFunctionName = {ST.name};
    parentFunctionName = parentFunctionName{2};
    [optionsTmp] = getSettings(parentFunctionName);
    if isempty(optionsTmp)
        % Do nothing, don't use defaults if not present
    else
        options = optionsTmp;
    end

    % Get list of available options
    validOptions = fieldnames(options);

    % Loop over all input arguments, overwrite default/input options
    for i = 1:2:length(inputArgs)
        % inputArgs = inputArgs{1};
        val = inputArgs{i};
        if ischar(val)
            %display([inputArgs{i} ': ' num2str(inputArgs{i+1})]);
            if strcmp('options',val)
                % Special options struct, only add field names defined by the user. Keep all original field names that are not input by the user.
                inputOptions = inputArgs{i+1};
                options = mirrorRightStruct(inputOptions,options,goptions,val);
            elseif sum(strcmp(val,validOptions))>0&isstruct(options.(val))&goptions.recursiveStructs==1
                % If struct name-value, add users field name changes only, keep all original field names in the struct intact, struct-recursion ON
                inputOptions = inputArgs{i+1};
                options.(val) = mirrorRightStruct(inputOptions,options.(val),goptions,val);
            elseif sum(strcmp(val,validOptions))>0
                % Non-options, non-struct value, struct-recursion OFF
                % elseif ~isempty(strcmp(val,validOptions))
                % Way more elegant, directly overwrite option
                options.(val) = inputArgs{i+1};
                % eval(['options.' val '=' num2str(inputArgs{i+1}) ';']);
            else
                if goptions.showWarnings==1
                    try
                        [ST,~] = dbstack;
                        callingFxn = ST(2).name;
                        callingFxnPath = which(ST(2).file);
                        callingFxnLine = num2str(ST(2).line);
%                         disp(['<strong>WARNING</strong>: <a href="">' val '</a> is not a valid option for <a href="matlab: opentoline(' callingFxnPath ',' callingFxnLine ')">' callingFxn '</a> on line ' callingFxnLine])
                    catch err
                        callingFxn = 'UNKNOWN FUNCTION';
%                         disp(['<strong>WARNING</strong>: <a href="">' val '</a> is not a valid option for "' callingFxn '"'])
                    end
                end
            end
        else
            if goptions.showWarnings==1
                try
                    [ST,~] = dbstack;
                    callingFxn = ST(2).name;
                    callingFxnPath = which(ST(2).file);
                    callingFxnLine = num2str(ST(2).line);
                    disp(['<strong>WARNING</strong>: enter the parameter name before its associated value in <a href="matlab: opentoline(' callingFxnPath ',' callingFxnLine ')">' callingFxn '</a> on line ' callingFxnLine])
                catch err
                    callingFxn = 'UNKNOWN FUNCTION';
                    disp(['<strong>WARNING</strong>: enter the parameter name before its associated value in "' callingFxn '"'])
                end
            end
            continue;
        end
    end
    %display(options);
end

function [toStruct] = mirrorRightStruct(fromStruct,toStruct,goptions,toStructName)
    % Overwrites fields in toStruct with those in fromStruct, other toStruct fields remain intact.
    % More generally, copies fields in fromStruct into toStruct, if there is an overlap in field names, fromStruct overwrites.
    % Fields present in toStruct but not fromStruct are kept in toStruct output.
    fromNames = fieldnames(fromStruct);
    for name = 1:length(fromNames)
        fromField = fromNames{name};
        % if a field name is a struct, recursively grab user options from it
        if isfield(toStruct, fromField)
            if isstruct(fromStruct.(fromField))&goptions.recursiveStructs==1
                % safety check: field exist in toStruct and is also a structure
                if isfield(toStruct, fromField)&isstruct(toStruct.(fromField))
                    toStruct.(fromField) = mirrorRightStruct(fromStruct.(fromField),toStruct.(fromField),goptions,[toStructName '.' fromField]);
                else
                    [ST,~] = dbstack;
                    callingFxn = ST(3).name;
                    callingFxnPath=which(ST(3).file);
                    callingFxnLine = num2str(ST(3).line);
                    disp(['<strong>WARNING</strong>: <a href="">' toStructName '.' fromField '</a> is not a originally a STRUCT, ignoring...'])
                end
            else
                toStruct.(fromField) = fromStruct.(fromField);
            end
        else
            if goptions.showWarnings==1
                try
                    [ST,~] = dbstack;
                    callingFxn = ST(3).name;
                    callingFxnPath=which(ST(3).file);
                    callingFxnLine = num2str(ST(3).line);
%                     disp(['<strong>WARNING</strong>: <a href="">' toStructName '.' fromField '</a> is not a valid option for <a href="matlab: opentoline(' callingFxnPath ',' callingFxnLine ')">' callingFxn '</a> on line ' callingFxnLine])
                catch err
                    callingFxn = 'UNKNOWN FUNCTION';
%                     disp(['<strong>WARNING</strong>: <a href="">' toStructName '.' fromField '</a> is not a valid option for "' callingFxn '"'])
                end
            end
        end
    end
end