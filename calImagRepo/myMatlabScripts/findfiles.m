
function flist = findfiles(pattern,basedir)
    fullpatt= [basedir,filesep,pattern];
    d = cellstr(ls(basedir));
    d(1:2) = [];
    d(~cellfun(@isdir,strcat(basedir,filesep,d))) = [];
    f = ls(fullpatt);
    if isempty(f)
        flist = {};
    else
        flist = strcat(basedir,filesep,cellstr(f));
    end
    for k = 1:length(d)
        flist = [flist;findfiles(pattern,[basedir,filesep,d{k}])];
    end
end
