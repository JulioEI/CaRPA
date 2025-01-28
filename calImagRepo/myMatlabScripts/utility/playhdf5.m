function playhdf5(hdf5)
if isstr(hdf5)
    hinf = hdf5info(hdf5);
    dataset = hdf5read(hinf.GroupHierarchy.Datasets);
else
    dataset = hdf5;
end
playMovie(dataset)
end

