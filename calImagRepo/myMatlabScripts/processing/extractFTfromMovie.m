function [ FT ] = extractFTfromMovie(movie,filters)


if isstr(movie)
hinf = hdf5info(movie);
movie = hdf5read(hinf.GroupHierarchy.Datasets);
end

if ~isa(movie, 'double')
    movie = double(movie);
end

movieR = movie;
movieR(isnan(movie)) = 0;

FT = reshape(filters,[size(filters,1)*size(filters,2),size(filters,3)])'*reshape(movieR,[size(movieR,1)*size(movieR,2),size(movieR,3)]);%/(size(movie,1)*size(movie,2));

end

