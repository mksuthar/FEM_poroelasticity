function [Melements] = extract_mesh2_elem(filename,lineBegin, lineEnd)
	% Read elements data, assume we are reading data x amount of column,
	% if data is not found, replace it with NaN object. Since we can read 
	% 2nd colum, we can figure out type of element, and how many points node
	% it has to extract it later. 
    fid = fopen(filename,'r');
    formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f';
    Melements = textscan(fid, formatSpec, lineEnd-lineBegin-2, 'Delimiter', ' ', 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN,'HeaderLines', lineBegin+1, 'ReturnOnError', false);
end