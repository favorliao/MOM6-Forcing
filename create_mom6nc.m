function create_mom6nc(filename,varname,varunits,varlongname,jm,im,time_origin,modulo_beg,modulo_end);
% This file create netcdf required by MOM6 forcing
% By Enhui Liao from Laure Resplandy Group in Princeton 
% 2018-10-10

nccreate(filename,'LAT','Dimensions',  {'LAT', jm});
ncwriteatt(filename, 'LAT', 'units', 'degrees_north');
ncwriteatt(filename, 'LAT', 'point_spacing', 'uneven');
ncwriteatt(filename, 'LAT', 'axis', 'Y');

nccreate(filename,'LON','Dimensions',  {'LON', im});
ncwriteatt(filename, 'LON', 'units', 'degrees_east');
ncwriteatt(filename, 'LON', 'modulo', '360.');
ncwriteatt(filename, 'LON', 'point_spacing', 'even');
ncwriteatt(filename, 'LON', 'axis', 'X');

nccreate(filename,varname,'Dimensions',{'LON',im,'LAT',jm,'TIME', inf});
ncwriteatt(filename, varname, 'missing_value', -1.e+34);
ncwriteatt(filename, varname, 'FillValue', -1.e+34);
ncwriteatt(filename, varname, 'long_name', varlongname);
ncwriteatt(filename, varname, 'units', varunits);

nccreate(filename,'TIME','Dimensions',{'TIME', inf});
ncwriteatt(filename, 'TIME', 'units', 'days since 1900-01-01 00:00:00');
ncwriteatt(filename, 'TIME', 'axis', 'T');
ncwriteatt(filename, 'TIME', 'bounds', 'TIME_bnds');
ncwriteatt(filename, 'TIME', 'time_origin', time_origin);
ncwriteatt(filename, 'TIME', 'calendar', 'gregorian');
ncwriteatt(filename, 'TIME', 'modulo', ' ');
ncwriteatt(filename, 'TIME', 'modulo_beg', modulo_beg);
ncwriteatt(filename, 'TIME', 'modulo_end', modulo_end);

nccreate(filename,'TIME_bnds','Dimensions',{'bnds', 2, 'TIME', inf});


