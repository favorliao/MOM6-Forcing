clear
close all

%u10 v10 t2 q2 precip snow radsw radlw

var_all=     {'u10','v10','t2','q2'   ,'precip' ,'snow'    ,'radsw','radlw','slp'};
varunits_all={'m/s','m/s','K' ,'kg/kg','kg/m2/s','kg/m2/s' ,'w/m2' ,'w/m2' ,'Pa'};
varlongname_all={'Zonal wind speed at 10m',...
                 'Meridional wind speed at 10m',...
                 'ERA-interim arctic correction (Brodeau et al. 2007) temperature at 2m and change to 10m by bulk_ncar'...
                 'ERA-interim corrected for DFS5 (Brodeau et al. 2007) specific humidity at 2m and change to 10m by bulk_ncar'...
                 'Total precipitation',...
                 'Snowfall (convective + stratiform)',...
                 'Surface thermal radiation downwards',...
                 'Surface solar radiation downwards',...
                 'Mean sea level pressure'};
dirnew='/tigress/GEOCLIM/LRGROUP/datasets/dfs_forcing_drowned/data_dfs_mom6';

for varnum=6:6
    varname=char(var_all(varnum));
    varunits=char(varunits_all(varnum));
    varlongname=char(varlongname_all(varnum));
    fall=dir(['/tigress/GEOCLIM/LRGROUP/datasets/dfs_forcing_drowned/data_dfs/drowned_',varname,'_DFS5.2_y*.nc']);
    if varnum==9;
       fall=dir('/tigress/GEOCLIM/LRGROUP/datasets/dfs_forcing_drowned/slp_ecmwf/regrid_yearly_slp_*.nc');
    end
    fm=length(fall);
    for fnum=1:fm
        clear time_bnds
        filename=[fall(fnum).folder,'/',fall(fnum).name];
        year=1958+fnum-1;

        disp(['Reading dfs file: ', filename]);
        disp(['Year: ',num2str(year),' Var: ',varname,' Units: ',varunits,' Longname: ',varlongname]);tic;
        time0=double(ncread(filename,'time'));
        if fnum==1;
           time_constant=double(ncread(filename,'time'));
        end
        if varnum==9;
           %time in slp nc is different with dfs forcing
           time=double(ncread(filename,'time'))/24.0+datenum(1900,1,1);
           data=double(ncread(filename,'msl'));
           lon=double(ncread(filename,'LON'));
           lat=double(ncread(filename,'LAT'));
           dt_bnds=1.5;
        else
           if length(time0)==365;
              time=time_constant+datenum(year,1,1)-0.5; %date of daily data is 1,2,3...,365
           elseif length(time0)==366;
              time=[time_constant;366]+datenum(year,1,1)-0.5;
           elseif length(time0)==2920;
              time=time_constant+datenum(year,1,1); %date of daily data is 0.125,0.25...,365
           elseif length(time0)==2928;
              time=[time_constant;365.125;365.25;365.375;365.5;365.625;365.75;365.875;366]+datenum(year,1,1);
           end  
           year2=str2num(filename(end-6:end-3));
           if datenum(year2,1,1)>time(1) | datenum(year2+1,1,1)<time(end)
              disp('time in the nc is wrong! Please check!'); break;
           end
           data=double(ncread(filename,varname));
           if varnum==6; data(data<0)=0; end %snow in dfs in some regions has negative value which induce MOM6 blow up
           if length(time)>=1000;
              time(2:end)=time(1:end-1);
              time(1)=time(2)-3/24;
              data(:,:,2:end)=data(:,:,1:end-1);
              data(:,:,1)    =data(:,:,2);
              dt_bnds=1.5;
           else
              dt_bnds=12;
           end
           lon=double(ncread(filename,'lon0'));
           lat=double(ncread(filename,'lat0'));
           if lat(1)>lat(end); %latitude is from 90-->-90 in dfs, I reverse it here
               lat=flip(lat);
               data=flip(data,2);
           end
        end 
        lm=length(time);
        [im,jm,lm]=size(data);
        time_bnds(1,1:lm)=time-dt_bnds/24;
        time_bnds(2,1:lm)=time+dt_bnds/24;
        time_origin=datestr(time_bnds(1,1),'dd-mmm-yyyy');
        modulo_beg=datestr(time_bnds(1,1),'yyyy-mm-dd HH:MM:SS');
        modulo_end=datestr(time_bnds(2,end),'yyyy-mm-dd HH:MM:SS');

        filename_new=[dirnew,'/dfs5.2_mom6_',varname,'_',num2str(year),'_10oct2018.nc'];
        toc;disp(['Creating and writing dfs file for MOM6: ', filename_new]);tic;
        disp(['Time_orgin: ',time_origin,' modulo_beg: ',modulo_beg,' modulo_end: ',modulo_end])
        delete(filename_new)
        create_mom6nc(filename_new,varname,varunits,varlongname,jm,im,time_origin,modulo_beg,modulo_end);
        ncwrite(filename_new,'LON', lon);
        ncwrite(filename_new,'LAT', lat);
        ncwrite(filename_new,'TIME',time-datenum(1900,1,1));
        ncwrite(filename_new,'TIME_bnds', time_bnds-datenum(1900,1,1));
        ncwrite(filename_new,varname, data);
        ncwriteatt(filename_new,'/','data_source','The file is from DFS5.2 (https://ige-meom-opendap.univ-grenoble-alpes.fr/thredds/catalog/meomopendap/extract/FORCING_ATMOSPHERIQUE/DFS5.2/ALL/catalog.html). The file is made for MOM6 input forcing');
        ncwriteatt(filename_new,'/','data_creator','Enhui Liao from Laure Resplandy Group in Princeton');
        ncwriteatt(filename_new,'/','creation_date',datestr(now));
        toc;
     end
end



