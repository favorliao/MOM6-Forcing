clear
close all

fdir='/tigress/GEOCLIM/LRGROUP/datasets/dfs_forcing_drowned/data_dfs_mom6/';
fall_u=dir([fdir,'dfs5.2_mom6_u10*.nc']);
fall_v=dir([fdir,'dfs5.2_mom6_v10*.nc']);
fall_q=dir([fdir,'dfs5.2_mom6_q2*.nc']);
fall_t=dir([fdir,'dfs5.2_mom6_t2*.nc']);
fall_slp=dir([fdir,'dfs5.2_mom6_slp*.nc']);
fall_sst=dir([fdir,'regrid_3hrs_sst_COBESST_yearly*.nc']);

fm_u=length(fall_u);
altu=10;
altt=2;
altq=2;
for i=1:fm_u
    file_u=[fall_u(i).folder,'/',fall_u(i).name];
    file_v=[fall_v(i).folder,'/',fall_v(i).name];
    file_q=[fall_q(i).folder,'/',fall_q(i).name];
    file_t=[fall_t(i).folder,'/',fall_t(i).name];
    file_slp=[fall_slp(i).folder,'/',fall_slp(i).name];
    file_sst=[fall_sst(i).folder,'/',fall_sst(i).name];
    disp('Reading input data ...');tic;
    disp([file_u]); disp([file_v]); disp([file_q]); disp([file_t]); disp([file_slp]); disp([file_sst]);tic
    
    u10=double(ncread(file_u,'u10'));
    time_u=double(ncread(file_u,'TIME'))+datenum(1900,1,1);
    time_bnds=double(ncread(file_u,'TIME_bnds'))+datenum(1900,1,1);
    lm_u10=length(time_u);
    if i==1;
       lon=double(ncread(file_u,'LON'));
       lat=double(ncread(file_u,'LAT'));
    end

    v10=double(ncread(file_v,'v10'));
    q2=double(ncread(file_q,'q2'));
    t2=double(ncread(file_t,'t2'));
    slp=double(ncread(file_slp,'slp'));
    sst=double(ncread(file_sst,'tos'));
    q10=nan(size(q2));
    t10=nan(size(t2));
 
    toc;disp('Computing the q10 and t10 ...');tic;
    for l=1:lm_u10
        if (mod(l,100)==0);toc;disp([num2str(l),'/',num2str(lm_u10)]);tic;end
        us=u10(:,:,l);vs=v10(:,:,l);
        sat=t2(:,:,l)-273.15;qar=q2(:,:,l);
        slp2d=slp(:,:,l)/100; %slp units is hPa;
        sst2d=sst(:,:,l);
        [wsx,wsy,qla,qsn,evp,tu,qu,dtu,dqu,w10n]=bulk_ncar(us,vs,sat,qar,slp2d,sst2d,altu,altt,altq);
        qu(isnan(sst2d))=qar(isnan(sst2d));
        tu(isnan(sst2d))=sat(isnan(sst2d));
        q10(:,:,l)=qu;
        t10(:,:,l)=tu+273.15;
    end
    clear u10 v10 t2 q2 slp sst
    year=1958+i-1;jm=length(lat);im=length(lon);
    time_origin=datestr(time_u(1),'dd-mmm-yyyy');
    modulo_beg=datestr(time_bnds(1,1),'yyyy-mm-dd HH:MM:SS');
    modulo_end=datestr(time_bnds(2,end),'yyyy-mm-dd HH:MM:SS');
%write q10
    varname='q10';varunits='kg/kg';varlongname='ERA corrected for DFS5 specific humidity at 10m (land is still 2m)';
    filename_new=[fdir,'dfs5.2_mom6_',varname,'_',num2str(year),'_10oct2018.nc'];
    toc;disp(['Creating and writing dfs file for MOM6: ', filename_new]);tic;
    disp(['Time_orgin: ',time_origin,' modulo_beg: ',modulo_beg,' modulo_end: ',modulo_end])
    delete(filename_new)
    create_mom6nc(filename_new,varname,varunits,varlongname,jm,im,time_origin,modulo_beg,modulo_end);
    ncwrite(filename_new,'LON', lon);
    ncwrite(filename_new,'LAT', lat);
    ncwrite(filename_new,'TIME',time_u-datenum(1900,1,1));
    ncwrite(filename_new,'TIME_bnds', time_bnds-datenum(1900,1,1));
    ncwrite(filename_new,varname, q10);
    ncwriteatt(filename_new,'/','data_source','The file is from DFS5.2 (https://ige-meom-opendap.univ-grenoble-alpes.fr/thredds/catalog/meomopendap/extract/FORCING_ATMOSPHERIQUE/DFS5.2/ALL/catalog.html). The file is made for MOM6 input forcing');
    ncwriteatt(filename_new,'/','data_creator','Enhui Liao from Laure Resplandy Group in Princeton');
    ncwriteatt(filename_new,'/','creation_date',datestr(now));
    toc;

%write t10
    varname='t10';varunits='K';varlongname='ERA corrected for DFS5 air temperature at 10m (land is till 2m)';
    filename_new=[fdir,'dfs5.2_mom6_',varname,'_',num2str(year),'_10oct2018.nc'];
    toc;disp(['Creating and writing dfs file for MOM6: ', filename_new]);tic;
    disp(['Time_orgin: ',time_origin,' modulo_beg: ',modulo_beg,' modulo_end: ',modulo_end])
    delete(filename_new)
    create_mom6nc(filename_new,varname,varunits,varlongname,jm,im,time_origin,modulo_beg,modulo_end);
    ncwrite(filename_new,'LON', lon);
    ncwrite(filename_new,'LAT', lat);
    ncwrite(filename_new,'TIME',time_u-datenum(1900,1,1));
    ncwrite(filename_new,'TIME_bnds', time_bnds-datenum(1900,1,1));
    ncwrite(filename_new,varname, t10);
    ncwriteatt(filename_new,'/','data_source','The file is from DFS5.2 (https://ige-meom-opendap.univ-grenoble-alpes.fr/thredds/catalog/meomopendap/extract/FORCING_ATMOSPHERIQUE/DFS5.2/ALL/catalog.html). The file is made for MOM6 input forcing');
    ncwriteatt(filename_new,'/','data_creator','Enhui Liao from Laure Resplandy Group in Princeton');
    ncwriteatt(filename_new,'/','creation_date',datestr(now));
    toc;

end



