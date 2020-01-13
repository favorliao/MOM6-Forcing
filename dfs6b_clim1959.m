clear 
close all

addpath(genpath('/tigress/enhuil/matlab'));
%varall=("precip" "q10" "radlw" "radsw" "slp" "snow" "t10" "u10" "v10")
fvarname={'precip','q10','radlw','radsw','slp','snow','t10','u10','v10'};
 varname={'precip','q10','radlw','radsw','slp','snow','t10','u10','v10'};
foridir='/tigress/GEOCLIM/LRGROUP/datasets/dfs_forcing_drowned/res_dfs_mom6_3yrs';
fcopydir='/tigress/GEOCLIM/LRGROUP/datasets/dfs_forcing_drowned/res_dfs_mom6_clim/res_step1sh';
foutdir='/tigress/GEOCLIM/LRGROUP/datasets/dfs_forcing_drowned/res_dfs_mom6_clim';
%fcopydir='/tigress/GEOCLIM/LRGROUP/datasets/JRA55-do-v1.3merge/data_origin';

tb=datenum(1959,1,1);
te=datenum(1960,1,1);

for vi=1:1
    fvar=char(fvarname(vi));var=char(varname(vi));
    file=dir([foridir,'/*',fvar,'*1959*.nc']);%dfs5.2_mom6_precip_1959_10oct2018.nc
    filename=[file.folder,'/',file.name];
    disp(['Reading file: ',filename]);tic;
    data=ncread(filename,var);toc;
    [im,jm,lm]=size(data);
    time=ncread(filename,'TIME')+datenum(1900,1,1);
    time_bnds=ncread(filename,'TIME_bnds')+datenum(1900,1,1);
    if lm==1096;
         lm_new=365;
         %time=daynoleap2datenum(ncread(filename,'Time'),1900);
         %time_bnds=daynoleap2datenum(ncread(filename,'Time_bnds'),1900);
    elseif lm==8768
         lm_new=2920;
         %time=daynoleap2datenum(ncread(filename,'time'),1900);
         %time_bnds=daynoleap2datenum(ncread(filename,'time_bnds'),1900);
    else
        disp(['time length is wrong']);break;
    end
    ind_decb=find(time<tb & time>=tb-15);
    ind_mid=find(time>=tb & time<te-15);
    ind_dece=find(time>=te-15 & time<te);
    data_new=nan(im,jm,lm_new);
%weight
    lwn=length(ind_dece);
    lwn2=length(ind_decb);
    if lwn==lwn2
       lw=repmat(1,[1,1,lwn]);lwr=lw;
       lw(1,1,:)=[0:1/(lwn-1):1]; %linear weight
       lwr(1,1,:)=[1:-1/(lwn-1):0]; %linear weight reverse
       lw2=repmat(lw,[im,jm,1]);
       lwr2=repmat(lwr,[im,jm,1]);
    else
       disp(['lwn and lwn2 not equal']);break;
    end

    data_decb=data(:,:,ind_decb);
    data_mid =data(:,:,ind_mid);
    data_dece=data(:,:,ind_dece);
    data_decenew=data_decb.*lw2+data_dece.*lwr2;
    data_new(:,:,1:length(ind_mid))=data_mid;
    data_new(:,:,length(ind_mid)+1:lm_new)=data_decenew;
    time_new=[time(ind_mid);time(ind_dece)];
    time_bnds_new=[time_bnds(:,ind_mid),time_bnds(:,ind_dece)];
%write new netcdf
    disp(['Writing netcdf ...']);tic;
    file_copy=dir([fcopydir,'/',fvar,'.DFS.clim.1959.nc']);%q_10.1959.18Oct2017.nc
    filename_copy=[file_copy.folder,'/',file_copy.name];
    filename_out=[foutdir,'/',fvar,'.DFS.clim.1959.nc'];
    copyfile(filename_copy,filename_out)
    time_new=time_new-datenum(1959,1,1);
    time_bnds_new=time_bnds_new-datenum(1959,1,1);
    ncwrite(filename_out,'TIME',time_new);
    ncwrite(filename_out,'TIME_bnds',time_bnds_new);
    ncwrite(filename_out,var,data_new);toc;
%    clear lw lwr lw2 lwr2 data_new time_new data_decb data_mid data_dece
end

i0=320;j0=160;
d0=squeeze(data(i0,j0,:));
time0=time;
d1=squeeze(data_new(i0,j0,:));
time1=time(1:2920);
time2=time_new;
tbf=datenum(1958,12,15);
tef=datenum(1959,1,5);
figure(1);clf
plot(time0,d0);hold on
plot(time1,d1,'--.');
plot(time2,d1,'--.');
tk=datelist(tbf,tef,1:12,1:1:31);
set(gca,'xlim',[tbf tef],'xtick',tk)
datetick('x','mm-dd','keeplimits','keepticks')
grid on

%a(2920-119:2920).*coefr+a(1:120).*coef;a(2920-119:2920).*coef+a(1:120).*coefr];


