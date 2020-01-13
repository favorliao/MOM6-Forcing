function [wsx,wsy,qla,qsn,evp,tu,qu,dtu,dqu,w10n]=bulk_ncar(us,vs,sat,qar,slp,sst,altu,altt,altq);
% This function is used to change the air temperature and specific humidity from 
% 2m (Observation height of air temperature) to 10m (Observation height of wind)
% Meanwhile, to compute the wind stress, latent and sensible heat flux
% rewriten from bulk_ncar.f90 by Enhui Liao from Laure Resplandy Group in Princeton University
% bulk_ncar.f90 is provided by Hiroyuki from JRA; 
% http://amaterasu.ees.hokudai.ac.jp/~tsujino/JRA55-do-suppl/programs/bulk-ncar.F90
%!   Subroutine bulk
%!
%!     Calculate Latent/Sensible heat flux based on
%!        the bulk formula by Large and Yeager (2004; 2009).
%!
%!   About CPP directives
%!     - LYCOEF    : Compute properties of moist air based on Large and
%!     Yeager (2004;2009)
%!     - CALHEIGHT : Estimate the height of the bottom level (where
%!     surface variables are defined)
%!                    of JMA global spectral models (This option will not
%!q_surface-q_2m!                    be usually used.)
%!
%!   Input:
%!     us   : zonal wind speed (m s-1)
%!     vs   : meridional wind speed (m s-1)
%!     sat  : Surface air temperature (degree Celsius)
%!     qar  : Specific humidity of surface air (kg kg-1)
%!     slp  : Sea level pressure (hPa) %slp units is hPa
%!     sst  : Sea surface temperature (degree Celsius)
%!     altu : Observation height of wind (m) %1D
%!     altt : Observation height of air temperature (m) %1D
%!     altq : Observation height of specific humidity (m) %1D
%!   Output:
%!     wsx  : Zonal wind stress (N m-2)
%!     wsy  : Meridional wind stress (N m-2)
%!     qla  : Latent heat flux (W m-2)
%!     qsn  : Sensible heat flux (W m-2)
%!     evp  : Evaporation (kg m-2 s-1)
%!     tu   : Air temperature at the height where wind is observed
%!     (degree Celsius)
%!     qu   : Specific humidity at the height where wind is observed (kg
%!     kg-1)
%!     dtu  : Difference of temperature between air temperature and sea
%!     surface temperature (degree Celsius)
%!     dqu  : Difference of specific humidity between air specific
%!     humidity
%!          :  and saturation specific humidity at the sea surface  (kg
%!          kg-1)
%!     w10n : Equivalent neutral wind speed at 10 m (m s-1)
if nargin ~=9
   error('bulk_ncar.m: Must pass 9 parameters')
end %if

ro=1.036;         %rho ocean, sea water density of MRI.com; g/cm3
grav=981.0;       %acceleration due to gravity
rhoa=1.22;        %L-Y air density
tab=273.15;       %absolute temperature
sll_cnst=2.5e6;   %L-Y, constant
mwwater=18.016;   %molecular weight of water vapor
mwair=28.966;     %molecular weight of air
rhoa_mks=rhoa;    %air density (kg/m3)
ro0mks=ro*1e3;    %ocean water density, 1036kg/m3
grav_mks=grav*1e-2; %gravity
gasr=287.04;
wvmin=0.3;           %floor on wind speed m/s
agill1 = 0.7859;   agill2 = 0.03477; agill3 = 0.00412; %! saturation vapor pressure
bgill0 = 0.98;     bgill1 = 1.0e-6;  bgill2 = 4.5; bgill3 = 0.0006; %vapor pressure
% reduction and correction for saturation vapor pressure over sea water
cgill1 = 2.5008e6; cgill2 = -2.3e3; %for latend heat of vaporization
dgill1 = 1004.6;   dgill2 = 0.8735; %for specific heat (Gill 1982, P.43)
karman=0.4;                          %von Karman constant
bgsm=0.995;
n_itts=5;
inc_ratio=1e-4;
eps_air=mwwater/mwair;
tvq_air=1.0/eps_air-1;

wdv=(us.^2+vs.^2).^0.5; %scalar wind speed m/s
slpres=slp;     %make sure slp units is hPa
qatmos=qar;
satmos=sat+tab; %t_2m=t_2m+273.15

calheight=0;
if calheight==1;
   pfull=slp*1e2; %pfull: sea-level pressure with units of Pa
   bgsmr=1-bgsm;
   hl1=(bgsmr*(1+log(pfull))+bgsm*log(bgsm))/bgsmr;
   psurf=exp(hl1); %psurf: pressure at the first level (Pa)
   height=(log(psurf)-log(pfull)).*(gasr*(1+tvq_air*qatmos)).*satmos/grav_mks;
   altu=height;altt=height;altq=height; %1D array, height
% satmos: potential temperature referred to sea level   
   satmos=satmos.*((psurf./pfull).^(2.0/7.0));
   slpres=psurf*1e-2; %units: hPa
end
tssurf=sst;
dtemp=tssurf+tab-satmos;
% calculation of saturation specific humidity at the sea surface
sll=cgill1+cgill2*tssurf;
ogcm_lycoef=1;
if ogcm_lycoef==1;
   rhoair=((slpres*1e2/gasr)./satmos)./(1.0+tvq_air*qatmos);
   qs=0.98*6.40380e5*exp(-5107.4./(tssurf+tab))./rhoair;
%!qs specific humidity at sea surface, L-Y 2004, eq. 5, q=rho_air^-1*q1*exp(q2/sst);
else
   hl1=10.0.^((agill1+agill2*tssurf)./(1+agill3*tssurf));
   es=hl1*bgill0.*(1.0+bgill1*slpres.*(bgill2+bgill3*tssurf.^2));
   qs=eps_air*es./(slpres-(1-eps_air)*es);
end
dqr=qs-qatmos;

wv=wdv;
wv(wv<wvmin)=wvmin; % wv = max(wv, wvmin)
tv=satmos.*(1.0+tvq_air*qatmos);
qatmosu=qatmos;
wv10n=wv;

hl1=(2.7./wv10n+0.142+0.0764*wv10n-3.14807e-10*(wv10n.^6))/1e3; %LY2009 eqn. 11a
cdn10=nan(size(wv10n));
cdn10(wv10n<33) =hl1(wv10n<33);
cdn10(wv10n>=33)=2.34e-3;
cdn10_rt=cdn10.^0.5;

cen10=34.6*cdn10_rt/1e3;

ctn10=nan(size(wv10n));
ctn10(-dtemp>=0)=(18/1e3).*cdn10_rt(-dtemp>=0);
ctn10(-dtemp<0)=(32.7/1e3).*cdn10_rt(-dtemp<0);
cdn=cdn10;
ctn=ctn10;
cen=cen10;
cdn_prv=cdn;

%loop_adjust
for n=1:n_itts
    cd_rt=cdn.^0.5;
    ustar=cd_rt.*wv;
    tstar=(ctn./cd_rt).*(-dtemp);
    qstar=(cen./cd_rt).*(-dqr);
% wind velocity
    bstar=grav_mks.*(tstar./tv+qstar./(qatmosu+1.0/tvq_air));
    zetau=(karman*bstar*altu)./(ustar.*ustar);
    zetau_t=abs(zetau);
    zetau_t(zetau_t>=10)=10;
    zetau_t(zetau<0)=-1*zetau_t(zetau<0);
    zetau=zetau_t;
    x2=(abs(1-16*zetau)).^0.5;
    x2(x2<1)=1;
    x=x2.^0.5;
    psi_mu=nan(size(zetau));
    psi_mu(zetau>0)=-5.0*zetau(zetau>0);
    psi_mu(zetau<0)=log((1.0+2.0*x(zetau<0)+x2(zetau<0)).*(1.0+x2(zetau<0))/8.0)-2.0*(atan(x(zetau<0))-atan(1.0));
    psi_hu=nan(size(zetau));
    psi_hu(zetau>0)=-5.0*zetau(zetau>0);
    psi_hu(zetau<0)=2.0*log((1.0+x2(zetau<0))/2.0);
% air temperature    
    zetat=(karman*bstar*altt)./(ustar.*ustar);
    zetat_t=abs(zetat);
    zetat_t(zetat_t>=10)=10;
    zetat_t(zetat<0)=-1*zetat_t(zetat<0);
    zetat=zetat_t;
    x2=(abs(1-16*zetat)).^0.5;
    x2(x2<1)=1;
    x=x2.^0.5;
    psi_mt=nan(size(zetat));
    psi_mt(zetat>0)=-5.0*zetat(zetat>0);
    psi_mt(zetat<0)=log((1.0+2.0*x(zetat<0)+x2(zetat<0)).*(1.0+x2(zetat<0))/8.0)-2.0*(atan(x(zetat<0))-atan(1.0));
    psi_ht=nan(size(zetat));
    psi_ht(zetat>0)=-5.0*zetat(zetat>0);
    psi_ht(zetat<0)=2.0*log((1.0+x2(zetat<0))/2);
% air humidity
    zetaq=(karman*bstar*altq)./(ustar.*ustar);
    zetaq_t=abs(zetaq);
    zetaq_t(zetaq_t>=10)=10;
    zetaq_t(zetaq<0)=-1*zetaq_t(zetaq<0);
    zetaq=zetaq_t;
    x2=(abs(1-16*zetaq)).^0.5;
    x2(x2<1)=1;
    x=x2.^0.5;
    psi_mq=nan(size(zetaq));
    psi_mq(zetaq>0)=-5.0*zetaq(zetaq>0);
    psi_mq(zetaq<0)=log((1.0+2.0*x(zetaq<0)+x2(zetaq<0)).*(1.0+x2(zetaq<0))/8.0)-2.0*(atan(x(zetaq<0))-atan(1.0));
    psi_hq=nan(size(zetaq));
    psi_hq(zetaq>0)=-5.0*zetaq(zetaq>0);
    psi_hq(zetaq<0)=2.0*log((1.0+x2(zetaq<0))/2);

% re-evaluation
    wv10n=wv./(1.0+cdn10_rt.*(log(altu/10.0)-psi_mu)/karman);
    wv10n(wv10n<wvmin)=wvmin; % 0.3m/s floor on wind
    satmosu=satmos-tstar.*(log(altt/altu)+psi_hu-psi_ht)/karman;
    qatmosu=qatmos-qstar.*(log(altq/altu)+psi_hu-psi_hq)/karman;
    tv=satmosu.*(1.0+tvq_air*qatmosu);

    hl1=(2.7./wv10n+0.142+0.0764*wv10n-3.14807e-10*(wv10n.^6))/1e3; %LY2009 eqn. 11a
    cdn10=nan(size(wv10n));
    cdn10(wv10n<33) =hl1(wv10n<33);
    cdn10(wv10n>=33)=2.34e-3;
    cdn10_rt=cdn10.^0.5;
    
    cen10=34.6*cdn10_rt/1e3;

    ctn10=nan(size(wv10n));
    ctn10(zetau>0)=(18/1e3).*cdn10_rt(zetau>0);
    ctn10(zetau<0)=(32.7/1e3).*cdn10_rt(zetau<0);

    zrough=10.0*exp(-karman./cdn10_rt);

    xx=(log(altu/10.0)-psi_mu)/karman;
    cdn=cdn10./((1.0+cdn10_rt.*xx).^2);
    xx=(log(altu/10.0)-psi_hu)/karman;
    ctn=ctn10./(1.0+ctn10.*xx./cdn10_rt).*((cdn./cdn10).^0.5); %ctn transfer coefficient 
    cen=cen10./(1.0+cen10.*xx./cdn10_rt).*((cdn./cdn10).^0.5);
    
    dtemp=tssurf+tab-satmosu;
    dqr=qs-qatmosu;

    rhoair=((slpres*1e2)./tv)/gasr;
    test=abs(cdn-cdn_prv)./(cdn+1e-8);
    cdn_prv=cdn;
    ind=find(test<inc_ratio);
    if n==1;
       ind_sum=zeros(size(test));
       cen_out=nan(size(test));
       cdn_out=cen_out;ctn_out=cen_out;
       rhoair_out=cen_out;satmosu_out=cen_out;
       qatmosu_out=cen_out;dtemp_out=cen_out;
       dqr_out=cen_out;wv10n_out=cen_out;
    end
    if ~isempty(ind) & n<=4
       cen_out(ind_sum==0 & test<inc_ratio)=cen(ind_sum==0 & test<inc_ratio);
       cdn_out(ind_sum==0 & test<inc_ratio)=cdn(ind_sum==0 & test<inc_ratio);
       ctn_out(ind_sum==0 & test<inc_ratio)=ctn(ind_sum==0 & test<inc_ratio);
       rhoair_out(ind_sum==0 & test<inc_ratio)=rhoair(ind_sum==0 & test<inc_ratio);
       satmosu_out(ind_sum==0 & test<inc_ratio)=satmosu(ind_sum==0 & test<inc_ratio);
       qatmosu_out(ind_sum==0 & test<inc_ratio)=qatmosu(ind_sum==0 & test<inc_ratio);
       dtemp_out(ind_sum==0 & test<inc_ratio)=dtemp(ind_sum==0 & test<inc_ratio);
       dqr_out(ind_sum==0 & test<inc_ratio)=dqr(ind_sum==0 & test<inc_ratio);
       wv10n_out(ind_sum==0 & test<inc_ratio)=wv10n(ind_sum==0 & test<inc_ratio);
       ind_sum(test<inc_ratio)=ind_sum(test<inc_ratio)+1;
    else
       cen_out(ind_sum==0)=cen(ind_sum==0);
       cdn_out(ind_sum==0)=cdn(ind_sum==0);
       ctn_out(ind_sum==0)=ctn(ind_sum==0);
       rhoair_out(ind_sum==0)=rhoair(ind_sum==0);
       satmosu_out(ind_sum==0)=satmosu(ind_sum==0);
       qatmosu_out(ind_sum==0)=qatmosu(ind_sum==0);
       dtemp_out(ind_sum==0)=dtemp(ind_sum==0);
       dqr_out(ind_sum==0)=dqr(ind_sum==0);
       wv10n_out(ind_sum==0)=wv10n(ind_sum==0);
    end
end

tu=satmosu_out-tab;
qu=qatmosu_out;    

dtu=dtemp_out;
dqu=dqr_out;
w10n=wv10n_out;

%evaluate: evaporation, latent heat, sensible heat
if ogcm_lycoef==1;
   cpa=1000.5;
   evp=rhoair.*wv.*cen.*dqr; %kg/m2/s
   qla=(-rhoair*sll_cnst).*cen.*dqr.*wv; %w/m2
   qsn=-rhoair.*cpa.*ctn.*dtemp.*wv;
   cdt=rhoair.*cdn.*wv;
else
   cpa=dgill1*(1.0+dgill2*qu);
   evp=rhoair.*wv.*cen.*dqr;
   qla=(-rhoair.*sll).*cen.*dqr.*wv;
   qsn=(-rhoair.*cpa).*ctn.*dtemp.*wv;
   cdt=rhoair.*cdn.*wv;
end

%calculate wind stress
cdn=cdt;
wsx=cdn.*us;
wsy=cdn.*vs;

















