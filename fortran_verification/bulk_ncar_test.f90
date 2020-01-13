! -*-F90-*-
!------------------------ bulk-ncar.F90 ----------------------------
!   Subroutine bulk
!
!     Calculate Latent/Sensible heat flux based on 
!        the bulk formula by Large and Yeager (2004; 2009).
!
!   About CPP directives
!     - LYCOEF    : Compute properties of moist air based on Large and
!     Yeager (2004;2009)
!     - CALHEIGHT : Estimate the height of the bottom level (where
!     surface variables are defined)
!                    of JMA global spectral models (This option will not
 !q_surface-q_2m!                    be usually used.)
!
!   Input:
!     us   : zonal wind speed (m s-1)
!     vs   : meridional wind speed (m s-1)
!     sat  : Surface air temperature (degree Celsius)
!     qar  : Specific humidity of surface air (kg kg-1)
!     wdv  : Scalar wind speed (m s-1)
!     slp  : Sea level pressure (hPa)
!     sst  : Sea surface temperature (degree Celsius)
!     imx,jmx : Size of the arrays (integer)
!     atexl: land-sea mask (0:land, 1:sea)
!     altu : Observation height of wind (m)
!     altt : Observation height of air temperature (m)
!     altq : Observation height of specific humidity (m)
!
!   Output:
!     wsx  : Zonal wind stress (N m-2)
!     wsy  : Meridional wind stress (N m-2)
!     qla  : Latent heat flux (W m-2)
!     qsn  : Sensible heat flux (W m-2)
!     evp  : Evaporation (kg m-2 s-1)
!     tu   : Air temperature at the height where wind is observed
!     (degree Celsius)
!     qu   : Specific humidity at the height where wind is observed (kg
!     kg-1)
!     dtu  : Difference of temperature between air temperature and sea
!     surface temperature (degree Celsius)
!     dqu  : Difference of specific humidity between air specific
!     humidity
!          :  and saturation specific humidity at the sea surface  (kg
!          kg-1)
!     w10n : Equivalent neutral wind speed at 10 m (m s-1)
!
program test_bulk
integer,parameter:: imx=5,jmx=1
real(8):: us(imx,jmx),vs(imx,jmx),sat(imx,jmx),qar(imx,jmx)
real(8):: wdv(imx,jmx),slp(imx,jmx),atexl(imx,jmx),sst(imx,jmx)
real(8):: altu,altt,altq
real(8):: wsx(imx,jmx),wsy(imx,jmx),qla(imx,jmx),qsn(imx,jmx)
real(8):: evp(imx,jmx),tu(imx,jmx),qu(imx,jmx),dtu(imx,jmx)
real(8):: dqu(imx,jmx),w10n(imx,jmx)

altu=10
altt=2
altq=2
do i=1,imx
   us(i,1)=0.0+(i+1)*2.0
   vs(i,1)=0.0+(i+1)*2.0
   sat(i,1)=(i+1)*2
   qar(i,1)=0.0+(i+1)*0.1
   slp(i,1)=(1004.0+10.0*(i+1))
   sst(i,1)=-20.0+(i+1)*16.0
   atexl(i,1)=1
   wdv(i,1)=(us(i,1)**2+vs(i,1)**2)**0.5
enddo
!print*,us,vs,sat,qar,slp,sst,atexl,wdv
call bulk(wsx,wsy,qla,qsn,evp,tu,qu,dtu,dqu,w10n,&
    & us,vs,sat,qar,wdv,slp,sst,imx,jmx,atexl,altu,altt,altq)

print*,'us= ',us,' vs= ',vs
print*,'wsx= ',wsx,' wsy= ',wsy

print*,'sat_2m= ',sat,' sat_10m= ',tu
print*,'qar_2m= ',qar,' qar_10m= ',qu

print*,'wdv= ',wdv,' w10n= ',w10n

end





subroutine bulk(wsx,wsy,qla,qsn,evp, &
     & tu,qu,dtu,dqu,w10n,&
     & us,vs,sat,qar,wdv,slp,sst,&
     & imx,jmx,atexl,altu,altt,altq)

  implicit none

  integer(4), intent(in) :: imx, jmx
  real(8), intent(out) :: wsx(imx,jmx), wsy(imx,jmx)
  real(8), intent(out) :: qla(imx,jmx), qsn(imx,jmx)
  real(8), intent(out) :: evp(imx,jmx)
  real(8), intent(out) :: tu(imx,jmx), qu(imx,jmx)
  real(8), intent(out) :: dtu(imx,jmx), dqu(imx,jmx)
  real(8), intent(out) :: w10n(imx,jmx)
  real(8), intent(in)  :: us (imx,jmx), vs (imx,jmx)
  real(8), intent(in)  :: sat(imx,jmx), qar(imx,jmx)
  real(8), intent(in)  :: wdv(imx,jmx), slp(imx,jmx)
  real(8), intent(in)  :: sst(imx,jmx)
  real(8), intent(in)  :: atexl(imx,jmx)
  real(8), intent(in)  :: altu, altt, altq 

  ! Constants
  ! [cgs]

  real(8), parameter :: ro   = 1.036d0  ! sea water density of MRI.COM
  real(8), parameter :: grav = 981.0d0  ! acceleration due to gravity 

  ! ...

  real(8), parameter :: rhoa = 1.22d0  ! L-Y !air density
  real(8), parameter :: tab = 273.15d0
!#ifdef OGCM_LYCOEF
  real(8), parameter :: sll_cnst = 2.5d6 ! L-Y
!#endif /* OGCM_LYCOEF */

  real(8), parameter :: mwwater = 18.016d0  ! molecular weight of water vapor
  real(8), parameter :: mwair = 28.966d0    ! molecular weight of air

  ! [MKS]

  real(8), parameter :: rhoa_mks = rhoa     ! density of air (kg/m^3)
  real(8), parameter :: ro0mks = ro * 1.d3
  real(8), parameter :: grav_mks = grav * 1.0d-2
  real(8), parameter :: gasr = 287.04d0

  real(8), parameter :: wvmin = 0.3d0 ! floor on wind speed (m/s)
  real(8), save :: eps_air   ! = 0.62197
  real(8), save :: tvq_air   ! = 0.6078

  ! Gill (1982) Appendix 4
  real(8), parameter :: agill1 = 0.7859d0, agill2 = 0.03477d0, agill3 = 0.00412d0 ! saturation vapor pressure
  real(8), parameter :: bgill0 = 0.98d0 ! 2 % reduction for saturation vapor pressure over sea water
  real(8), parameter :: bgill1 = 1.d-6, bgill2 = 4.5d0, bgill3 = 0.0006d0 ! correction of vapor pressure
  real(8), parameter :: cgill1 = 2.5008d6, cgill2 = -2.3d3   ! for latent heat of vaporization
  real(8), parameter :: dgill1 = 1004.6d0, dgill2 = 0.8735d0 ! for specific heat (Gill 1982, P.43)

  !   Variables

  !   sll: latent heat of vaporization (J/Kg)

  real(8) :: sll, es, qs
  real(8) :: cpa  ! specific heat of air
  !
  real(8) :: dqr, dtemp
  !
  real(8) :: qatmos, satmos, slpres, tssurf
  real(8) :: rhoair
  !
  ! defined at 10 [m]
  !
  real(8) :: cdn10, cen10, ctn10
  real(8) :: cdn10_rt
  real(8) :: wv10n ! neutral 10m wind for evaluating cdn10
  !
  ! defined at velocity level
  !
  real(8) :: qatmosu, satmosu
  real(8) :: cdn, cen, ctn
  real(8) :: cd_rt
  real(8) :: wv
  !
  real(8), parameter :: karman = 0.4    ! von Karman constant
  real(8) :: zrough                     ! roughness length
  real(8) :: stab
  !
  real(8) :: ustar, tstar, qstar, bstar
  real(8) :: zetau, zetat, zetaq
  real(8) :: psi_mu, psi_hu
  real(8) :: psi_mt, psi_ht
  real(8) :: psi_mq, psi_hq
  !
  real(8) :: x, xx, x2, tv
  !
!#ifdef OGCM_CALHEIGHT
!  ! this is for GSM of JMA
!  real(8), parameter :: bgsm = 0.995
!  real(8)            :: bgsmr
!  real(8)            :: psurf, pfull, height
!#endif /* OGCM_CALHEIGHT */
  !
!#ifdef OGCM_BULKITER
  integer(4) :: n
  integer(4), parameter :: n_itts = 5
  real(8), parameter :: inc_ratio = 1.0d-4
  real(8) :: test, cdn_prv
!#endif /* OGCM_BULKITER */
  !
  !  work variables
  !
  integer(4) :: i, j
  real(8) :: hl1
  !
!#ifdef OGCM_TAUBULK
  real(8) :: cdt(imx,jmx)
!#endif /* OGCM_TAUBULK */
  !
  !----------------------------------------------------------------------

  eps_air = mwwater / mwair !molecular weight of water vapor/air~=18/29
  tvq_air = 1.0d0/eps_air - 1.0d0

  wsx(:,:) = 0.0d0
  wsy(:,:) = 0.0d0
  qla(:,:) = 0.0d0
  qsn(:,:) = 0.0d0
  evp(:,:) = 0.0d0
  dtu(:,:) = 0.0d0
  dqu(:,:) = 0.0d0
  w10n(:,:) = 0.0d0
  tu(:,:) = 0.0d0
  qu(:,:) = 0.0d0

!$omp parallel
!$omp do private(j,i,n,slpres,qatmos,satmos,&
!$omp & tssurf, dtemp, sll, es, qs, dqr, wv, tv, wv10n, hl1, &
!$omp & cdn10, cdn10_rt, cen10, stab, ctn10, cdn, ctn, cen, &
!$omp & cd_rt, ustar, tstar, qstar, bstar, zetau, x2, x, xx, &
!$omp & psi_mu, psi_hu, zetat, psi_mt, psi_ht, &
!$omp & zetaq, psi_mq, psi_hq, satmosu, qatmosu, zrough, rhoair, cpa,
!cdn_prv, test)
  do j = 1, jmx
    do i = 1, imx

      if (atexl(i,j) == 1.0d0) then !atexl land-sea mask
        !
        slpres = slp(i,j)
        qatmos = qar(i,j)       !qar=q2m
        satmos = sat(i,j) + tab !sat=t2m
        !
!#ifdef OGCM_CALHEIGHT
        !
        ! !!!!! do not change the line order below !!!!!
        !
        ! psurf : Sea level pressure [Pa]
        ! pfull : Pressure at the first level [Pa]
        !
!        pfull = slpres * 1.0d2
!        bgsmr = 1.0d0 - bgsm
!        hl1 = (bgsmr * (1.0d0 + log(pfull)) + bgsm * log(bgsm)) / bgsmr
!        psurf = exp(hl1)
!        !
!        height = (log(psurf) - log(pfull)) &
!             &  * gasr * (1.0d0 + tvq_air * qatmos) * satmos / grav_mks
!        !
!        altu = height
!        altt = height
!        altq = height
        !
        ! satmos : potential temperature referred to sea level
        !
!        satmos = satmos * ((psurf / pfull)**(2.0d0/7.0d0))
!        slpres = psurf * 1.0d-2
        !
!#endif /* OGCM_CALHEIGHT */
        !
        tssurf = sst(i,j)
        dtemp = tssurf + tab - satmos
        !
        ! Calculation of saturation specific humidity at the sea surface
        !
        sll = cgill1 + cgill2 * tssurf
!#ifdef OGCM_LYCOEF
        rhoair = (slpres * 1.0d2)/gasr/satmos/(1.0d0+tvq_air*qatmos)
        qs = 0.98d0 * 6.40380d5 * exp(-5107.4d0 / (tssurf + tab)) / rhoair
        !qs specific humidity at sea surface, L-Y 2004, eq. 5, q=rho_air-1*q1*exp(q2/sst);
!#else /* OGCM_LYCOEF */
!        hl1 = 10.D0**((agill1 + agill2 * tssurf) / (1.0d0 + agill3 * tssurf))  ! [hPa]
!        es = hl1 * bgill0 * (1.0d0 + bgill1 * slpres * (bgill2 + bgill3 * tssurf**2)) ! [hPa]
!        qs = eps_air * es / (slpres - (1.0d0 - eps_air) * es)
!#endif /* OGCM_LYCOEF */
        dqr = qs - qatmos !q_surface-q_2m
        !
        wv = wdv(i,j) ! wv: wind velocity, when input data is in MKS
        wv = max(wv, wvmin)                  ! 0.3 [m/s] floor on wind
        !
        tv = satmos * (1.0d0 + tvq_air * qatmos)
        qatmosu = qatmos

        wv10n = wv    !10n: neutral stability 10 meter, first guess 10m wind
        !
        !!!!!cdn10 = (2.7d0 / wv10n + 0.142d0 + 0.0764d0 * wv10n) &
        !!!!!     &      / 1.0d3                                    !
        !L-Y eqn. 6a (LY2004)
        hl1 = (2.7d0 / wv10n + 0.142d0 + 0.0764d0 * wv10n - 3.14807d-10 * (wv10n**6)) &
             &      / 1.0d3                                    ! LY2009 eqn. 11a
        cdn10 = (0.5d0 - sign(0.5d0,wv10n-33.0d0)) * hl1 &
             & + (0.5d0 + sign(0.5d0,wv10n-33.0d0)) * 2.34d-3   ! LY2009 eqn. 11b
        cdn10_rt = sqrt(cdn10)
        cen10 = 34.6d0 * cdn10_rt / 1.0d3                      ! L-Y 2004 eqn. 6b
        stab = 0.5d0 + sign(0.5d0,-dtemp)
        ctn10 = (18.0d0 * stab + 32.7d0 * (1.0d0 - stab)) &
             &      * cdn10_rt / 1.0d3                         ! L-Y 2004 eqn. 6c

        cdn = cdn10 !first guess for exchange coeff's at z
                    !transfer coefficient for wind
        ctn = ctn10 !transfer coefficient for air temperature
        cen = cen10 !transfer coefficient for humidity

        cdn_prv = cdn

!!!!!#ifdef OGCM_BULKITER ! Always iterate
        !
        LOOP_ADJUST: do n = 1, n_itts                          !Monin-Obukhov iteration
          !
          cd_rt = sqrt(cdn)
          ustar = cd_rt * wv                                   ! L-Y 2004 eqn. 7a
          tstar = (ctn / cd_rt) * (-dtemp)                     ! L-Y 2004 eqn. 7b
          qstar = (cen / cd_rt) * (-dqr)                       ! L-Y 2004 eqn. 7c, dqr=q_surf-q_2m
          !ustar: turbulent scale for wind velocity
          !tstar,qstar: same for air temperature t and humidity q
          !bstar = grav_mks * &
          !     & ( tstar / tv + qstar / (qatmos + 1.0d0 / 0.6078d0))
          bstar = grav_mks * &
               & ( tstar / tv + qstar / (qatmosu + 1.0d0 / tvq_air))
          !
          ! velocity !!!  level
          !
          zetau = karman * bstar * altu / (ustar * ustar)      ! L-Y 2004 eqn. 8a
          zetau = sign(min(abs(zetau),10.d0),zetau)            !undocumented NCAR
          x2 = sqrt(abs(1.0d0 - 16.0d0 * zetau))               ! L-Y eqn. 8b
          x2 = max(x2,1.0d0)                                   !undocumented NCAR
          x = sqrt(x2)
          if (zetau > 0.0d0) then
            psi_mu = - 5.0d0 * zetau                           ! L-Y eqn. 8c
            psi_hu = - 5.0d0 * zetau                           ! L-Y eqn. 8c
          else
            psi_mu = log((1.0d0 + 2.0d0 * x + x2) &
                 & * (1.0d0 + x2) / 8.0d0) &
                 & - 2.0d0 * (atan(x) - atan(1.0d0))           ! L-Y eqn. 8d
            psi_hu = 2.0d0 * log((1.0d0 + x2) / 2.0d0)         ! L-Y eqn. 8e
          end if
          !
          ! temperature level
          !
          zetat = karman * bstar * altt / (ustar * ustar)      ! L-Y eqn. 8a
          zetat = sign(min(abs(zetat),10.d0),zetat)            ! undocumented NCAR
          x2 = sqrt(abs(1.0d0 - 16.0d0 * zetat))               ! L-Y eqn. 8b
          x2 = max(x2,1.0d0)                                   ! undocumented NCAR
          x = sqrt(x2)
          if (zetat > 0.0d0) then
            psi_mt = - 5.0d0 * zetat                            ! L-Y eqn. 8c
            psi_ht = - 5.0d0 * zetat                            ! L-Y eqn. 8c
          else
            psi_mt = log((1.0d0 + 2.0d0 * x + x2) &
                 & * (1.0d0 + x2) / 8.0d0) &
                 & - 2.0d0 * (atan(x) - atan(1.0d0))            ! L-Y eqn. 8d
            psi_ht = 2.0d0 * log((1.0d0 + x2) / 2.0d0)          ! L-Y eqn. 8e
          end if
          !
          ! humidity !!!water vapor level
          !
          zetaq = karman * bstar * altq / (ustar * ustar)      ! L-Y eqn. 8a
          zetaq = sign(min(abs(zetaq),10.d0),zetaq)            ! undocumented NCAR
          x2 = sqrt(abs(1.0d0 - 16.0d0 * zetaq))               ! L-Y eqn. 8b
          x2 = max(x2,1.0d0)                                   ! undocumented NCAR
          x = sqrt(x2)
          if (zetaq > 0.0d0) then
            psi_mq = - 5.0d0 * zetaq                            ! L-Y eqn. 8c
            psi_hq = - 5.0d0 * zetaq                            ! L-Y eqn. 8c
          else
            psi_mq = log((1.0d0 + 2.0d0 * x + x2) &
                 & * (1.0d0 + x2) / 8.0d0) &
                 & - 2.0d0 * (atan(x) - atan(1.0d0))           ! L-Y eqn. 8d
            psi_hq = 2.0d0 * log((1.0d0 + x2) / 2.0d0)          ! L-Y eqn. 8e
          end if
          !
          ! re-evaluation
          !
          wv10n = wv / (1.0d0 + &
               & cdn10_rt * (log(altu / 10.0d0) - psi_mu) &
               & / karman)                                     ! L-Y 2004 eqn. 9a
          wv10n = max(wv10n, wvmin)             ! 0.3 [m/s] floor on wind
          satmosu = satmos - tstar * &
               & (log(altt / altu) + psi_hu - psi_ht) / karman ! L-Y 2004 eqn. 9b
          !the air temperature at 10m, observed wind u level
          qatmosu = qatmos - qstar * &
               & (log(altq / altu) + psi_hu - psi_hq) / karman ! L-Y 2004 eqn. 9c
          !the specific humidity at 10m, 
          tv = satmosu * (1.0d0 + tvq_air * qatmosu)
        
          !!!!!cdn10 = (2.7d0 / wv10n + 0.142d0 + 0.0764d0 * wv10n) & !
          !!!!!     & / 1.0d3                                       !
          !L-Y eqn. 6a again  (LY2004)
          hl1 = (2.7d0 / wv10n + 0.142d0 + 0.0764d0 * wv10n - 3.14807d-10 * (wv10n**6)) &
               &      / 1.0d3                                    !LY2009 eqn. 11a
          cdn10 = (0.5d0 - sign(0.5d0,wv10n-33.0d0)) * hl1 &
               & + (0.5d0 + sign(0.5d0,wv10n-33.0d0)) * 2.34d-3  !LY2009 eqn. 11b
          cdn10_rt = sqrt(cdn10)                               !
          cen10 = 34.6d0 * cdn10_rt / 1.0d3                    ! L-Y eqn. 6b again
          stab = 0.5d0 + sign(0.5d0,zetau)
          ctn10 = (18.0d0 * stab + 32.7d0 * (1.0d0 - stab)) &
               & * cdn10_rt / 1.0d3                            ! L-Y eqn. 6c again
          zrough = 10.d0 * exp(- karman / cdn10_rt)            ! diagnostic
        
          xx = (log(altu / 10.0d0) - psi_mu) / karman
          cdn = cdn10 / (1.0d0 + cdn10_rt * xx)**2             ! L-Y 2004 10a
          xx = (log(altu / 10.0d0) - psi_hu) / karman
          ctn = ctn10 / (1.0d0 + ctn10 * xx / cdn10_rt) * sqrt(cdn / cdn10)  ! L-Y 2004 10b
          cen = cen10 / (1.0d0 + cen10 * xx / cdn10_rt) * sqrt(cdn / cdn10)  ! L-Y 2004 10c
          !ctn: transfer coefficient for evaporation
          !cen: 
        
          dtemp = tssurf + tab - satmosu
          dqr = qs - qatmosu
          rhoair = (slpres * 1.0d2)/gasr/tv

          test = abs(cdn - cdn_prv) / (cdn + 1.0d-8)

          if (test < inc_ratio) exit LOOP_ADJUST

          cdn_prv = cdn
        end do LOOP_ADJUST
print*,n
        tu(i,j) = satmosu - tab
        qu(i,j) = qatmosu
!#endif /* OGCM_BULKITER */
        dtu(i,j) = dtemp
        dqu(i,j) = dqr
        w10n(i,j) = wv10n

        ! Evaluate
        ! evaporation
        ! latent heat
        ! sensible heat

!#ifdef OGCM_LYCOEF
        cpa = 1000.5d0
        evp(i,j) = rhoair * wv * cen * dqr ! [kg/m2/s] 
        qla(i,j) = - rhoair * sll_cnst * cen * dqr * wv ! [W/m2]
        qsn(i,j) = - rhoair * cpa * ctn * dtemp * wv ! [W/m2]
!#ifdef OGCM_TAUBULK
        cdt(i,j) = rhoair * cdn * wv ! MKS same as new MRI.COM
!#endif /* OGCM_TAUBULK */
!#else /* OGCM_LYCOEF */
!        cpa = dgill1 * (1.0d0 + dgill2 * qu(i,j))
!        evp(i,j) = rhoair * wv * cen * dqr           ! [kg / m2 / s]
!        qla(i,j) = - rhoair * sll * cen * dqr * wv   ! [W/m2]
!        qsn(i,j) = - rhoair * cpa * ctn * dtemp * wv ! [W/m2]
!#ifdef OGCM_TAUBULK
        cdt(i,j) = rhoair * cdn * wv 
!#endif /* OGCM_TAUBULK */

!#endif /* OGCM_LYCOEF */

      else

        evp(i,j) = 0.0d0
        qla(i,j) = 0.0d0
        qsn(i,j) = 0.0d0

      end if

    end do
  end do
!$omp end parallel
  !
!#ifdef OGCM_TAUBULK ! Always calculate wind stress

  do j = 1, jmx
    do i = 1, imx
      if (atexl(i,j) == 1.0d0) then
        cdn = cdt(i,j)
        wsx(i,j) = cdn * us(i,j) ! [N/m2], us [m/s]
        wsy(i,j) = cdn * vs(i,j) ! [N/m2], vs [m/s]
      else
        wsx(i,j) = 0.0d0
        wsy(i,j) = 0.0d0
      end if
    end do
  end do

!#endif /* OGCM_TAUBULK */

end subroutine bulk
