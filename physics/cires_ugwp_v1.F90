!>  \file cires_ugwp_v1.F90
!! This file contains the Version 1 Unified Gravity Wave Physics (UGWP) scheme by Valery Yudin (University of Colorado, CIRES)
!! See Valery Yudin's presentation at 2017 NGGPS PI meeting:
!! Gravity waves (GWs): Mesoscale GWs transport momentum, energy (heat), and create eddy mixing in the whole atmosphere domain; Breaking and
!dissipating GWs deposit: (a) momentum; (b) heat (energy); and create (c) turbulent mixing of momentum, heat, and tracers
!! To properly incorporate GW effects (a-c) unresolved by DYCOREs we need GW physics
!! "Unified": a) all GW effects due to both dissipation/breaking; b) identical GW solvers for all GW sources; c) ability to replace solvers.
!! Unified Formalism:
!! 1. GW Sources: Stochastic and physics based mechanisms for GW-excitations in the lower atmosphere, calibrated by the high-res analyses/forecasts, and observations (3 types of GW sources: orography, convection, fronts/jets).
!! 2. GW Propagation: Unified solver for "propagation, dissipation and breaking" excited from all type of GW sources.
!! 3. GW Effects: Unified representation of GW impacts on the "resolved" flow for all sources (energy-balanced schemes for momentum, heat and mixing).
!!
!! https://www.weather.gov/media/sti/nggps/Presentations%202017/02%20NGGPS_VYUDIN_2017_.pdf

module cires_ugwp_v1

    use machine, only: kind_phys

    use cires_ugwp_module_v1, only: knob_ugwp_version, cires_ugwp_init_v1, &
                                 cires_ugwp_finalize_v1

    use gwdps, only: gwdps_run

    implicit none

    private

    public cires_ugwp_v1_run

    logical :: is_initialized = .False.


contains



! -----------------------------------------------------------------------
! CCPP entry points for CIRES Unified Gravity Wave Physics (UGWP) scheme v0
! -----------------------------------------------------------------------
!    originally from cires_ugwp_driver_v1.f
!    driver of cires_ugwp   (_driver)
! -----------------------------------------------------------------------
!   driver is called after pbl & before chem-parameterizations
! -----------------------------------------------------------------------
!  order = dry-adj=>conv=mp-aero=>radiation -sfc/land- chem -> vertdiff-> [rf-gws]=> ion-re
! -----------------------------------------------------------------------
!>@brief The subroutine initializes the CIRES UGWP
!> \section arg_table_cires_ugwp_init Argument Table
!! \htmlinclude cires_ugwp_v1_init.html
!!
! -----------------------------------------------------------------------
!


subroutine cires_ugwp_v1_init(Model,Init_parm,Grid)
                ! NOTE:  This code is taken from Valery's code in subroutine
                ! GFS_initialize in GFS_driver.F90

use cires_ugwp_module,   only: cires_ugwp_init_emc, cires_ugwp_init_v1, cires_indx_ugwp

implicit none

real(kind=kind_phys), parameter   :: p_ref = 101325.0d0


!----  initialization of cires_ugwp .
!     if ( Model%me == Model%master) print *,  ' VAY-nml ',  Model%fn_nml
!     if ( Model%me == Model%master) print *,  ' VAY-nml2 ', Model%input_nml_file
    if (Model%cdmbgwd(3) > 0.0) then 
!     if ( Model%me == Model%master) print *,  ' VAY-nml ',  Model%fn_nml,
!     Model%input_nml_file
     call cires_ugwp_init_emc(Model%me,      Model%master, Model%nlunit,  Init_parm%logunit, &
                           Model%fn_nml,  Model%lonr,   Model%latr,    Model%levs,          &    
                           Init_parm%ak,  Init_parm%bk, p_ref,         Model%dtp,           &    
                           Model%cdmbgwd(1:4), Model%cgwf,   Model%prslrd0, Model%ral_ts)
                           
    endif

   if (Model%do_ugwp .and. Model%cdmbgwd(3) == 0.0 ) then 
!

    
     call cires_ugwp_init_v1(Model%me,      Model%master, Model%nlunit,  Init_parm%logunit, &
                Model%Jdat, Model%fn_nml,  Model%lonr,   Model%latr,    Model%levs,         &    
                            Init_parm%ak,  Init_parm%bk, p_ref,         Model%dtp,          &    
                           Model%cdmbgwd(1:2), Model%cgwf,   Model%prslrd0, Model%ral_ts)  
                           
!                          

     do nb = 1, nblks
     
       call cires_indx_ugwp (Init_parm%blksz(nb),Model%me,  Model%master, Grid(nb)%xlat_d,  &
         Grid(nb)%jindx1_tau, Grid(nb)%jindx2_tau, Grid(nb)%ddy_j1tau, Grid(nb)%ddy_j2tau,  &
         Grid(nb)%jindx1_qbo, Grid(nb)%jindx2_qbo, Grid(nb)%ddy_j1qbo, Grid(nb)%ddy_j2qbo,  &
         Grid(nb)%dexp_latqbo)  
          
        if ( Model%me == Model%master) print *,  ' VAY- cires_indx_ugwp NB=',   nb,  Init_parm%blksz(nb)     
      enddo     
    
   endif

end subroutine cires_ugwp_v1_init




! -----------------------------------------------------------------------
!>@brief The subroutine initializes the CIRES UGWP
!> \section arg_table_cires_ugwp_finalize Argument Table
!!
! -----------------------------------------------------------------------
!

  subroutine cires_ugwp_finalize
!
! deallocate sources/spectra & some diagnostics need to find where "deaalocate
! them"
! before "end" of the FV3GFS
!
    implicit none
!
!   deallocate arrays employed in:
!     cires_ugwp_advance / cires_ugwp_driver / cires_ugwp_init
!
    deallocate( kvg,   ktg  )
    deallocate( krad,  kion )
    deallocate( zkm,   pmb  )
    deallocate( rfdis, rfdist)
    deallocate(ugwp_taulat,  ugwp_qbolat)
    deallocate(tau_limb, uzmf_merra)
    deallocate(days_limb,  days_merra,  pmb127)

   end subroutine cires_ugwp_finalize




!>@brief These subroutines and modules execute the CIRES UGWP Version 1
!>\defgroup cires_ugwp_run Unified Gravity Wave Physics General Algorithm
!> @{
!! The physics of NGWs in the UGWP framework (Yudin et al. 2018 \cite yudin_et_al_2018) is represented by four GW-solvers, which is introduced in Lindzen (1981) \cite lindzen_1981, Hines (1997) \cite hines_1997, Alexander and Dunkerton (1999) \cite alexander_and_dunkerton_1999, and Scinocca (2003) \cite scinocca_2003. The major modification of these GW solvers is represented by the addition of the background dissipation of temperature and winds to the saturation criteria for wave breaking. This feature is important in the mesosphere and thermosphere for WAM applications and it considers appropriate scale-dependent dissipation of waves near the model top lid providing the momentum and energy conservation in the vertical column physics (Shaw and Shepherd 2009 \cite shaw_and_shepherd_2009). In the UGWP-v0, the modification of Scinocca (2003) \cite scinocca_2003 scheme for NGWs with non-hydrostatic and rotational effects for GW propagations and background dissipation is represented by the subroutine \ref fv3_ugwp_solv2_v0. In the next release of UGWP, additional GW-solvers will be implemented along with physics-based triggering of waves and stochastic approaches for selection of GW modes characterized by horizontal phase velocities, azimuthal directions and magnitude of the vertical momentum flux (VMF).
!!
!! In UGWP-v0, the specification for the VMF function is adopted from the GEOS-5 global atmosphere model of GMAO NASA/GSFC, as described in Molod et al. (2015) \cite molod_et_al_2015 and employed in the MERRRA-2 reanalysis (Gelaro et al., 2017 \cite gelaro_et_al_2017). The Fortran subroutine \ref slat_geos5_tamp describes the latitudinal shape of VMF-function as displayed in Figure 3 of Molod et al. (2015) \cite molod_et_al_2015. It shows that the enhanced values of VMF in the equatorial region gives opportunity to simulate the QBO-like oscillations in the equatorial zonal winds and lead to more realistic simulations of the equatorial dynamics in GEOS-5 operational and MERRA-2 reanalysis products. For the first vertically extended version of FV3GFS in the stratosphere and mesosphere, this simplified function of VMF allows us to tune the model climate and to evaluate multi-year simulations of FV3GFS with the MERRA-2 and ERA-5 reanalysis products, along with temperature, ozone, and water vapor observations of current satellite missions. After delivery of the UGWP-code, the EMC group developed and tested approach to modulate the zonal mean NGW forcing by 3D-distributions of the total precipitation as a proxy for the excitation of NGWs by convection and the vertically-integrated  (surface - tropopause) Turbulent Kinetic Energy (TKE). The verification scores with updated NGW forcing, as reported elsewhere by EMC researchers, display noticeable improvements in the forecast scores produced by FV3GFS configuration extended into the mesosphere.
!!
!! In UGWP-v1, (add description here....)
!!
!> \section arg_table_cires_ugwp_run Argument Table
!! \htmlinclude cires_ugwp_v1_run.html
!!
!> \section gen_cires_ugwp CIRES UGWP Scheme General Algorithm
!! @{



!       subroutine cires_ugwp_driver_v1 (me,  master,              &
       subroutine cires_ugwp_v1_run (me,  master,              &
         im,  levs, nmtvr, dtp, kdt, imx, do_ugwp, do_tofd,       & 
         cdmbgwd,  jdat, xlat, xlatd, sinlat, coslat, spgrid,     &
         j1_tau, j2_tau, ddy_j1tau, ddy_j2tau,                    &
         j1_qbo, j2_qbo, ddy_j1qbo, ddy_j2qbo, dexpy,             &	 
         ugrs, vgrs, tgrs, qgrs, prsi, prsl, prslk,               &
         phii, phil, del, oro_stat,  sgh30, kpbl,                 &
         dusfcg, dvsfcg, gw_dudt, gw_dvdt, gw_dtdt, gw_kdis,      &
         tau_tofd, tau_mtb, tau_ogw, tau_ngw,                     &
         zmtb, zlwb, zogw, du3dt_mtb,du3dt_ogw, du3dt_tms,        &
         uqbo, ax_qbo, tau_sat, tau_qbo, rdxzb, fhzero, lprnt, ipr)
	 	 
!-----------------------------------------------------------
! Part 1 "revised" LM97 oro-scheme (if do_ugwp=.false.)
! Part 2  non-stationary multi-wave GWs FV3GFS-v0
! Part 3  Dissipative version of UGWP-tendency application
!         (similar to WAM-2017)
!----------------------------------------------------------- 
       use machine,          only : kind_phys
 
       use ugwp_common ,     only : rgrav
       use cires_ugwp_module,only : calendar_ugwp
       
       use cires_ugwp_module,only  : vert_qbo, latqbo, widqbo, taurel, kz1, kz2  
       use cires_ugwp_module, only : ugwp_qbolat
       use cires_ugwp_module, only : uqboe, u2 => uzmf_merra         
       use ugwp_wmsdis_init, only  : tamp_mpa, ilaunch
!   
       implicit none
!input
       real(kind=kind_phys), parameter :: pogw=1.0, pngw=1.0, pked=1.0
       real(kind=kind_phys), parameter :: fw1_tau=1.0
        
       integer, intent(in) :: me,  master, jdat(8)
       integer, intent(in) :: im, levs, kdt, imx, nmtvr, ipr

       real(kind=kind_phys), intent(in) :: dtp, fhzero
       real(kind=kind_phys), intent(in) :: cdmbgwd(4)

       logical             :: do_ugwp, do_tofd, lprnt
       integer, intent(in) :: kpbl(im)
       real(kind=kind_phys), intent(in), dimension(im) :: xlat, xlatd,   &
                                 sgh30, sinlat, coslat, spgrid               ! spgrid = tile-area
     
       real(kind=kind_phys), intent(in) :: oro_stat(im,nmtvr)
       
       real(kind=kind_phys), intent(inout), dimension(im,levs) :: ugrs
       
       real(kind=kind_phys), intent(in), dimension(im,levs) ::   &
                vgrs, tgrs, qgrs, prsi, prsl, prslk, phii, phil, del
		
! Temporary changes (removed intent(in))
       integer, dimension(im) ::j1_tau, j2_tau, j1_qbo, j2_qbo
       real , dimension(im) :: ddy_j1qbo, ddy_j2qbo, ddy_j2tau, ddy_j1tau 
       real , dimension(im) :: dexpy(im)     
!
!out
!
       real(kind=kind_phys), dimension(im,levs) :: gw_dudt, gw_dvdt,    &
                                                   gw_dTdt, gw_kdis

!-----locals	+ diagnostics output

       real(kind=kind_phys), dimension(im,levs)   :: Pdvdt, Pdudt,      &
                             Pdtdt, Pkdis, ed_dudt, ed_dvdt, ed_dTdt

       real(kind=kind_phys), dimension(im)      :: dusfcg, dvsfcg

       real(kind=kind_phys), dimension(im)      :: rdxzb, zmtb,          &
           zlwb, zogw, tau_mtb, tau_ogw, tau_tofd,  tau_ngw
       real(kind=kind_phys), dimension(im,levs) :: du3dt_mtb, du3dt_ogw, &
                                                   du3dt_tms
      
      real(kind=kind_phys), dimension(im,levs)   ::  zmet
      real(kind=kind_phys), dimension(im,levs+1) ::  zmeti
! 
! extra-diag for NGWs and forcing
!      
      real :: tauabs(im,levs), wrms(im,levs),  trms(im,levs)
      
!  
      
      real, dimension(im, levs) ::  uqbo, ax_qbo
      real, dimension(im)       ::  tau_sat,  tau_qbo   
      real                      ::  dexpz(levs)            
      
      real :: wqbo, dforc, sdexpz
                   
      
      integer :: y4, month, day,  ddd_ugwp, curdate, curday
      integer :: hour     
      real    :: hcurdate, hcurday, fhour, fhrday 
      integer :: kdtrest  
      integer :: curday_ugwp 
      integer :: curday_save=20150101       
      logical :: first_qbo=.true.
      real    ::  hcurday_save =20150101.00   
      save first_qbo, curday_save, hcurday_save
!      , stau_sat, stau_qbo, suqbo,  sax_qbo, dexpy, dexpz               
! locals
      
       integer              :: i, j, k, ix, iz
!
       if (me == master .and. kdt < 2) then
         print *
         write(6,*) 'FV3GFS executes ugwp_driver_v1 '
!        write(6,*) 'FV3GFS execute ugwp_driver_v1 nmtvr=', nmtvr
         write(6,*) ' COORDE EXPER pogw = ' , pogw
!         write(6,*) ' COORDE EXPER pgwd = ' , pgwd
!         write(6,*) ' COORDE EXPER pgwd4 = ', pgwd4
         print *
       endif
 
!
! zlwb is not activated
!       
         zlwb(:) = 0.
         zmeti  = phii*rgrav
         zmet   = phil*rgrav
!
! 1) ORO stationary GWs
!    ------------------
       
 
         call gwdps_oro_v1 (im, levs,  imx,   do_tofd,                 &
                      Pdvdt, Pdudt, Pdtdt, Pkdis,                      &
                      ugrs , vgrs, tgrs, qgrs,KPBL, prsi,del,prsl,     &
                      prslk, zmeti, zmet, dtp, kdt, nmtvr, oro_stat,   &  
                      sgh30,   DUSFCg, DVSFCg, xlatd, sinlat, coslat,  &  
                      spgrid,cdmbgwd(1:2), me, master, rdxzb,          &  
                      zmtb, zogw, tau_mtb, tau_ogw, tau_tofd,          &
                      du3dt_mtb, du3dt_ogw, du3dt_tms)
!
          if (me == master .and. kdt < 2) then
           print *
           write(6,*) 'FV3GFS finished GWDPS_ORO_V1 ugwp_driver_v1 '
	   print *, 'xlatd=', xlatd(1:im)
           print *
          endif  

! --------    
! 2) non-stationary GWs with GEOS-5/MERRA GW-forcing
!    ----------------------------------------------
!--------	
! GMAO GEOS-5/MERRA GW-forcing	lat-dep
!--------
         call slat_geos5_tamp(im, tamp_mpa, xlatd, tau_ngw)
	 
	 y4 = jdat(1); month = jdat(2); day = jdat(3) ; hour = jdat(5)
	 
	 fhour = float(hour)+float(jdat(6))/60. + float(jdat(7))/3600.
	 fhour = (kdt-1)*dtp/3600.
	 fhrday  = fhour/24.  - nint(fhour/24.)
	 fhour = fhrday*24.	
	  	 
	 call calendar_ugwp(y4, month, day, ddd_ugwp)
	 curdate = y4*1000 + ddd_ugwp
	 curday = y4*10000 + month*100 + day
         hcurdate = float(curdate) + fhrday
	 hcurday  = float(curday)  + fhrday
!
     if (mod(fhour,fhzero) == 0 .or. first_qbo) then
     
         call tau_limb_advance(me, master, im, levs, ddd_ugwp, curdate,   &
	    j1_tau, j2_tau, ddy_j1tau, ddy_j2tau,  tau_sat,            kdt ) 
	    
	    if (first_qbo) kdtrest = kdt	  
	    first_qbo = .false.
	    curday_save = curday	 
	    hcurday_save= hcurday	    
      endif
      	    
	 tau_ngw = fw1_tau*tau_ngw + tau_sat*(1.-fw1_tau)   
	    	
     goto 111	              
     if (mod(fhour,fhzero) == 0 .or. first_qbo) then
        	 
         call tau_qbo_advance(me, master, im, levs, ddd_ugwp, curdate, &
	    j1_tau, j2_tau, ddy_j1tau, ddy_j2tau,  j1_qbo, j2_qbo,      &
	    ddy_j1qbo, ddy_j2qbo, tau_sat, tau_qbo,  uqbo, ax_qbo, kdt     )  
	   
	 	  
	 if (me == master) then
	 print *, ' curday_save first_qbo ', curday, curday_save, kdt	   
	 print *, ' hcurdays ', hcurdate,  float(hour)/24.
	 print *, jdat(5), jdat(6),  jdat(7), (kdt-1)*dtp/3600., ' calendar '	 	 
!	   print *, ' curday curday_ugwp first_qbo ', hcurday, first_qbo	 
!	   print *, ' vay_tau-limb U' ,  maxval(uqbo), minval(uqbo) 
!	   print *, ' vay_tau-limb TS' , maxval(tau_sat), minval(tau_sat) 
!	   print *, ' vay_tau-limb TQ' , maxval(tau_qbo), minval(tau_qbo) 	   	   
	 endif	


	    if (first_qbo) kdtrest = kdt	  
	    first_qbo = .false.
	    curday_save = curday	 
	    hcurday_save= hcurday	    
      endif
      	  
               
 		
        
	 if (mod(kdt, 720) == 0 .and. me == master ) then	 
	   print *, ' vay_qbo_U' , maxval(uqbo), minval(uqbo) , kdt 		   	      	   
	 endif     
	 
             wqbo = dtp/taurel     
             do k =1, levs 
!	     	sdexpz =  wqbo*vert_qbo(k)  
		sdexpz =  0.25*vert_qbo(k)          
	       do i=1, im
!	        if (dexpy(i) > 0.0) then
		 dforc = 0.25
!	         ugrs(i,k)  = ugrs(i,k)*(1.-dforc) + dforc*uqbo(i,levs+1-k)
!	         tgrs(i,k)  = tgrs(i,k)*(1.-dforc) + dforc*tqbo(i,levs+1-k)		 
!		endif 
	       enddo
	     enddo	    
! 
111      continue


         call cires_ugwp_solv2_v1(im,   levs,  dtp,                     &
                        tgrs, ugrs,  vgrs,   qgrs, prsl, prsi,          &
                        zmet, zmeti,prslk, xlatd, sinlat, coslat,       &
                        gw_dudt, gw_dvdt, gw_dTdt, gw_kdis,             &
                        tauabs, wrms, trms,   tau_ngw, me, master, kdt)

         if (me == master .and. kdt < 2) then
           print *
           write(6,*)'FV3GFS finished fv3_ugwp_solv2_v1  in ugwp_driver_v0 '
           write(6,*) ' non-stationary GWs with GMAO/MERRA GW-forcing '
           print *
         endif

         do k=1,levs
           do i=1,im
             gw_dtdt(i,k) = pngw*gw_dtdt(i,k) + pogw*Pdtdt(i,k)
             gw_dudt(i,k) = pngw*gw_dudt(i,k) + pogw*Pdudt(i,k)   !+(uqbo(i,levs+1-k)-ugrs(i,k))/21600.
             gw_dvdt(i,k) = pngw*gw_dvdt(i,k) + pogw*Pdvdt(i,k)
             gw_kdis(i,k) = pngw*gw_kdis(i,k)                     !+ pogw*Pkdis(i,k)
           enddo
         enddo
       

       
 
       if (pogw == 0.0) then
!        zmtb = 0.;  zogw =0.
         tau_mtb   = 0.0  ; tau_ogw   = 0.0 ;  tau_tofd = 0.0
         du3dt_mtb = 0.0  ; du3dt_ogw = 0.0 ;  du3dt_tms= 0.0
       endif

       return
 
!=============================================================================
! make "ugwp eddy-diffusion" update for gw_dtdt/gw_dudt/gw_dvdt by solving
! vert diffusion equations & update "Statein%tgrs, Statein%ugrs, Statein%vgrs"
!=============================================================================
!
! 3) application of "eddy"-diffusion to "smooth" UGWP-related tendencies
!------------------------------------------------------------------------------
       
       ed_dudt(:,:) = 0.0 ; ed_dvdt(:,:) = 0.0 ; ed_dtdt(:,:) = 0.0
         
       

!       call edmix_ugwp_v1(im,   levs, dtp,                       &
!                        tgrs, ugrs, vgrs, qgrs, del,             &
!                         prsl, prsi, phil, prslk,                &
!                         gw_dudt, gw_dvdt, gw_dTdt, gw_kdis,     &
!                        ed_dudt, ed_dvdt, ed_dTdt,
!                         me, master, kdt )

      do k=1,levs
        do i=1,im
          gw_dtdt(i,k) = gw_dtdt(i,k) + ed_dtdt(i,k)*pked
          gw_dvdt(i,k) = gw_dvdt(i,k) + ed_dvdt(i,k)*pked
          gw_dudt(i,k) = gw_dudt(i,k) + ed_dudt(i,k)*pked
        enddo
      enddo

!      end subroutine cires_ugwp_driver_v1
      end subroutine cires_ugwp_v1_run
!! @}


!
      subroutine tau_limb_advance(me, master, im, levs, ddd, curdate,   &
	    j1_tau, j2_tau, ddy_j1tau, ddy_j2tau,  tau_sat,            kdt )  
	    



      use machine , only : kind_phys
      
      use cires_ugwp_module, only : ntau_d1y, ntau_d2t       
      use cires_ugwp_module, only : ugwp_taulat, days_limb,  tau_limb
      
!      use cires_ugwp_module, only : ugwp_qbolat,  days_merra, pmb127, days_y4md, days_y4ddd
!      use cires_ugwp_module, only :  tau_qbo,  stau_qbo,  uqboe, u2 => uzmf_merra  

      implicit none

      integer, intent(in) ::    me, master, im, levs, ddd, curdate, kdt    
      integer, intent(in), dimension(im) :: j1_tau, j2_tau
      
      real , intent(in),  dimension(im) :: ddy_j1tau, ddy_j2tau  
           
      real, intent(out) ::  tau_sat(im)
      
      integer           :: i, j1, j2, k, it1, it2, iday
      real              :: tem,  tx1, tx2, w1, w2, day2, day1, ddx
      integer           :: yr1, yr2  
!
      integer           ::  iqbo1=1      
!

	 
	 
            it1 = 2
         do iday=1, ntau_d2t
	    if (float(ddd) .lt. days_limb(iday) ) then
	    it2 = iday
	    exit
	    endif
	 enddo
	 it2 = min(it2,ntau_d2t)	 
	 it1 = max(it2-1,1)
	 if (it2 > ntau_d2t ) then
	  print *, ' it1, it2, ntau_d2t ', it1, it2, ntau_d2t
	  stop
	 endif
	 w2 = (float(ddd)-days_limb(it1))/(days_limb(it2)-days_limb(it1))
	 w1 = 1.0-w2
      do i=1, im	 
	 j1 = j1_tau(i)
	 j2 = j2_tau(i)
	 tx1 = tau_limb(j1, it1)*ddy_j1tau(i)+tau_limb(j2, it1)*ddy_j2tau(i)
	 tx2 = tau_limb(j1, it2)*ddy_j1tau(i)+tau_limb(j2, it2)*ddy_j2tau(i)	 
	 tau_sat(i) =  tx1*w1 + w2*tx2 
      enddo
      
         if (me == master ) then	    
	    print*, maxval(tau_limb), minval(tau_limb), ' tau_limb '
	    print*, ntau_d2t
	    print*, days_limb(1) ,  days_limb(ntau_d2t)	, ddd,  ' days-taulimb '
	    print*, 'curdate  ', curdate
	    print*, maxval(tau_sat), minval(tau_sat), ' tau_sat_fv3 '	        	    
	 endif 
      return

      end subroutine tau_limb_advance



      subroutine tau_qbo_advance(me, master, im, levs, ddd, curdate,   &
	    j1_tau, j2_tau, ddy_j1tau, ddy_j2tau,  j1_qbo, j2_qbo,      &
	    ddy_j1qbo, ddy_j2qbo, tau_sat, tau_qbofv3, uqbo, ax_qbo, kdt )  
	    



      use machine , only : kind_phys
      use cires_ugwp_module, only : ntau_d1y, ntau_d2t, nqbo_d1y, nqbo_d2z, nqbo_d3t  
      use cires_ugwp_module, only : ugwp_taulat, days_limb,  tau_limb
      use cires_ugwp_module, only : ugwp_qbolat,  days_merra, pmb127, days_y4md, days_y4ddd
      use cires_ugwp_module, only :  tau_qbo,  stau_qbo,  uqboe, u2 => uzmf_merra  

      implicit none

      integer, intent(in) ::    me, master, im, levs, ddd, curdate, kdt    
      integer, intent(in), dimension(im) :: j1_tau, j2_tau, j1_qbo, j2_qbo
      
      real , intent(in),  dimension(im) :: ddy_j1qbo, ddy_j2qbo
      real , intent(in),  dimension(im) :: ddy_j1tau, ddy_j2tau  
           
      real, intent(out) ::  tau_sat(im), tau_qbofv3(im)
      real, intent(out) ::  uqbo(im, levs), ax_qbo(im, levs)
      integer           :: i, j1, j2, k, it1, it2, iday
      real              :: tem,  tx1, tx2, w1, w2, day2, day1, ddx
      integer           :: yr1, yr2  
!      integer           ::  itau1/1/
      integer           ::  iqbo1=1 
      save  iqbo1      
!
         if (me == master ) then
	    
	    print*, maxval(u2), minval(u2),  ' uqbo ', kdt
	    print*, maxval(uqboe), minval(uqboe),  ' uqboe '
	   	    print*, maxval(tau_limb), minval(tau_limb), ' tau_limb '
	    print*, ntau_d2t, 	nqbo_d3t
	    print*, days_limb(1) ,  days_limb(ntau_d2t)	, ddd,  ' days-taulimb '
	    print*, days_y4ddd(1) , days_y4ddd(nqbo_d3t)	, ddd, ' days-qbo ', curdate	        	    
	 endif 
	 
	 
            it1 = 2
         do iday=1, ntau_d2t
	    if (float(ddd) .lt. days_limb(iday) ) then
	    it2 = iday
	    exit
	    endif
	 enddo
	 it2 = min(it2,ntau_d2t)	 
	 it1 = max(it2-1,1)
	 if (it2 > ntau_d2t ) then
	  print *, ' it1, it2, ntau_d2t ', it1, it2, ntau_d2t
	  stop
	 endif
	 w2 = (float(ddd)-days_limb(it1))/(days_limb(it2)-days_limb(it1))
	 w1 = 1.0-w2
      do i=1, im	 
	 j1 = j1_tau(i)
	 j2 = j2_tau(i)
	 tx1 = tau_limb(j1, it1)*ddy_j1tau(i)+tau_limb(j2, it1)*ddy_j2tau(i)
	 tx2 = tau_limb(j1, it2)*ddy_j1tau(i)+tau_limb(j2, it2)*ddy_j2tau(i)	 
	 tau_sat(i) =  tx1*w1 + w2*tx2 
      enddo
!      
!   
!
            it2=iqbo1
         do iday=iqbo1, nqbo_d3t
	    if (float(curdate) .lt. days_y4ddd(iday) ) then
	    it2 = iday
	    exit
	    endif
	 enddo
	 it2 = min(it2, nqbo_d3t)	 
	 it1 = max(it2-1, 1)
	 if (it2 > nqbo_d3t ) then
	  print *, ' it1, it2, nqbo_d3t ', it1, it2, nqbo_d3t
	  it2 = nqbo_d3t 
	  it1 = it2-1
	 endif	 
	 
	 
	 yr1 = (nint(days_y4ddd(it1)) -ddd)/1000
	 yr2 = (nint(days_y4ddd(it2)) -ddd)/1000
	
	 day1 = days_y4ddd(it1)	 - yr1*1000.
	 day2 = days_y4ddd(it2)	 - yr2*1000.
	 if (yr1 < yr2) then
	   ddx = 365. 
	   if ( mod(yr2,4).ne.0 ) ddx= 366.
	   day2 = day2+ddx
	 endif  
	 
	 w2 = (float(ddd)-day1) !/(day2-day1)
	 w1 = 1.-w2
!	 	 
        do i=1, im	
	 j1 = j1_qbo(i)
	 j2 = j2_qbo(i)
	 
	 tx1 = stau_qbo(it1)*ddy_j1qbo(i)+stau_qbo(it1)*ddy_j2qbo(i)
	 tau_qbofv3(i) =  tx1
	 
!	 tx2 = stau_qbo(it2)*ddy_j1qbo(i)+stau_qbo(it2)*ddy_j2qbo(i)	 
!	 tau_qbofv3(i) =  tx1*w1 + w2*tx2 

	 do k=1, levs
	   tx1 = u2(j1,k, it1)*ddy_j1qbo(i)+ u2(j2,k, it1)*ddy_j2qbo(i)
	   uqbo(i,k) = tx1
!	   tx2 = u2(j1,k, it2)*ddy_j1qbo(i)+ u2(j2,k, it2)*ddy_j2qbo(i)	   
!	   uqbo(i,k) = tx1*w1 + w2*tx2 
	 enddo	 

        enddo
	iqbo1 = it1
        if (me == master .and. kdt < 1500 ) print *,  'vay_advance ', iqbo1, maxval(uqbo), kdt
!       	 
      return
      end subroutine tau_qbo_advance



!>@}
end module cires_ugwp_v1
