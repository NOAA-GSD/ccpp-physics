
module  cires_ugwp_module

!
!   driver is called after pbl & before chem-parameterizations
!    it uses ugwp_common (like phys_cons) and some module-param od solvers/sources init-modules
!....................................................................................
!  order = dry-adj=>conv=mp-aero=>radiation -sfc/land- chem -> vertdiff-> [rf-gws]=> ion-re
!...................................................................................
!
!
    use ugwp_common, only : arad, pi, pi2, hpscale, rhp, rhp2, rh4 
    implicit none
    logical            :: module_is_initialized
!logical               :: do_ugwp   = .false.              ! control => true - ugwp false old gws + rayeleigh friction
    character(len=8)   :: strsolver='pss-1986'
    logical            :: do_physb_gwsrcs = .false.        ! control for physics-based GW-sources
    logical            :: do_rfdamp       = .false.        ! control for Rayleigh friction inside ugwp_driver
    integer, parameter :: idebug_gwrms=1                   ! control for diag computaions pw wind-temp GW-rms and MF fluxs   
    logical, parameter :: do_adjoro = .false.
    real, parameter    ::  max_kdis = 250.              ! 400 m2/s
    real, parameter    ::  max_axyz = 250.e-5           ! 400 m/s/day
    real, parameter    ::  max_eps =  max_kdis*4.e-7    ! ~16   K/day  max_kdis*BN2/cp 
    real, parameter    ::  maxdudt =  max_axyz
    real, parameter    ::  maxdtdt =  max_eps    
    real, parameter    ::  dked_min = 0.01
    real, parameter    ::  dked_max = max_kdis  
  
 
    real, parameter    :: hps   = hpscale
    real, parameter    :: hpskm = hps/1000.
!

    real, parameter    :: ricrit = 0.25
    real, parameter    :: frcrit = 0.50
    real, parameter    :: linsat = 1.00
    real, parameter    :: linsat2 = linsat*linsat
!
!    integer               ::    curday_ugwp                !  yyyymmdd    20150101 
!    integer               ::    ddd_ugwp                   !  ddd of year from 1-366
    
    integer               :: knob_ugwp_solver=1            ! 1, 2, 3, 4 - (linsat, ifs_2010, ad_gfdl, dsp_dis)
    integer, dimension(4) :: knob_ugwp_source              ! [1,1,1,0]  - (oro, fronts, conv, imbf-owp]
    integer, dimension(4) :: knob_ugwp_wvspec              !  number of waves for- (oro, fronts, conv, imbf-owp]
    integer, dimension(4) :: knob_ugwp_azdir               !   number of wave azimuths for- (oro, fronts, conv, imbf-owp]
    integer, dimension(4) :: knob_ugwp_stoch               !  1 - deterministic ; 0 - stochastic
    real,    dimension(4) :: knob_ugwp_effac               !  efficiency factors for- (oro, fronts, conv, imbf-owp]

    integer               :: knob_ugwp_doaxyz=1            ! 1 -gwdrag
    integer               :: knob_ugwp_doheat=1            ! 1 -gwheat
    integer               :: knob_ugwp_dokdis=0            ! 1 -gwmixing
    integer               :: knob_ugwp_ndx4lh = 2          ! n-number  of  "unresolved" "n*dx" for lh_gw
    integer               :: knob_ugwp_nslope = 1          ! spectral"growth" S-slope of GW-energy spectra mkz^S    
 
    real                  :: knob_ugwp_palaunch = 500.e2   ! fixed pressure layer in Pa for "launch" of NGWs
    real                  :: knob_ugwp_lzmax = 12.5e3      ! 12.5 km max-VERT-WL of GW-spectra
    real                  :: knob_ugwp_lzstar = 2.0e3      ! UTLS  mstar = 6.28/lzstar  2-2.5 km 
    real                  :: knob_ugwp_lzmin = 1.5e3       ! 1.5 km  min-VERT-WL of GW-spectra       
    real                  :: knob_ugwp_taumin = 0.25e-3
    real                  :: knob_ugwp_tauamp = 7.75e-3    ! range from 30.e-3 to 3.e-3 ( space-borne values)
    real                  :: knob_ugwp_lhmet  = 200.e3     ! 200 km
    real                  :: knob_ugwp_tlimb  = .true.    
!
    real                  :: kxw = pi2/200.e3              ! single horizontal wavenumber of ugwp schemes
!
! tune-ups for qbo
!	    
    real                  :: knob_ugwp_qbolev = 500.e2   ! fixed pressure layer in Pa for "launch" of conv-GWs
    real                  :: knob_ugwp_qbosin = 1.86     ! semiannual cycle of tau_qbo_src in radians 
    real                  :: knob_ugwp_qbotav = 2.285e-3 ! additional to "climate" for QBO-sg forcing 
    real                  :: knob_ugwp_qboamp = 1.191e-3 ! additional to "climate" QBO 
    real                  :: knob_ugwp_qbotau = 10.      !  relaxation time scale in days
    real                  :: knob_ugwp_qbolat = 15.      !  qbo-domain for extra-forcing
    real                  :: knob_ugwp_qbowid = 7.5      !  qbo-attenuation for extra-forcing
    character(len=8)      :: knob_ugwp_orosolv='pss-1986'

    character(len=255)    :: ugwp_qbofile =  'qbo_zmf_2009_2018.nc'
    character(len=255)    :: ugwp_taufile =  'ugwp_limb_tau.nc' 
    
!    character(len=250)     :: knob_ugwp_qbofile='qbo_zmf_2009_2018.nc'!
!    character(len=250)     :: knob_ugwp_amffile='mern_zmf_amf_12month.nc'
!     character(len=255)    :: file_limb_tab='ugwp_limb_tau.nc'

!     integer, parameter    :: ny_tab=73, nt_tab=14
!     real,    parameter    :: rdy_tab = 1./2.5,  rdd_tab = 1./30. 
!     real                  :: days_tab(nt_tab), lat_tab(ny_tab)
!     real                  :: abmf_tab(ny_tab,nt_tab)   

    integer  :: ugwp_azdir
    integer  :: ugwp_stoch

    integer  :: ugwp_src
    integer  :: ugwp_nws
    real     :: ugwp_effac

!
    data knob_ugwp_source / 1,0, 1, 0 /                    !  oro-conv-fjet-okw-taub_lat:      1-active 0-off
    data knob_ugwp_wvspec /1,32,32,32/                     !  number of waves for- (oro, fronts, conv, imbf-owp, taulat]
    data knob_ugwp_azdir  /2, 4, 4,4/                      !  number of wave azimuths for- (oro, fronts, conv, imbf-okwp]
    data knob_ugwp_stoch  /0, 0, 0,0/                      !  0 - deterministic ; 1 - stochastic, non-activated option
    data knob_ugwp_effac  /1.,1.,1.,1./                    !  efficiency factors for- (oro, fronts, conv, imbf-owp]
    integer  :: knob_ugwp_version = 0
    integer  :: launch_level = 55
!
    namelist /cires_ugwp_nml/ knob_ugwp_solver, knob_ugwp_source,knob_ugwp_wvspec, knob_ugwp_azdir,  &
            knob_ugwp_stoch,  knob_ugwp_effac,knob_ugwp_doaxyz,  knob_ugwp_doheat, knob_ugwp_dokdis, &
            knob_ugwp_ndx4lh, knob_ugwp_version, knob_ugwp_palaunch, knob_ugwp_nslope,  knob_ugwp_lzmax, &
	    knob_ugwp_lzmin,  knob_ugwp_lzstar, knob_ugwp_lhmet, knob_ugwp_tauamp, knob_ugwp_taumin,  &
	    knob_ugwp_qbolev, knob_ugwp_qbosin, knob_ugwp_qbotav, knob_ugwp_qboamp, knob_ugwp_qbotau, &
	    knob_ugwp_qbolat, knob_ugwp_qbowid, knob_ugwp_tlimb, knob_ugwp_orosolv  

!&cires_ugwp_nml
! knob_ugwp_solver=2
! knob_ugwp_source=1,1,1,0
! knob_ugwp_wvspec=1,32,32,32
! knob_ugwp_azdir =2, 4, 4,4
! knob_ugwp_stoch =0, 0, 0,0
! knob_ugwp_effac=1, 1, 1,1
! knob_ugwp_doaxyz=1
! knob_ugwp_doheat=1
! knob_ugwp_dokdis=0
! knob_ugwp_ndx4lh=4
!/
!
! allocatable arrays, initilized during "cires_ugwp_init" &
!                     released   during "cires_ugwp_finalize"
!
   real, allocatable :: kvg(:), ktg(:), krad(:), kion(:)
   real, allocatable :: zkm(:), pmb(:)
   real, allocatable :: rfdis(:), rfdist(:)
   integer           :: levs_rf
   real              :: pa_rf, tau_rf
!
! tabulated GW-sources
!
   integer           :: ntau_d1y, ntau_d2t, nqbo_d1y, nqbo_d2z, nqbo_d3t  
   real, allocatable :: ugwp_taulat(:),  ugwp_qbolat(:)
   real, allocatable :: tau_limb(:,:), days_limb(:)
   real, allocatable :: uzmf_merra(:,:,:),  days_merra(:),  pmb127(:)
   real, allocatable :: uqboe(:,:)  
   real, allocatable ::  days_y4ddd(:),  zkm127(:)
   real, allocatable ::  tau_qbo(:),   stau_qbo(:) 
   integer,allocatable ::  days_y4md(:) 
   real, allocatable   ::   vert_qbo(:)
   
!
! limiters
!
   real, parameter    :: latqbo =20., widqbo=15., taurel = 21600.        
   integer, parameter :: kz2 = 127-7, kz1= 127-49, kz5=5    !  64km - 18km   
!   

!======================================================================
   real, parameter :: F_coriol=1                    ! Coriolis effects
   real, parameter :: F_nonhyd=1                    ! Nonhydrostatic waves
   real, parameter :: F_kds   =0                    ! Eddy mixing due to GW-unstable below
   real, parameter :: iPr_ktgw =1./3., iPr_spgw=iPr_ktgw 
   real, parameter :: iPr_turb =1./3., iPr_mol =1.95
   real, parameter :: rhp1=1./hps, rh2=0.5*rhp1, rhp4 = rh2*rh2
   real, parameter :: khp =  0.287*rhp1             ! R/Cp/Hp
   real, parameter :: cd_ulim = 1.0                 ! critical level precision or Lz ~ 0 ~dz of model

   contains
!
! -----------------------------------------------------------------------
!
! init  of cires_ugwp   (_init)  called from GFS_driver.F90
!
! -----------------------------------------------------------------------



  subroutine cires_ugwp_init_v1 (me, master, nlunit, logunit, jdat_gfs, fn_nml2, &
              lonr, latr, levs, ak, bk, pref, dtp, cdmvgwd, cgwf,                &
              pa_rf_in, tau_rf_in)
!
!  input_nml_file ='input.nml'=fn_nml   ..... OLD_namelist and cdmvgwd(4) Corrected Bug Oct 4
!
    use  netcdf
    use  ugwp_oro_init,     only :  init_oro_gws
    use  ugwp_conv_init,    only :  init_conv_gws
    use  ugwp_fjet_init,    only :  init_fjet_gws
    use  ugwp_okw_init,     only :  init_okw_gws
    use  ugwp_wmsdis_init,  only :  initsolv_wmsdis
    
    use  ugwp_lsatdis_init, only :  initsolv_lsatdis
       
    
    use  ugwp_wmsdis_init,  only : ilaunch, nslope, lhmet, lzmax, lzmin, lzstar   
    use  ugwp_wmsdis_init,  only : tau_min, tamp_mpa    
    implicit none

    integer, intent (in) :: me
    integer, intent (in) :: master
    integer, intent (in) :: nlunit
    integer, intent (in) :: logunit
    integer, intent (in) :: lonr
    integer, intent (in) :: levs
    integer, intent (in) :: latr
    integer, intent (in) :: jdat_gfs(8)
    real,    intent (in) :: ak(levs+1), bk(levs+1), pref
    real,    intent (in) :: dtp
    real,    intent (in) :: cdmvgwd(2), cgwf(2)             ! "scaling" controls for "old" GFS-GW  dims(2) !!!
    real,    intent (in) :: pa_rf_in, tau_rf_in
 
    character(len=64), intent (in) :: fn_nml2
    character(len=64), parameter   :: fn_nml='input.nml'

!    character,  intent (in) :: input_nml_file
!    integer, parameter :: logunit =  6
    integer :: ios
    logical :: exists
    real    :: dxsg
    
    integer :: ncid,  iernc, vid, dimid, status         
    integer :: k
    integer :: ddd_ugwp,    curday_ugwp 
    real    :: avqbo(6)
    avqbo = [0.05, 0.1, 0.25, 0.5, 0.75, 0.95]
!
    if (me == master) print *, trim (fn_nml), ' GW-namelist file '
    inquire (file =trim (fn_nml) , exist = exists)
!
    if (.not. exists) then
       if (me == master) &
        write (6, *) 'separate ugwp :: namelist file: ', trim (fn_nml), ' does not exist'
    else
        open (unit = nlunit, file = trim(fn_nml), action = 'read', status = 'old', iostat = ios)
    endif
    rewind (nlunit)
    read   (nlunit, nml = cires_ugwp_nml)
    close  (nlunit)
!
    strsolver= knob_ugwp_orosolv     
    pa_rf   = pa_rf_in
    tau_rf  = tau_rf_in
    
    curday_ugwp = jdat_gfs(1)*10000 + jdat_gfs(2)*100 +jdat_gfs(3)    
    call calendar_ugwp(jdat_gfs(1), jdat_gfs(2), jdat_gfs(3), ddd_ugwp)
        
! write version number and namelist to log file
    if (me == master) then
        write (logunit, *) " ================================================================== "
        write (logunit, *) "cires_ugwp_namelist_extended_v1"
        write (logunit, nml = cires_ugwp_nml)
        write (logunit, *) " ================================================================== "
	
        write (6, *) " ================================================================== "
        write (6, *) "cires_ugwp_namelist_extended_v1"
        write (6, nml = cires_ugwp_nml)		
        write (6, *) " ================================================================== "	
        write (6, *) "calendar_ugwp ddd_ugwp=", ddd_ugwp
        write (6, *) "calendar_ugwp curday_ugwp=", curday_ugwp
        write (6, *) " ================================================================== "	
        write (6, *) 	 ddd_ugwp, ' jdat_gfs ddd of year '	
    endif
!
! effective kxw - resolution-aware
!
    dxsg =  pi2*arad/float(lonr) * knob_ugwp_ndx4lh
    kxw  = pi2/knob_ugwp_lhmet
!
!   kxw  =  pi2/dxsg
!
! init global background dissipation for ugwp -> 4d-variable for fv3wam linked with pbl-vert_diff
!

!   allocate(fcor(latr), fcor2(latr)  )
!
    allocate( kvg(levs+1),   ktg(levs+1)  )
    allocate( krad(levs+1),  kion(levs+1) )        
    allocate( zkm(levs),   pmb(levs) )
    allocate( rfdis(levs), rfdist(levs) )
    
    allocate (vert_qbo(levs))    
    
!
! ak -pa  bk-dimensionless  from surf => tol_lid_pressure =0
!
    
    do k=1, levs
       pmb(k) = 1.e0*(ak(k) + pref*bk(k))    ! Pa -unit  Pref = 1.e5, pmb = Pa
       zkm(k) = -hpskm*alog(pmb(k)/pref)
    enddo   
    vert_qbo(1:levs) = 0.   
     
   do k=kz1, kz2        		  
      vert_qbo(k)=1.
      if (k.le.(kz1+kz5))  vert_qbo(k)  = avqbo(k+1-kz1)
      if (k.ge.(kz2-kz5))  vert_qbo(k) =  avqbo(kz2+1-k)
      if (me == master) print *, 'vertqbo', vert_qbo(k), zkm(k) 		   
    enddo
    
!
! find ilaunch
!
    
   do k=levs, 1, -1
     if (pmb(k) .gt. knob_ugwp_palaunch ) exit
   enddo

      launch_level = max(k-1, 5)    ! above 5-layers from the surface
    
!
! Part-1 :init_global_gwdis
!
    call init_global_gwdis(levs, zkm, pmb, kvg, ktg, krad, kion, pa_rf, tau_rf,  me,  master)
    call rf_damp_init     (levs, pa_rf, tau_rf, dtp, pmb, rfdis, rfdist, levs_rf)
!
! Part-2 :init_SOURCES_gws
!
    
!    
! call init-solver for "stationary" multi-wave spectra and sub-grid oro
!
    call init_oro_gws( knob_ugwp_wvspec(1), knob_ugwp_azdir(1), &
         knob_ugwp_stoch(1), knob_ugwp_effac(1), lonr, kxw, cdmvgwd )
!
! call init-sources for "non-sationary" multi-wave spectra
!
    do_physb_gwsrcs=.true.

    IF (do_physb_gwsrcs) THEN

      if (me == master) print *, ' do_physb_gwsrcs ',  do_physb_gwsrcs, ' in cires_ugwp_init '
      if (knob_ugwp_wvspec(4) > 0) then
! okw
        call init_okw_gws(knob_ugwp_wvspec(4), knob_ugwp_azdir(4), &
                          knob_ugwp_stoch(4), knob_ugwp_effac(4), lonr, kxw )
        if (me == master) print *, ' init_okw_gws '
      endif

      if (knob_ugwp_wvspec(3) > 0) then
! fronts
        call init_fjet_gws(knob_ugwp_wvspec(3), knob_ugwp_azdir(3), &
                           knob_ugwp_stoch(3), knob_ugwp_effac(3), lonr, kxw )
        if (me == master) print *, ' init_fjet_gws '
      endif

      if (knob_ugwp_wvspec(2) > 0) then
! conv
        call init_conv_gws(knob_ugwp_wvspec(2), knob_ugwp_azdir(2), &
                           knob_ugwp_stoch(2), knob_ugwp_effac(2), lonr, kxw, cgwf )
        if (me == master)   &
           print *, ' init_convective GWs cgwf', knob_ugwp_wvspec(2), knob_ugwp_azdir(2)

      endif

     ENDIF   !IF (do_physb_gwsrcs)
!
!
! Tabulated sources 
!
!      goto 121
      
      iernc=NF90_OPEN(trim(ugwp_taufile), nf90_nowrite, ncid)
     
       if(iernc.ne.0) then         
          write(6,*) 
          write(6,*) ' cannot open file_limb_tab data-file',  trim(ugwp_taufile)  
          write(6,*) 
          stop	  
        else


       status = nf90_inq_dimid(ncid, "lat", DimID)
!      if (status /= nf90_noerr) call handle_err(status)
!
       status = nf90_inquire_dimension(ncid, DimID,  len =ntau_d1y )
       
       status = nf90_inq_dimid(ncid, "days", DimID)
       status = nf90_inquire_dimension(ncid, DimID,  len =ntau_d2t )
           if (me == master)  print *, ntau_d1y, ntau_d2t, ' dimd-tlimb '
        allocate (ugwp_taulat(ntau_d1y ), days_limb(ntau_d2t))
        allocate ( tau_limb (ntau_d1y, ntau_d2t ))         
                  
	iernc=nf90_inq_varid( ncid, 'DAYS', vid )
        iernc= nf90_get_var( ncid, vid, days_limb)
	iernc=nf90_inq_varid( ncid, 'LATS', vid )
        iernc= nf90_get_var( ncid, vid, ugwp_taulat)
	iernc=nf90_inq_varid( ncid, 'ABSMF', vid )
        iernc= nf90_get_var( ncid, vid, tau_limb)
			
	iernc=nf90_close(ncid)
	
	endif
!	
       iernc=NF90_OPEN(trim(ugwp_qbofile), nf90_nowrite, ncid)
     
       if(iernc.ne.0) then         
          write(6,*) 
          write(6,*) ' cannot open qbofile data-file',  trim(ugwp_qbofile)  
          write(6,*) 
          stop	  
        else
	
       status = nf90_inq_dimid(ncid, "lat", DimID)
       status = nf90_inquire_dimension(ncid, DimID, len =nqbo_d1y )
       status = nf90_inq_dimid(ncid, "lev", DimID)
       status = nf90_inquire_dimension(ncid, DimID, len =nqbo_d2z)       
       status = nf90_inq_dimid(ncid, "days", DimID)
       status = nf90_inquire_dimension(ncid, DimID, len =nqbo_d3t )
       if (me == master) print *, nqbo_d1y, nqbo_d2z, nqbo_d3t, ' dims tauqbo '
        allocate (ugwp_qbolat(nqbo_d1y ), days_merra(nqbo_d3t) )
        allocate (zkm127(nqbo_d2z), pmb127(nqbo_d2z))
        allocate ( uzmf_merra (nqbo_d1y, nqbo_d2z, nqbo_d3t )) 
        allocate ( uqboe (nqbo_d2z, nqbo_d3t )) 	        
        allocate (days_y4ddd(nqbo_d3t), days_y4md(nqbo_d3t) )
	allocate (tau_qbo(nqbo_d3t), stau_qbo(nqbo_d3t) )  
	               
	iernc=nf90_inq_varid( ncid, 'DAYS', vid )	
        iernc= nf90_get_var( ncid, vid, days_merra)
	
	iernc=nf90_inq_varid( ncid, 'Y4MD', vid )	
        iernc= nf90_get_var( ncid, vid, days_y4md)
	
	iernc=nf90_inq_varid( ncid, 'Y4DDD', vid )	
        iernc= nf90_get_var( ncid, vid, days_y4ddd)		

	iernc=nf90_inq_varid( ncid, 'LATS', vid )
        iernc= nf90_get_var( ncid, vid, ugwp_qbolat)

	iernc=nf90_inq_varid( ncid, 'LEVS', vid )
        iernc= nf90_get_var( ncid, vid, zkm127)	

	
	iernc=nf90_inq_varid( ncid, 'UQBO', vid )
        iernc= nf90_get_var( ncid, vid, uzmf_merra)
	
	iernc=nf90_inq_varid( ncid, 'TAUQBO', vid )
        iernc= nf90_get_var( ncid, vid, tau_qbo)	
	
	iernc=nf90_inq_varid( ncid, 'STAUQBO', vid )
        iernc= nf90_get_var( ncid, vid, stau_qbo)
	iernc=nf90_inq_varid( ncid, 'UQBOE', vid )
        iernc= nf90_get_var( ncid, vid, uqboe)				
	iernc=nf90_close(ncid)	
      endif	
      
     if (me == master)	then
        print *
        print *, ' ugwp_tabulated files input '
	print *, ' ugwp_taulat ', ugwp_taulat
	print *, ' days ', days_limb
        print *, ' TAU-limb ', maxval(tau_limb)*1.e3, minval(tau_limb)*1.e3
        print *, ' TAU-qbo ', maxval(stau_qbo)*1.e3, minval(stau_qbo)*1.e3	
	print *, ' YMD-qbo ', maxval(days_y4md), minval(days_y4md)
	print *, ' YDDD-qbo ', maxval(days_y4ddd), minval(days_y4ddd)
	print *, ' uzmf_merra ',maxval(uzmf_merra), minval(uzmf_merra)
	print *, ' uEq_merra ',maxval(uqboe), minval(uqboe)
	print *	
     endif
     
!
121  continue
!    endif   ! tabulated sources SABER/HIRDLS/QBO
    
!======================
! Part-3 :init_SOLVERS
! =====================
!
! call init-solvers for "broad" non-stationary multi-wave spectra
!
    if   (knob_ugwp_solver==1) then
!
      call initsolv_lsatdis(me, master, knob_ugwp_wvspec(2), knob_ugwp_azdir(2), &
                            knob_ugwp_stoch(2), knob_ugwp_effac(2), do_physb_gwsrcs, kxw )
    endif
     if   (knob_ugwp_solver==2) then 
!
! re-assign from namelists 
!     
       nslope = knob_ugwp_nslope             ! the GW sprctral slope at small-m             
       lzstar = knob_ugwp_lzstar
       lzmax  = knob_ugwp_lzmax 
       lzmin  = knob_ugwp_lzmin      
       lhmet  = knob_ugwp_lhmet
       tamp_mpa =knob_ugwp_tauamp           !amplitude for GEOS-5/MERRA-2
       tau_min  =knob_ugwp_taumin             ! min of GW MF 0.25 mPa                        
       ilaunch = launch_level
       kxw = pi2/lhmet
       call initsolv_wmsdis(me, master, knob_ugwp_wvspec(2), knob_ugwp_azdir(2), &
                   knob_ugwp_stoch(2), knob_ugwp_effac(2), do_physb_gwsrcs, kxw)
     endif
!
! other solvers not yet tested for fv3gfs
!
!<    if   (knob_ugwp_solver==3) call init_dspdis
!<    if   (knob_ugwp_solver==4) call init_adodis
! 

!======================
    module_is_initialized = .true.
    if (me == master)   print *, ' CIRES-ugwp-V1 is initialized ', module_is_initialized

    end subroutine cires_ugwp_init_v1
! -----------------------------------------------------------------------
!
!    driver  of cires_ugwp   (_driver)
!    called from GFS_physics_driver.F90
!
! -----------------------------------------------------------------------
!      call cires_ugwp_driver                                        &
!       (im, levs, dtp, kdt, me, lprnt,  Model%lonr,                 &
!       Model%prslrd0, Model%ral_ts,  Model%cdmbgwd,                 &
!       Grid%xlat, Grid%xlat_d, Grid%sinlat,  Grid%coslat,           &
!       Statein, delp_gws,  Oro_stat,                                & 
!       dusfcg, dvsfcg, gw_dudt,  gw_dvdt, gw_dtdt, gw_kdis,         &       
!       Diag%gwp_ax, Diag%gwp_axo, Diag%gwp_axc, Diag%gwp_axf,       &
!       Diag%gwp_ay, Diag%gwp_ayo, Diag%gwp_ayc, Diag%gwp_ayf,       &
!       Diag%gwp_dtdt,   Diag%gwp_kdis, Diag%gwp_okw, Diag%gwp_fgf,  &
!       Diag%gwp_dcheat, Diag%gwp_precip, Diag%gwp_klevs,            &
!       Diag%zmtb,   Diag%gwp_scheat, dlength, cldf,                 & 
!       Diag%tau_tofd, Diag%tau_mtb, Diag%tau_ogw, Diag%tau_ngw,     &
!       Diag%zmtb, Diag%zlwb, Diag%zogw, Diag%du3dt_mtb,             &
!       Diag%du3dt_ogw, Diag%du3dt_tms )

    subroutine cires_ugwp_driver                                     &
       (im, levs, dtp, kdt, me, lprnt,   lonr,                       & 
       pa_rf, tau_rf, cdmbgwd,  xlat, xlatd, sinlat,  coslat,        &
       ugrs, vgrs, tgrs, qgrs, prsi, prsl, prslk, phii, phil,        &
       delp, orostat, kpbl,                                          & 
       dusfc, dvsfc, dudt,  dvdt, dtdt, kdis,                        &
       axtot, axo, axc, axf,  aytot, ayo, ayc, ayf,                  &
       eps_tot,  ekdis,  trig_okw, trig_fgf,                         &
       dcheat, precip, cld_klevs, zmtb, scheat, dlength, cldf,       & 
       taus_sso, taus_ogw, tauf_ogw, tauf_ngw,                       &
       ugw_zmtb, ugw_zlwb, ugw_zogw, ugw_axmtb, ugw_axlwb, ugw_axtms ) 

!
     use machine,       only: kind_phys
     use physcons,      only: con_cp, con_fvirt, con_g, con_rd
     use ugwp_common,   only: omega2 
!
!
     use ugwp_okw_init,  only : &
      eff_okw, nstokw, nwokw, ch_okwp, nazokw, spf_okwp, xaz_okwp, yaz_okwp
     use ugwp_conv_init, only : &
      eff_con, nstcon, nwcon, ch_conv, nazcon, spf_conv, xaz_conv, yaz_conv
     use ugwp_fjet_init, only : &
      eff_fj,  nstfj,  nwfj,  ch_fjet, nazfj,  spf_fjet, xaz_fjet, yaz_fjet         

!     
     implicit none
!
 
     logical         :: lprnt
     integer         :: me, im, levs, kdt, lonr
     real(kind_phys) :: dtp
     real(kind_phys) :: pa_rf, tau_rf
     real(kind_phys) :: cdmbgwd(2)

     integer,         intent(in)        ::  kpbl(im)
     real(kind_phys)                    ::  hpbl(im)
     real(kind_phys), intent(in)        ::  orostat(im, 14)
     real(kind_phys), intent(in), dimension(im,levs) :: ugrs, vgrs, &
                   tgrs, qgrs, prsi, prsl, prslk, phii, phil, delp
!
     real(kind_phys), dimension(im)       :: xlat, xlatd, sinlat, coslat
     real(kind_phys), dimension(im, levs) :: trig_okw, trig_fgf
     real(kind_phys), dimension(im)       :: precip              ! precip-n rates and
     integer        , dimension(im, 3)    :: cld_klevs           ! indices fo cloud top/bot/?
     real(kind_phys), dimension(im, levs) :: dcheat,  scheat     ! deep and shal conv heat tend.


     real(kind_phys), dimension(im)       :: dlength           ! tail-grid  box scale in meters
     real(kind_phys), dimension(im)       :: cldf              ! "bizzard" old cgwd-tuning knobs dimensionless
!===================
! tendency + kdis
!===================
     real(kind_phys), dimension(im, levs) :: dudt,  dvdt, dtdt, kdis
     real(kind_phys), dimension(im, levs) :: axtot, axo, axc, axf
     real(kind_phys), dimension(im, levs) :: aytot, ayo, ayc, ayf
     real(kind_phys), dimension(im, levs) :: eps_tot,  ekdis
 
!
     real(kind_phys), dimension(im, levs) :: eds_o, kdis_o
     real(kind_phys), dimension(im, levs) :: eds_c, kdis_c
     real(kind_phys), dimension(im, levs) :: eds_f, kdis_f
     real(kind_phys), dimension(im, levs) :: ax_rf, ay_rf, eps_rf
!
!==================================================================================
! diagnostics for OGW & NGW  +  SSO effects axmtb, axlwb, axtms
!==================================================================================
     real(kind_phys), dimension(im)       ::  dusfc, dvsfc
     real(kind_phys), dimension(im)       ::  taus_sso, taus_ogw, tauf_ogw, tauf_ngw
     real(kind_phys), dimension(im)       ::  ugw_zmtb, ugw_zlwb, ugw_zogw
     real(kind_phys), dimension(im, levs) ::  ugw_axmtb,ugw_axlwb, ugw_axtms
     real(kind_phys), dimension(im, levs) ::  tauz_ogw, tauz_ngw,  wtauz

!
!    knob_ugwp_source=[ 1,    1,    1,       0 ]
!                       oro conv    nst   imbal-okw
!    locals
!
    integer :: i, j, k, istype, ido
!
! internal diagnostics for oro-waves, lee waves, and mtb :
!
     real(kind_phys), dimension(im) :: dusfc_mb, dvsfc_mb, dusfc_ogw, dvsfc_ogw
     real(kind_phys), dimension(im) :: dusfc_lwb, dvsfc_lwb
     real(kind_phys), dimension(im) :: zmtb, zlwb, zogw      ! GW-launch levels in "meters"
!
     real(kind_phys), dimension(im) ::  fcor,  c2f2
!
! three sources with different: a) spectra-content/azimuth; b) efficiency ;c) spectral shape
!
     real(kind_phys), dimension(im) ::  taub_con,  taub_fj, taub_okw
     integer        , dimension(im) ::  klev_okw,  klev_fj, klev_con
     integer        , dimension(im) ::  if_okw,    if_con,  if_fj
     integer                        ::  nf_okw,    nf_con,  nf_fj
!
     dudt  = 0.
     dvdt  = 0.
     dtdt  = 0.
     kdis  = 0.
     axo   = 0. ; axc    = 0. ;  axf   = 0.
     ayo   = 0. ; ayc    = 0. ;  ayf   = 0.
     eds_o = 0. ; kdis_o = 0. ; eds_f  = 0.  ; kdis_f = 0. ;  eds_c = 0. ; kdis_c = 0.
     ax_rf = 0. ; ay_rf  = 0. ; eps_rf = 0
 
     hpbl(:) = 2000.    ! hpbl (1:im) = phil(1:im, kpbl(1:im))
!

     do  i=1, im
       fcor(i) = omega2*sinlat(i)
       c2f2(i) = fcor(i)*fcor(i)/(kxw*kxw)
     enddo

!        i=im
!        print *, i, fcor(i), 6.28e-3/kxw, sqrt(c2f2(i))
!     print *, maxval(statein%prsl/statein%tgrs)/287. , ' density '
 
!
!
! What can be computed for ALL types of GWs? =>
!  "Br-Vi frequency"with "limits" in case of "conv-unstable" layers
!   Background dissipation:  Molecular + Eddy
!   Wind projections may differ from GW-sources/propagation azimuths
!
     do istype=1, size(knob_ugwp_source)

       ido        = knob_ugwp_source(istype)                    ! 0 or 1 off or active

       ugwp_azdir = knob_ugwp_azdir(istype)
       ugwp_stoch = knob_ugwp_stoch(istype)
       ugwp_nws   = knob_ugwp_wvspec(istype)
       ugwp_effac = knob_ugwp_effac(istype)
 
!
! oro-gw effects
!
       if (ido == 1 .and. istype ==1 ) then
!
! 1. solve for OGW effects on the mean flow
! 2. all parts of ORO effexra inside: MTB TOFD LeeWB OGW-drag
!
      call ugwp_oro(im, levs, dtp, kdt, me, lprnt,              &
       fcor, c2f2, ugrs, vgrs, tgrs,                            &
       qgrs, prsi, delp, prsl, prslk, phii, phil,               &
       orostat,  hpbl, axo, ayo, eds_o, kdis_o,                 &
       dusfc, dvsfc, dusfc_mb, dvsfc_mb, dusfc_ogw, dvsfc_ogw,  &
       dusfc_lwb, dvsfc_lwb, zmtb, zlwb, zogw,tauf_ogw,tauz_ogw,&
       ugw_axmtb,ugw_axlwb, ugw_axtms)
!
!       taus_sso, taus_ogw, tauz_ogw, tauz_ngw, tauf_ogw, tauf_ngw,   &
!       ugw_zmtb, ugw_zlwb, ugw_zogw, ugw_axmtb,ugw_axlwb, ugw_axtms 
! collect column-integrated  "dusfc, dvsfc" only for oro-waves
! 
        taus_sso =  dusfc_mb + dusfc_lwb + dusfc_ogw
        taus_ogw =  dusfc_ogw
        ugw_zmtb =  zmtb
        ugw_zlwb =  zlwb
        ugw_zogw =  zogw

!        tauz_ogw/tauf_ogw => output
! ugwp_azdir, ugwp_stoch, ugwp_nws ..... "multi-wave + stochastic"
!
! stationary gw-mode ch=0, with "gw_solver_linsat"
! compute column-integrated  "dusfc, dvsfc" only for oro-waves
!
          dudt = dudt + axo   * ugwp_effac
          dvdt = dvdt + ayo   * ugwp_effac
          dtdt = dtdt + eds_o * ugwp_effac
          kdis = kdis + kdis_o* ugwp_effac   
!	  print *,  ' ido istype ORO=1 ', ido, istype, ' ugwp_oro  as a solver '
       endif

       if (ido == 1 .and. istype ==2 ) then
!
! convective gw effects
!
! 1. specify spectra + forcing   nstcon, nwcon, ch_conv, nazcon, spf_conv
!
         call get_spectra_tau_convgw                                &
             (nwcon, im, levs, dcheat,  scheat, precip, cld_klevs,  &
              xlatd, sinlat, coslat, taub_con, klev_con, if_con, nf_con)
!
! 2. solve for GW effects on the mean flow 
!
         if ( nf_con > 0) then 
 
           klev_con(:) = 52    ! ~5 km
!
!eff_con, nstcon, nwcon, ch_conv, nazcon, spf_conv, xaz_conv, yaz_conv
!
           if (knob_ugwp_solver == 1) call gw_solver_linsatdis                 &
           (im, levs, dtp, kdt, me, taub_con, klev_con, if_con, nf_con,        &  
            nwcon,  ch_conv, nazcon,  spf_conv, xaz_conv, yaz_conv,            &                                 
            fcor, c2f2,    ugrs, vgrs, tgrs, qgrs, prsi, delp,                 &
            prsl, prslk, phii, phil,                                           &
            axc, ayc, eds_c, kdis_c, wtauz) 
 

           if (knob_ugwp_solver == 2) then
!   	     print *,  ' before CONV-2 ', ido, istype, ' gw_solver_wmsdis ', knob_ugwp_solver
             call gw_solver_wmsdis                                               &
             (im, levs, dtp, kdt, me, taub_con, klev_con, if_con, nf_con,        &
              nwfj,  ch_fjet, nazfj,  spf_fjet, xaz_fjet, yaz_fjet,              &                                 
              fcor, c2f2,    ugrs, vgrs, tgrs,                                   &
              qgrs, prsi, delp, prsl, prslk, phii, phil,                         &
              axc, ayc, eds_c, kdis_c, wtauz)
!	      print *,  ' after ido istype CONV-2 ', ido, istype, ' gw_solver_wmsdis ', knob_ugwp_solver
           endif
 
           dudt = dudt + axc     * ugwp_effac
           dvdt = dvdt + ayc     * ugwp_effac
           dtdt = dtdt + eds_c   * ugwp_effac
           kdis = kdis + kdis_c  * ugwp_effac
 
           tauz_ngw = wtauz

         endif

       endif

       if (ido == 1 .and. istype ==3 ) then
!
! nonstationary gw effects
!
! 1. specify spectra + forcing
!
         call get_spectra_tau_nstgw (nwfj, im, levs, &
                    trig_fgf, xlatd, sinlat, coslat, taub_fj, klev_fj, if_fj, nf_fj)
!
! 2. solve for GW effects on the mean flow 
!
         print *, ' tau_nstgw nf_fj-GW triggers ',  nf_fj, ' ugwp_solver = ', knob_ugwp_solver
         if ( nf_fj > 0) then 

           if (knob_ugwp_solver == 1) call gw_solver_linsatdis                 &
           (im, levs, dtp, kdt, me, taub_fj, klev_fj, if_fj, nf_fj,            &  
            nwfj,  ch_fjet, nazfj,  spf_fjet, xaz_fjet, yaz_fjet,              &                                 
            fcor, c2f2,    ugrs, vgrs, tgrs,                                   &
            qgrs, prsi, delp, prsl, prslk, phii, phil,                         &
            axf, ayf, eds_f, kdis_f, wtauz)
 
 

           if (knob_ugwp_solver == 2) call gw_solver_wmsdis                    &
           (im, levs, dtp, kdt, me, taub_fj, klev_fj, if_fj, nf_fj,            &  
            nwfj,  ch_fjet, nazfj,  spf_fjet, xaz_fjet, yaz_fjet,              &                                 
            fcor, c2f2,    ugrs, vgrs, tgrs,                                   &
            qgrs, prsi, delp, prsl, prslk, phii, phil,                         &
            axf, ayf, eds_f, kdis_f, wtauz)

           dudt = dudt + axf    * ugwp_effac
           dvdt = dvdt + ayf    * ugwp_effac
           dtdt = dtdt + eds_f  * ugwp_effac
           kdis = kdis + kdis_f * ugwp_effac  
           tauz_ngw = wtauz
           print *,  ' ido istype for FJ 1-4 ', ido, istype, ' gw_solver_wmsdis ', knob_ugwp_solver

         endif
       endif
!            print *,  ' ido istype for okw 1-4 ', ido, istype
       if (ido == 1 .and. istype == 4 ) then
!
! nonstationary gw effects due to both "convection +fronts/jets " = imbalance of rs-flow
!
! 1. specify spectra + forcing
!
         call get_spectra_tau_okw (nwokw, im, levs,&
              trig_okw, xlatd, sinlat, coslat, taub_okw, klev_okw, if_okw, nf_okw)
!
! 2. solve for GW effects on the mean flow
!
         if ( nf_okw > 0) then
!
           if (knob_ugwp_solver == 1) call gw_solver_linsatdis                 &
           (im, levs, dtp, kdt, me, taub_okw, klev_okw, if_okw, nf_okw,        &  
            nwfj,  ch_fjet, nazfj,  spf_fjet, xaz_fjet, yaz_fjet,              &                                 
            fcor, c2f2,    ugrs, vgrs, tgrs,                                   &
            qgrs, prsi, delp, prsl, prslk, phii, phil,                         &
            axf, ayf, eds_f, kdis_f, wtauz)


           if (knob_ugwp_solver == 2) call gw_solver_wmsdis                    &
           (im, levs, dtp, kdt, me, taub_okw, klev_okw, if_okw, nf_okw,        &  
            nwfj,  ch_fjet, nazfj,  spf_fjet, xaz_fjet, yaz_fjet,              &                                 
            fcor, c2f2,    ugrs, vgrs, tgrs,                                   &
            qgrs, prsi, delp, prsl, prslk, phii, phil,                         &
            axf, ayf, eds_f, kdis_f, wtauz)

           dudt = dudt + axf    * ugwp_effac
           dvdt = dvdt + ayf    * ugwp_effac
           dtdt = dtdt + eds_f  * ugwp_effac
           kdis = kdis + kdis_f * ugwp_effac
           tauz_ngw = wtauz
         endif
       endif
!
! broad gw-spectra
!
 356 continue
     enddo
!
! gw-diag only
!
         axtot   = dudt
         aytot   = dvdt
         eps_tot = dtdt

!
! optional rf-damping
!
     if (do_rfdamp) then
!
!
       call rf_damp(im, levs, levs_rf, dtp, rfdis, rfdist, ugrs, vgrs, ax_rf, ay_rf, eps_rf)
!
!     gw-diag only + rf-damping ..... now orchestrate it with FV3-dycore RF-damping
!
       do k=levs_rf, levs

          dudt(:,k) = dudt(:,k) + ax_rf(:,k)
          dvdt(:,k) = dvdt(:,k) + ay_rf(:,k)
          dtdt(:,k) = dtdt(:,k) + eps_rf(:,k)

       enddo

     endif
!================================================================================
! To update U-V-T STATE by [dudt dvdt dtdt kdis+rf] => Solve 3-diag VD-equation
!================================================================================
! to do for fv3wam=> 
!                   joint eddy+molecular viscosity/conductivity/diffusion
!                   requires "dqdt" + dudt_vis, dvdt_vis.  dtdt_cond

!     print *, '   cires_ugwp_driver +++++++++++++++++ '
!
    end subroutine cires_ugwp_driver

    
!=============================================


     subroutine cires_ugwp_advance
!-----------------------------------------------------------------------
!
!   options for the day-to-day variable sources/spectra + diagnostics
!           for stochastic "triggers"
!   diagnose GW-source functions * FGF + OKWP + SGO/CONV from IAU-fields
!     or use for stochastic GWP-sources "memory"
!-----------------------------------------------------------------------
      implicit none
!
! update sources
!  a) physics-based triggers for multi-wave
!  b) stochastic-based spectra and amplitudes
!  c) use "memory" on GW-spectra from previous time-step
!  d) update "background" GW  dissipation as needed 
!
     end subroutine cires_ugwp_advance
 
!      
! -----------------------------------------------------------------------
! finalize  of cires_ugwp   (_finalize)
! -----------------------------------------------------------------------


  subroutine cires_ugwp_finalize
!
! deallocate sources/spectra & some diagnostics need to find where "deaalocate them"
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
!-----------------------------------------------------------------------------
!
!cires_ugwp_init_emc    CURRENT-style of GFS-v16 with additional diagnostics
!
!-----------------------------------------------------------------------------
   subroutine cires_ugwp_init_emc (me, master, nlunit, logunit, fn_nml2, &
              lonr, latr, levs, ak, bk, pref, dtp, cdmvgwd, cgwf,    &
              pa_rf_in, tau_rf_in)
!
!  input_nml_file ='input.nml'=fn_nml   ..... OLD_namelist and cdmvgwd(4) Corrected Bug Oct 4
!
    use  ugwp_oro_init,     only :  init_oro_gws
    use  ugwp_conv_init,    only :  init_conv_gws
    use  ugwp_fjet_init,    only :  init_fjet_gws
    use  ugwp_okw_init,     only :  init_okw_gws
    use  ugwp_wmsdis_init,  only :  initsolv_wmsdis
    use  ugwp_lsatdis_init, only :  initsolv_lsatdis
    use  ugwp_wmsdis_init,  only :  ilaunch, nslope, lhmet, lzmax, lzmin, lzstar   
    use  ugwp_wmsdis_init,  only :  tau_min, tamp_mpa        
    implicit none

    integer, intent (in) :: me
    integer, intent (in) :: master
    integer, intent (in) :: nlunit
    integer, intent (in) :: logunit
    integer, intent (in) :: lonr
    integer, intent (in) :: levs
    integer, intent (in) :: latr
    real,    intent (in) :: ak(levs+1), bk(levs+1), pref
    real,    intent (in) :: dtp
    real,    intent (in) :: cdmvgwd(4), cgwf(2)             ! "scaling" controls for "old" GFS-GW schemes BUG !!!
    real,    intent (in) :: pa_rf_in, tau_rf_in
 
    character(len=64), intent (in) :: fn_nml2
    character(len=64), parameter   :: fn_nml='input.nml'

!    character,  intent (in) :: input_nml_file
!    integer, parameter :: logunit =  6
    integer :: ios
    logical :: exists
    real    :: dxsg
    integer :: k
!
    if (me == master) print *, trim (fn_nml), ' GW-namelist file '
    inquire (file =trim (fn_nml) , exist = exists)
!
    if (.not. exists) then
       if (me == master) &
        write (6, *) 'separate ugwp :: namelist file: ', trim (fn_nml), ' does not exist'
    else
        open (unit = nlunit, file = trim(fn_nml), action = 'read', status = 'old', iostat = ios)
    endif
    rewind (nlunit)
    read   (nlunit, nml = cires_ugwp_nml)
    close  (nlunit)
!
    ilaunch = launch_level
    pa_rf   = pa_rf_in
    tau_rf  = tau_rf_in

! write version number and namelist to log file
    if (me == master) then
        write (logunit, *) " ================================================================== "
        write (logunit, *) "  cires_ugwp_namelist_emc_v0 short-list "
        write (logunit, nml = cires_ugwp_nml)
    endif
!
! effective kxw - resolution-aware
!
    dxsg =  pi2*arad/float(lonr) * knob_ugwp_ndx4lh
!
!   kxw  =  pi2/dxsg
!
! init global background dissipation for ugwp -> 4d-variable for fv3wam linked with pbl-vert_diff
!

!   allocate(fcor(latr), fcor2(latr)  )
!
    allocate( kvg(levs+1),   ktg(levs+1)  )
    allocate( krad(levs+1),  kion(levs+1) )        
    allocate( zkm(levs),   pmb(levs) )
    allocate( rfdis(levs), rfdist(levs) )
!
! ak -pa  bk-dimensionless  from surf => tol_lid_pressure =0
!
    do k=1, levs
       pmb(k) = 1.e0*(ak(k) + pref*bk(k))    ! Pa -unit  Pref = 1.e5
       zkm(k) = -hpskm*alog(pmb(k)/pref)
    enddo
!
! Part-1 :init_global_gwdis
!
    call init_global_gwdis(levs, zkm, pmb, kvg, ktg, krad, kion, pa_rf, tau_rf, me, master)
    call rf_damp_init     (levs, pa_rf, tau_rf, dtp, pmb, rfdis, rfdist, levs_rf)
!
! Part-2 :init_SOURCES_gws

!=====================================================
! a)    launch_level 
!=====================================================

   do k=levs, 1, -1
     if (pmb(k) .gt. knob_ugwp_palaunch ) exit
   enddo

      launch_level = max(k-1, 5)    ! above 5-layers from the surface
      ilaunch =   launch_level   
!    
! call init-solver for "stationary" multi-wave spectra and sub-grid oro
!
    call init_oro_gws( knob_ugwp_wvspec(1), knob_ugwp_azdir(1), &
         knob_ugwp_stoch(1), knob_ugwp_effac(1), lonr, kxw, cdmvgwd )
!
! call init-sources for "non-sationary" multi-wave spectra
!
    do_physb_gwsrcs=.true.

    IF (do_physb_gwsrcs) THEN

      if (me == master) print *, ' do_physb_gwsrcs ',  do_physb_gwsrcs, ' in cires_ugwp_init '
      if (knob_ugwp_wvspec(4) > 0) then
! okw
        call init_okw_gws(knob_ugwp_wvspec(4), knob_ugwp_azdir(4), &
                          knob_ugwp_stoch(4), knob_ugwp_effac(4), lonr, kxw )
        if (me == master) print *, ' init_okw_gws '
      endif

      if (knob_ugwp_wvspec(3) > 0) then
! fronts
        call init_fjet_gws(knob_ugwp_wvspec(3), knob_ugwp_azdir(3), &
                           knob_ugwp_stoch(3), knob_ugwp_effac(3), lonr, kxw )
        if (me == master) print *, ' init_fjet_gws '
      endif

      if (knob_ugwp_wvspec(2) > 0) then
! conv
        call init_conv_gws(knob_ugwp_wvspec(2), knob_ugwp_azdir(2), &
                           knob_ugwp_stoch(2), knob_ugwp_effac(2), lonr, kxw, cgwf )
        if (me == master)   &
           print *, ' init_convective GWs cgwf', knob_ugwp_wvspec(2), knob_ugwp_azdir(2)

      endif

    ENDIF   !IF (do_physb_gwsrcs)

!======================
! Part-3 :init_SOLVERS
! =====================
!
! call init-solvers for "broad" non-stationary multi-wave spectra
!
    if   (knob_ugwp_solver==1) then
!
      call initsolv_lsatdis(me, master, knob_ugwp_wvspec(2), knob_ugwp_azdir(2), &
                            knob_ugwp_stoch(2), knob_ugwp_effac(2), do_physb_gwsrcs, kxw )
    endif
     if   (knob_ugwp_solver==2) then 
       
       nslope = knob_ugwp_nslope             ! the GW sprctral slope at small-m             
       lzstar = knob_ugwp_lzstar
       lzmax  = knob_ugwp_lzmax 
       lzmin  = knob_ugwp_lzmin      
       lhmet  = knob_ugwp_lhmet
       tamp_mpa =knob_ugwp_tauamp           !amplitude for GEOS-5/MERRA-2
       tau_min  =knob_ugwp_taumin             ! min of GW MF 0.25 mPa                        
       ilaunch = launch_level
       kxw = pi2/lhmet
       
       call initsolv_wmsdis(me, master, knob_ugwp_wvspec(2), knob_ugwp_azdir(2), &
                            knob_ugwp_stoch(2), knob_ugwp_effac(2), do_physb_gwsrcs, kxw)
     endif
!
! other solvers not yet tested for fv3gfs
!
!<    if   (knob_ugwp_solver==3) call init_dspdis
!<    if   (knob_ugwp_solver==4) call init_adodis
! 

!======================
    module_is_initialized = .true.
    if (me == master)   print *, ' EMC-ugwp CDMB(4) is initialized ', module_is_initialized

    end subroutine cires_ugwp_init_emc
!
!
!
!
    subroutine calendar_ugwp(yr, mm, dd, ddd_ugwp)
!
! computes day of year to get tau_limb forcing written with 1-day precision
!    
    implicit none
    integer, intent(in) :: yr, mm, dd
    integer :: ddd_ugwp
    
    integer ::  iw3jdn
    integer :: jd1, jddd
    jd1     = iw3jdn(yr,1,1)
    jddd    = iw3jdn(yr,mm,dd)
    ddd_ugwp = jddd-jd1+1
    
    end subroutine calendar_ugwp
    
    
    subroutine cires_indx_ugwp (npts, me, master, dlat,j1_tau,j2_tau, w1_j1tau, w2_j2tau, &
                                j1_qbo,j2_qbo, w1_j1qbo, w2_j2qbo, dexp_latqbo  )   
	 
    implicit none
!    
! ntau_d1y, ntau_d2t, nqbo_d1y, nqbo_d2z, nqbo_d3t
! ugwp_taulat(:),  ugwp_qbolat(:), ugwp_merlat(:)
!
      integer    ::   npts, me, master
      integer, dimension(npts) ::  j1_tau,j2_tau, j1_qbo, j2_qbo
      real ,   dimension(npts) ::  dlat, w1_j1tau, w2_j2tau, w1_j1qbo, w2_j2qbo
      real ,   dimension(npts) ::  dexp_latqbo  
      real    :: widqbo2, xabs          
!
      integer i,j, j1, j2
!      
! weights for tau_limb  w1_j1tau, w2_j2tau
!
      do j=1,npts
        j2_qbo(j) = nqbo_d1y
        do i=1, nqbo_d1y
          if (dlat(j) < ugwp_qbolat(i)) then
            j2_qbo(j) = i
            exit
          endif
        enddo
	
        
        j2_qbo(j) = min(j2_qbo(j),nqbo_d1y)
        j1_qbo(j) = max(j2_qbo(j)-1,1)	
	
        if (j1_qbo(j) /= j2_qbo(j) ) then
          w2_j2qbo(j) = (dlat(j)  - ugwp_qbolat(j1_qbo(j))) &
                 / (ugwp_qbolat(j2_qbo(j))-ugwp_qbolat(j1_qbo(j)))
       	 
        else
          w2_j2qbo(j) = 1.0
        endif
          w1_j1qbo(j) = 1.0 - w2_j2qbo(j)	
	  
!
      enddo	 
!	   
! weights for tau_limb  w1_j1tau, w2_j2tau
!
      do j=1,npts
        j2_tau(j) = ntau_d1y
        do i=1,ntau_d1y
          if (dlat(j) < ugwp_taulat(i)) then
            j2_tau(j) = i
            exit
          endif
        enddo
	
      
        j2_tau(j) = min(j2_tau(j),ntau_d1y)
        j1_tau(j) = max(j2_tau(j)-1,1)	
	
        if (j1_tau(j) /= j2_tau(j) ) then
          w2_j2tau(j) = (dlat(j)  - ugwp_taulat(j1_tau(j))) &
                 / (ugwp_taulat(j2_tau(j))-ugwp_taulat(j1_tau(j)))
       	 
        else
          w2_j2tau(j) = 1.0
        endif
          w1_j1tau(j) = 1.0 -	w2_j2tau(j)	

      enddo
 		widqbo2 =1./widqbo/widqbo     
       do j=1,npts         
       	    dexp_latqbo(j) =0.		 		 
	    xabs =abs(dlat(j))
	if (xabs .le. latqbo) then
	    dexp_latqbo(j) = exp(-xabs*xabs*widqbo2)
	    if (xabs .le. 4.0 ) dexp_latqbo(j)  =1.
!	    print *, ' indx_ugwp dexp=', dexp_latqbo(j), nint(dlat(j))
        endif		 
       enddo     
      
      if (me == master ) then
222   format( 2x, 'vay-wqbo',  I4, 5(2x, F10.3))     
223   format( 2x, 'vay-limb',  I4, 5(2x, F10.3))
      print *, 'vay_indx_ugwp ',  size(dlat),  ' npts ', npts
        do j=1,npts
	  j1 =  j1_tau(j)
	  j2 =  j2_tau(j)	  
	 write(6,223) j,  ugwp_taulat(j1),  dlat(j),  ugwp_taulat(j2), w2_j2tau(j), w1_j1tau(j)
	enddo 
	print *
        do j=1,npts
	  j1 =  j1_qbo(j)
	  j2 =  j2_qbo(j)	  
	 write(6,222) j,  ugwp_qbolat(j1),  dlat(j),  ugwp_qbolat(j2), w2_j2qbo(j), w1_j1qbo(j)
	enddo          
      endif 	 
    end subroutine cires_indx_ugwp 

!    
 end module cires_ugwp_module

