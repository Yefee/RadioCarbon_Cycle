module c14_cycle
!-----------------------------------------------------------------------------------------
! Purpose : make a frame work for 
! version 1:
! C14 cycle in CAM, 2 functions are lost,
! 1. c14 production
! 2. C14 decay
! whether or not read in surface data flux is not determined.
! Author :CF
! 09/27/2016
!-----------------------------------------------------------------------------------------
! version 2:
! milestone!
! right now, the production has been added and chekced, but scale factor is one, we need to 
! negotiate with ucar team to determine the sf scheme;
! decay process added.
! coupling past is abscent.
! Anyway, it can be used as production now!!
! 2 initial condition settings:
! 1. co2 mixing ratio ==> c14 mixing ratio, that infers globally well mixed, one value.
! 2. run co2 cycle first ==> co2 mixing ratio varies spatially ==> c14 mixing ratio varies spatially.
! 
! 
!                       _oo0oo_
!                      o8888888o
!                      88" . "88
!                      (| -_- |)
!                      0\  =  /0
!                    ___/`---'\___
!                  .' \\|     |// '.
!                 / \\|||  :  |||// \
!                / _||||| -:- |||||- \
!               |   | \\\  -  /// |   |
!               | \_|  ''\---/''  |_/ |
!               \  .-\__  '-'  ___/-. /
!             ___'. .'  /--.--\  `. .'___
!          ."" '<  `.___\_<|>_/___.' >' "".
!         | | :  `- \`.;`\ _ /`;.`/ - ` : | |
!         \  \ `_.   \_ __\ /__ _/   .-` /  /
!     =====`-.____`.___ \_____/___.-`___.-'=====
!                       `=---='
!
!
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!               佛祖保佑         永无BUG
!

!------------------------------------------------------------------------------------------------

use shr_kind_mod,   only: r8 => shr_kind_r8
use spmd_utils,     only: masterproc ! master process
use ppgrid,         only: pver ! vertical level
use physics_types,  only: physics_state, physics_ptend, physics_ptend_init
use physconst,      only: mwdry, mwco2, gravit, cpair, avogad, pi, gravit 
use constituents,   only: cnst_add, cnst_get_ind, cnst_name, cnst_longname, sflxnam 
use chem_surfvals,  only: chem_surfvals_get !get the surface mixxing ratio
use cam_abortutils,  only: endrun

implicit none
private
save
 
! Public interfaces
public c14_cycle_readnl              ! read the namelist 
public c14_register                  ! register consituents
public c14_transport                 ! turn on co2 tracers transport
public c14_implements_cnst           ! returns true if consituent is implemented by this package
public c14_init_cnst                 ! initialize mixing ratios if not read from initial file
public c14_init                      ! initialize (history) variables
public c14_timestep_init		      ! place to perform per timestep initialization
!public findIndex					  ! find index in c14_lookup_table

! Public type
public table

! Public data								 
public c14_i                           ! global index for new constituents
public c14_table					   ! c14 production lookup table


! Namelist variables
logical :: c14_flag            = .false.      ! true => turn on c14 code, namelist variable
character(len=256) :: c14_global_production  = 'unset' ! c14 global production rate
character(len=256) :: c14_lookup_table = 'unset' ! c14 lookuptable  
character(len=256) :: Phi_and_VADM_and_C14Production = 'const'     ! if condst, Phi_time_series and M_time_series are not useful
real(r8) :: phi_val = 550.0	! Phi for snapshoot run, default 550
real(r8) :: vadm_val = 1.0	! normilized VADM for snapshoot run, default 1.0
real(r8) :: C14prod = 2.0   ! gcm^-2s^-1
character(len=256) :: Phi_time_series = 'unset' ! Phi time series                    
character(len=256) :: M_time_series = 'unset' ! normolized VADM time series                    


type :: table

	integer :: M_i != 9			! length of M dim
	integer :: Phi_i != 31		! length of Phi dim
	integer :: lat_i != 9			! length of lat dim
	integer :: depth_i != 34		! length of depth dim

!	real(r8) :: prodc14(34, 9, 31, 9) ! C14 production table
	real(r8), pointer :: prodc14(:, :, :, :) ! C14 production table
!	real(r8) :: depth(34)						 ! the depth dimension, 1st dim
	real(r8), pointer :: depth(:)						 ! the depth dimension, 1st dim
!	real(r8) :: lat(9)						 ! the lat dimension, 2nd dim
	real(r8), pointer :: lat(:)						 ! the lat dimension, 2nd dim
!	real(r8) :: phi(31)						 ! the Phi dimension, 3rd dim
	real(r8), pointer :: phi(:)						 ! the Phi dimension, 3rd dim
!	real(r8) :: m(9) 							 ! the M dimension, 4th dim
	real(r8), pointer :: m(:) 							 ! the M dimension, 4th dim

end type table

type (table) :: c14_table


!-----------------------------------------------------------------------
! new constituents
integer, parameter :: ncnst=3                      ! number of constituents implemented

character(len=8), dimension(ncnst), parameter :: & ! constituent names, ocn, fuel, land, atm
	c_names = (/'C14_abio', 'CO2_abio', 'D14Cabio'/)
real(r8), dimension(ncnst), parameter :: &         ! molecular weights
	c_mw = (/mwco2 + 2, mwco2, mwco2/)
real(r8), dimension(ncnst), parameter :: &         ! heat capacities
	c_cp = (/cpair, cpair, cpair/)
real(r8), dimension(ncnst), parameter :: &         ! minimum mmr (what is mmr??)
	c_qmin = (/1.e-20_r8, 1.e-20_r8, -1.e6_r8/)
integer, dimension(ncnst) :: c14_i                 ! global index


!================================================================================================
contains
!================================================================================================

subroutine c14_cycle_readnl(nlfile)

	! Read co2_cycle_nl namelist group.

	use namelist_utils,  only: find_group_name
	use units,           only: getunit, freeunit
	use mpishorthand
	implicit none
	character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

	! Local variables
	integer :: unitn, ierr, i
	character(len=*), parameter :: subname = 'c14_cycle_readnl'

	namelist /c14_cycle_nl/ c14_flag, c14_global_production, c14_lookup_table, &
							Phi_time_series, M_time_series, Phi_and_VADM_and_C14Production, phi_val, vadm_val, C14prod
	!-----------------------------------------------------------------------------

	if (masterproc) then
		unitn = getunit()
		open( unitn, file=trim(nlfile), status='old' )
		call find_group_name(unitn, 'c14_cycle_nl', status=ierr)
		if (ierr == 0) then
			read(unitn, c14_cycle_nl, iostat=ierr)
			if (ierr /= 0) then
				call endrun(subname // ':: ERROR reading namelist')
			end if
		end if
		close(unitn)
		call freeunit(unitn)
	end if

#ifdef SPMD
	! Broadcast namelist variables
	call mpibcast (c14_flag,                               1,   mpilog,  0, mpicom)
	call mpibcast (c14_global_production,   len(c14_global_production),   mpichar, 0, mpicom)
	call mpibcast (c14_lookup_table, len(c14_lookup_table),   mpichar, 0, mpicom)
	call mpibcast (Phi_time_series, len(Phi_time_series),   mpichar, 0, mpicom)
	call mpibcast (M_time_series, len(M_time_series),   mpichar, 0, mpicom)
	call mpibcast (Phi_and_VADM_and_C14Production, len(Phi_and_VADM_and_C14Production),   mpichar, 0, mpicom)
	call mpibcast (phi_val, 1,   mpir8, 0, mpicom)
	call mpibcast (vadm_val, 1,   mpir8, 0, mpicom)
	call mpibcast (C14prod, 1,   mpir8, 0, mpicom)
	
#endif

end subroutine c14_cycle_readnl

!================================================================================================

subroutine c14_register
!----------------------------------------------------------------------- 
! 
! Purpose: register advected constituents 
! 
!-----------------------------------------------------------------------
	implicit none
	integer  :: i

	if (.not. c14_flag) return
 
! CO2 as dry tracer
	do i = 1, ncnst
		call cnst_add(c_names(i), c_mw(i), c_cp(i), c_qmin(i), c14_i(i), longname=c_names(i), mixtype='dry')
	end do

end subroutine c14_register

!================================================================================================

function c14_transport()

!-----------------------------------------------------------------------
 
! Purpose: return true if this package is active

!-----------------------------------------------------------------------
	logical :: c14_transport
!-----------------------------------------------------------------------

	c14_transport = c14_flag

end function c14_transport

!================================================================================================

function c14_implements_cnst(name)

!----------------------------------------------------------------------- 
! 
! Purpose: return true if specified constituent is implemented by this package
! 
!-----------------------------------------------------------------------
	 implicit none
!-----------------------------Arguments---------------------------------

	 character(len=*), intent(in) :: name  ! constituent name
	 logical :: c14_implements_cnst        ! return value

	 integer :: m     
		
	 c14_implements_cnst = .false.
 
	 if (.not. c14_flag) return
 
	 do m = 1, ncnst
		 if (name == c_names(m)) then
			 c14_implements_cnst = .true.
			 return
		 end if
	 end do
  end function c14_implements_cnst

!===============================================================================  
subroutine c14_init

!----------------------------------------------------------------------- 
! 
! Purpose: initialize c14,
!          declare history variables,
!          read co2 flux form ocn,  as data_flux_ocn
!          read co2 flux form fule, as data_flux_fuel
!
!-----------------------------------------------------------------------

	 use cam_history, only: addfld, add_default, phys_decomp
	 use netcdf
	 use cam_logfile, only: iulog
	 use error_messages, only : alloc_err, handle_ncerr, handle_err
	 use cam_abortutils,   only: endrun
 	 implicit none

	 integer :: m, mm
	 integer :: ncid, pc14id, depid, latid, phiid, mid
	 
	
	 if (.not. c14_flag) return
 
	 ! Add constituents and fluxes to history file
	 do m = 1, ncnst
		
		 call cnst_get_ind(c_names(m), mm)
		if (m < 3) then
		 call addfld(trim(cnst_name(mm))//'_BOT', 'kg/kg',     1, 'A', trim(cnst_longname(mm))//', Bottom Layer', phys_decomp)
		 call addfld(cnst_name(mm),               'kg/kg',  pver, 'A', cnst_longname(mm), phys_decomp)
		 call addfld(sflxnam(mm),                 'kg/m2/s',   1, 'A', trim(cnst_name(mm))//' surface flux', phys_decomp)

		 call add_default(cnst_name(mm), 1, ' ')
		 call add_default(sflxnam(mm),   1, ' ')

		 ! The addfld call for the 'TM*' fields are made by default in the 
		 ! constituent_burden module.
		 call add_default('TM'//trim(cnst_name(mm)), 1, ' ')
		
		else
		 call addfld(trim(cnst_name(mm))//'_BOT', 'permil',     1, 'A', trim(cnst_longname(mm))//', Bottom Layer', phys_decomp)
		 call addfld(cnst_name(mm),               'permil',  pver, 'A', cnst_longname(mm), phys_decomp)
		 call addfld(sflxnam(mm),                 'permil/m2/s',   1, 'A', trim(cnst_name(mm))//' surface flux', phys_decomp)

		 call add_default(cnst_name(mm), 1, ' ')
		 call add_default(sflxnam(mm),   1, ' ')

		 ! The addfld call for the 'TM*' fields are made by default in the 
		 ! constituent_burden module.
		 call add_default('TM'//trim(cnst_name(mm)), 1, ' ')
		
		endif
	 end do
	
 	!---------------------------------------------------------------------------------------
	!
	! read c14 production table and its dimensions
	!
	!---------------------------------------------------------------------------------------
	! open netcdf file
	call handle_ncerr(nf90_open(trim(c14_lookup_table), 0, ncid), 'C14_lookup_table ==> c14_cycle.F90')
	write (iulog,*) 'C14_lookup_table_Read: ncid is', ncid, 'for file', trim(c14_lookup_table)
		
	! Get depth id
	call handle_ncerr(nf90_inq_varid(ncid, 'depth', depid), 'C14_lookup_table_get_depthid ==> c14_cycle.F90')
	! Get depth dim length
	call handle_ncerr( nf90_inquire_dimension( ncid, depid, len=c14_table%depth_i),'C14_lookup_table_get_depth_dim_length ==> c14_cycle.F90')
	allocate(c14_table%depth(c14_table%depth_i))
	! Get depth data
	call handle_ncerr(nf90_get_var(ncid, depid, c14_table%depth), 'C14_lookup_table_get_depthval ==> c14_cycle.F90')
	write (iulog,*) 'Successfully reind in depth'

		
	! Get lat id
	call handle_ncerr(nf90_inq_varid(ncid, 'lat', latid), 'C14_lookup_table_get_latid ==> c14_cycle.F90')
	! Get lat dim length
	call handle_ncerr( nf90_inquire_dimension( ncid, latid, len=c14_table%lat_i),'C14_lookup_table_get_lat_dim_length ==> c14_cycle.F90')
	allocate(c14_table%lat(c14_table%lat_i))
	! Get lat data
	call handle_ncerr(nf90_get_var(ncid, latid, c14_table%lat), 'C14_lookup_table_get_latval ==> c14_cycle.F90')
		
	! Get Phi id
	call handle_ncerr(nf90_inq_varid(ncid, 'Phi', phiid), 'C14_lookup_table_get_phiid ==> c14_cycle.F90')
	! Get  Phi length
	call handle_ncerr( nf90_inquire_dimension( ncid, phiid, len=c14_table%Phi_i),'C14_lookup_table_get_Phi_dim_length ==> c14_cycle.F90')
	allocate(c14_table%phi(c14_table%Phi_i))
	! Get Phi data
	call handle_ncerr(nf90_get_var(ncid, phiid, c14_table%phi), 'C14_lookup_table_get_phival ==> c14_cycle.F90')
		
	! Get normilized VADM id
	call handle_ncerr(nf90_inq_varid(ncid, 'M', mid), 'C14_lookup_table_get_mid ==> c14_cycle.F90')
	! Get VADM dim length
	call handle_ncerr( nf90_inquire_dimension( ncid, mid, len=c14_table%M_i),'C14_lookup_table_get_depth_dim_length ==> c14_cycle.F90')
	allocate(c14_table%m(c14_table%M_i))
	! Get normilized VADM data
	call handle_ncerr(nf90_get_var(ncid, mid, c14_table%m), 'C14_lookup_table_get_mval ==> c14_cycle.F90')
	
	! Get prod_14c id
	call handle_ncerr(nf90_inq_varid(ncid, 'prod_14C', pc14id), 'C14_lookup_table_get_C14id ==> c14_cycle.F90')
	write (iulog,*) 'C14_lookup_table_Read: pc14id is',pc14id, 'for variable', 'prod_14C' 
	allocate(c14_table%prodc14(c14_table%depth_i, c14_table%lat_i, c14_table%Phi_i, c14_table%M_i))
	! Get prod_14c data
	call handle_ncerr(nf90_get_var(ncid, pc14id, c14_table%prodc14), 'C14_lookup_table_get_C14val ==> c14_cycle.F90')
	write (iulog,*) 'C14_lookup_table_Read: sueeceesfully reading in c14table' 

	
	select case (Phi_and_VADM_and_C14Production)
	
		case('const')
			write (iulog,*) 'Sucessfully initializing snapshoot run for constant Phi and VADM'
		
		case ('transcient')
			write (iulog,*) 'Sucessfully initializing transcient run for varying Phi and VADM' 
			
		case default
			write (iulog,*) 'ERROR: Phi_and_VADM_and_C14Production type is not approved, only const and transcient are supported.' 
			call endrun('ERROR: Phi_and_VADM_and_C14Production type is not approved, only const and transcient are supported.')
			
	end select
	 
end subroutine c14_init

!===========================================================================================

subroutine c14_init_cnst(name, q, gcid)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Set initial values of CO2_OCN, CO2_FFF, CO2_LND, CO2
! Need to be called from process_inidat in inidat.F90
! (or, initialize co2 in co2_timestep_init)        
!
!-----------------------------------------------------------------------
	use cam_logfile,    only : iulog
	use cam_abortutils,   only: endrun
! Arguments
	implicit none
	character(len=*), intent(in) :: name         ! constituent name
	real(r8), intent(out) :: q(:,:)   !  mass mixing ratio
	integer, intent(in) :: gcid(:)    ! global column id
!-----------------------------------------------------------------------

	if (.not. c14_flag) return
 
	select case (name)
	case ('C14_abio')
        q =  chem_surfvals_get('CO2MMR') * 1.1_r8 ! 100 permil D14C
		
	case ('CO2_abio')
!		q = 1.4_r8 * 4.704e-16_r8 * 46._r8 / mwdry ! 4.704e-16_r8 is converted from 400ppm CO2, 1.4 mean D14C = 400 permil
		q = chem_surfvals_get('CO2MMR')
	case ('D14Cabio')
		q = 0._r8
	case default
		write (iulog,*) 'ERROR: Phi_and_VADM_and_C14Production type is not approved, only const and transcient are supported.' 
		call endrun('ERROR: Unsupported c_name in c14 cycle.')	
	end select

end subroutine c14_init_cnst
!===============================================================================
 
subroutine c14_timestep_init( phys_state )

!-----------------------------------------------------------------------
! Provides a place to reinitialize diagnostic constituents HORZ and VERT
!-----------------------------------------------------------------------

	use time_manager,   only: get_curr_date, get_step_size, get_curr_calday
	use ppgrid,         only: begchunk, endchunk
	use physics_types,  only: physics_state
	use cam_logfile,    only : iulog
	use cam_abortutils,   only: endrun
	implicit none
	type(physics_state), intent(inout), dimension(begchunk:endchunk), optional :: phys_state    
	integer c, i, k, ncol
	integer yr, mon, day, tod
    integer :: latindex, depindex, phiindex, mindex
	real(r8), allocatable :: tmp(:,:,:)  ! temperary prodc14, (chunk, col, pver)
	integer :: dtime ! time step length (seconds)
	real(r8) :: yrs, calday
	
	! if not
	if (.not. c14_flag) return
	call get_curr_date (yr,mon,day,tod)
!	write (iulog, *) 'Time now is:', yr, 'year;', mon, 'months;', day, 'day;', tod, 'tod.'

	calday = get_curr_calday()
!	write (iulog, *) 'calday now is:', calday

	! time step
	dtime = get_step_size()
	! seconds in year
	yrs = 86400._r8 * 365._r8
	
	select case (Phi_and_VADM_and_C14Production)
	
	case ('const')
  
	! find index for Phi and M_i
	call findIndex(c14_table%phi, phi_val, phiindex, c14_table%Phi_i)
	! write (iulog,*) 'successfully find index for Phi', phiindex 
    
	call findIndex(c14_table%m, vadm_val, mindex, c14_table%M_i)
	! write (iulog,*) 'successfully find index for M_i', mindex 

    allocate(tmp(endchunk-begchunk,ncol,pver)) ! this is used for acculumate c14 production, modify it later.
 
	do c = begchunk, endchunk
		ncol = phys_state(c)%ncol
		do i = 1, ncol				
			call findIndex(c14_table%lat, abs(phys_state(c)%lat(i) * 180._r8 / pi), latindex, c14_table%lat_i)
			do k = 1, pver					
				call findIndex(c14_table%depth, phys_state(c)%pmid(i,k)/10.0_r8/gravit, depindex, c14_table%depth_i) ! roughly devided by 100, modify it later
!				tmp(c,i,k) = c14_table%prodc14(depindex, latindex, phiindex, mindex)
!				phys_state(c)%q(i,k,c14_i(1)) = phys_state(c)%q(i,k,c14_i(1)) + c14_table%prodc14(depindex, latindex, phiindex, mindex) / (avogad / 1.0e3_r8) * 46.0_r8 * dtime
				phys_state(c)%q(i,k,c14_i(1)) = c14_table%prodc14(depindex, latindex, phiindex, mindex) / 1.176e-12_r8 / (avogad / 1.0e3_r8) * 46.0_r8 * dtime + phys_state(c)%q(i,k,c14_i(1)) * (1 - log(2._r8) / 5730._r8 / yrs * dtime)
				phys_state(c)%q(i,k,c14_i(3)) = (phys_state(c)%q(i,k,c14_i(1)) / phys_state(c)%q(i,k,c14_i(2)) - 1._r8) * 1000._r8
				!phys_state(c)%q(i,k,c14_i(1)) = phys_state(c)%q(i,k,c14_i(1)) + 1.2e-9_r8 * sin(phys_state(c)%lat(i))
!				phys_state(c)%q(i,k,c14_i(1)) = 0.0_r8
!			    phys_state(c)%q(i,k,c14_i(2)) = phys_state(c)%q(i,k,c14_i(2)) * (1 - log(2._r8) / 5730._r8 / yrs * dtime) + c14_table%prodc14(depindex, latindex, phiindex, mindex) / (avogad / 1.0e3_r8) * 46.0_r8 * dtime 
			    

			end do
		end do
!		write(iulog,*) 'C14O2 is', phys_state(c)%q(1,pver,c14_i(1))
!		write(iulog,*) 'CO2 is', phys_state(c)%q(1,pver,c14_i(2))
!		write(iulog,*) 'd14C is', phys_state(c)%q(1,pver,c14_i(3))
	end do

	
    deallocate(tmp)

	case ('transcient')
	! this block for trcient run, modify it later!
    write (iulog,*) 'successfully entering select case: transcient' 
	
	
	case default
    write (iulog,*) 'ERROR: Phi_and_VADM_and_C14Production type is not approved, only const and transcient are supported.' 
	call endrun('ERROR: Phi_and_VADM_and_C14Production type is not approved, only const and transcient are supported.')
	
	end select
	
end subroutine c14_timestep_init

subroutine findIndex(x, y, k, sizex)
!------------------------------------------------------------------------------
! This is used for finding index in C14 lookup table, shall not expose to other 
! modules.
!------------------------------------------------------------------------------
    implicit none
	integer, intent(in) :: sizex
	real(r8), intent(in) :: x(sizex)
	real(r8), intent(in) :: y
	integer, intent(out) :: k
	integer :: i

	if (x(1) > y) then
		k = 1
		return
	end if

	if (x(sizex) < y) then
		k = sizex
		return
	end if


	do i = 1, sizex - 1
		if ((x(i) <= y) .and. (y < x(i+1))) then
			if ((abs(x(i) - y) <= abs(x(i+1) - y))) then
				k = i
				return
			else
				k = i + 1
				return
			end if
		end if
	end do
end subroutine findIndex
 
end module c14_cycle
