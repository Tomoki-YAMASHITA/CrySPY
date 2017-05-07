  MODULE  LATTICEVECTOR
  use prec
  implicit none
  real(dp), dimension(3,3):: prim_lat, recip_lat
  real(dp), dimension(3,3):: inv_prim_lat, inv_recip_lat, tran_inv_prim_lat

  END MODULE

#if USE_M_CMD_ARG
  module m_cmd_arg
  use prec 
  implicit none

   contains

   subroutine  get_cmd_arg(rmin,rmax,sigma, npoints, nrepeat_min, nrepeat_max, poscarfilename) 
     use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                          stdout=>output_unit, &
                                          stderr=>error_unit
    implicit none
    real(dp) ::  rmin,  rmax, sigma
    integer :: npoints,nrepeat_min, nrepeat_max
    character(len=256) ::poscarfilename

    integer::  narg,iarg
    character(len=256) :: argv 
    logical :: showparam 
!default values
!     rmin = 0.5d0
!     rmax = 10.0d0
!     npoints = 100
!     sigma = 0.05
!     poscarfilename="POSCAR"

   showparam=.false.
! process command line arguments
    narg = command_argument_count()
    iarg=1
    do while (iarg<=narg)

     call get_command_argument(iarg, argv)

     if (argv(1:1)=="-") then

      select case(argv)
       case ('-rmin')
          iarg = iarg + 1
          call get_command_argument(iarg, argv)
          read(argv,*) rmin 
       case ('-rmax')
          iarg = iarg + 1
          call get_command_argument(iarg, argv)
          read(argv,*) rmax
       case ('-npoints')
          iarg = iarg + 1
          call get_command_argument(iarg, argv)
          read(argv,*) npoints
       case ('-sigma')
          iarg = iarg + 1
          call get_command_argument(iarg, argv)
          read(argv,*) sigma
       case ('-nrepeat_max')
          iarg = iarg + 1
          call get_command_argument(iarg, argv)
          read(argv,*) nrepeat_max
       case ('-nrepeat_min')
          iarg = iarg + 1
          call get_command_argument(iarg, argv)
          read(argv,*) nrepeat_min
       case ('-showparam')
          iarg = iarg + 1
          call get_command_argument(iarg, argv)
          ! always true if there exist the flag
          showparam=.true.
       case default
          write(stderr,*)'unkown argument', trim(argv)
          stop  200
      end select
    else 
       read(argv,*) poscarfilename
    endif 

     iarg = iarg + 1
    enddo

    if (showparam) then 
    write(stderr,*)'rmin=',rmin
    write(stderr,*)'rmax=',rmax
    write(stderr,*)'npoints=',npoints
    write(stderr,*)'sigma=',sigma
    write(stderr,*)'nrepeat_min=',nrepeat_min
    write(stderr,*)'nrepeat_max=',nrepeat_max
    write(stderr,*)'filename=',trim(poscarfilename)
    endif
   end subroutine  get_cmd_arg
  end module m_cmd_arg

#endif 

  program cal_fingerprint 
   use prec
   use constant
   use latticevector
   use strings 
   use convolute
#if USE_M_CMD_ARG
   use m_cmd_arg
#endif
   implicit none
     integer i, j, k,n, istat,ios, nions,ntype, numb, error_flag, ibuf
     integer ii, jj, kk, itype, jtype, iatom,jatom
     integer:: k2
     integer, allocatable, dimension(:):: number_atom 
     character(len=2), allocatable, dimension(:):: element_name
     real(dp), allocatable, dimension(:,:):: pos, pos_c_ang, pos_f 
     real(dp), dimension(3,3):: scaled_prim_lat, metric_tensor
     real(dp) :: ascale, volume, dtemp, halfcell
     real(dp), dimension(3):: vect_pos12
     real(dp), dimension(27):: dist_temp 
     real(dp), allocatable, dimension(:,:) :: dist_at
     character(len=256):: filename, cbuf,title,foutname
     character(len=1):: delims
     character(len=1):: ctype 
     integer, parameter:: StrMax=10, Nmax=30
     character(len=StrMax), dimension(Nmax):: args
     logical :: lexist, isv5, isfc, issd 
     real(dp):: a1m, a2m, a3m, a, b, c, alpha, beta, gama
     real(dp):: sigma
     real(dp):: rmin, rmax, rd, y0
     integer:: npoints, ipoints, npairs, ipairs, nrepeat
     real(dp), allocatable,dimension(:) :: rad, pcf, pcf_temp
     real(dp), allocatable,dimension(:,:) :: pair_pcf 
!!!  For half of maximual cell length
     real(dp) :: rdh
     real(dp), allocatable,dimension(:) :: radh, pcfh, pcfh_temp
     real(dp), allocatable,dimension(:,:) :: pair_pcfh 
    ! Declare local variables
     INTEGER :: AllocateStatus, DeAllocateStatus
     character(len=32) :: cmdarg
     integer :: ncmdarg

! supercell index 
     integer:: nrepeat_min, nrepeat_max 



! default values     
     rmin = 0.5d0
     rmax = 10.0d0
     npoints = 100
     sigma = 0.05 
     nrepeat_min = 1
     nrepeat_max = 10
     filename="POSCAR"

#if USE_M_CMD_ARG
     call get_cmd_arg(rmin,rmax, sigma, npoints,nrepeat_min, nrepeat_max,filename)
     rd= (rmax - rmin)/(npoints-1)
     call compact(filename)

#else

     rd= (rmax - rmin)/(npoints-1)
!
     ncmdarg = command_argument_count() 
     if ( ncmdarg .LT. 1) then
         print '(a)','I will read POSCAR by default.'
         filename = 'POSCAR'
     end if
     
     call get_command_argument(1, cmdarg)

     call compact(cmdarg)
     filename = cmdarg

#endif 
     inquire(file=filename,exist=lexist)
     if ( .not.  lexist) then
          write(*,'(3A)') "Error:",trim(filename), "does not exist"
          STOP
     end if

     open(10,file=filename,status='old')
     read(10,'(a)',iostat=istat) title 
     read(10,*) ascale
     do i = 1, 3
        read(10,*) (scaled_prim_lat(i,j), j = 1, 3)
        !! Unit of ascale in Angstrom
         prim_lat(i,:) = ascale * scaled_prim_lat(i,:)
     end do

     !determine the format of POSCAR or CONTCAR: vasp4.6 or vasp5.x
     read(10,'(a)',iostat=istat) cbuf 
     delims=' '
     
     call parse(cbuf, delims, args, n)
     ntype = n
     allocate(element_name(ntype), number_atom(ntype), & 
&             STAT = AllocateStatus)
     IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

     call value(trim(args(1)), numb, ios)
     if ( ios .NE. 0 ) then 
         isv5 = .TRUE.
         do j = 1, n
            element_name(j) = trim(args(j))
         end do
         read(10,'(a)',iostat=istat) cbuf 
         call parse(cbuf, delims, args, n)
     else
         isv5 = .FALSE.
         do j = 1, n
            ibuf= j
            call writenum(ibuf, cbuf,'I6')
            element_name(j)= trim(cbuf)
         end do
     end if

     nions=0
     do j = 1, n
        call value(trim(args(j)), numb, ios)
        if ( ios .EQ. 0 ) then
          number_atom(j) = numb 
          nions = nions + numb
        end if 
     end do

     npairs = ntype * (ntype+1)/2 
     allocate(pos(nions,3), pos_c_ang(nions,3), pos_f(nions,3),  &
&           dist_at(nions,nions),   &
&     rad(npoints), pcf(npoints), pcf_temp(npoints), pair_pcf(npairs, npoints), &
&     radh(npoints), pcfh(npoints), pcfh_temp(npoints), pair_pcfh(npairs, npoints), &
&     STAT = AllocateStatus)
     IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

     read(10,'(a1)', iostat=istat) ctype
     if ( (ctype .EQ. 'S') .OR. (ctype .EQ. 's')) then
          issd= .TRUE.
          read(10,'(a1)', iostat=istat) ctype
     else 
          issd = .FALSE.
     end if

     if ( (ctype .EQ. 'D') .OR. (ctype .EQ. 'd')) then
           isfc = .TRUE.
     else  
           isfc = .FALSE.
     end if

     do i = 1, nions
        read(10,*) (pos(i,j), j=1,3)
     end do

     close(10)

     call findinv(prim_lat,inv_prim_lat, 3, error_flag)
     tran_inv_prim_lat=TRANSPOSE(inv_prim_lat)

!   get the information of lattice parameters     
     a1m = sqrt(DOT_PRODUCT(prim_lat(1,:),prim_lat(1,:)))
     a2m = sqrt(DOT_PRODUCT(prim_lat(2,:),prim_lat(2,:)))
     a3m = sqrt(DOT_PRODUCT(prim_lat(3,:),prim_lat(3,:)))
     alpha = acos(DOT_PRODUCT(prim_lat(2,:),prim_lat(3,:))/(a2m * a3m))
     beta  = acos(DOT_PRODUCT(prim_lat(1,:),prim_lat(3,:))/(a1m * a3m))
     gama  = acos(DOT_PRODUCT(prim_lat(1,:),prim_lat(2,:))/(a1m * a2m))
     volume = a1m * a2m * a3m * sqrt(1.0d0 + 2.0d0 * cos(alpha) * cos(beta)  &
 &              * cos(gama) - (cos(alpha))**2 - (cos(beta))**2           &
 &              - (cos(gama))**2 ) 

      do i = 1, 3
         do j = 1, 3
          metric_tensor(i,j) = DOT_PRODUCT(prim_lat(i,:), prim_lat(j,:))
         end do
      end do

     vect_pos12 = prim_lat(1,:) + prim_lat(2,:) + prim_lat(3,:)
     halfcell = sqrt(DOT_PRODUCT(vect_pos12,vect_pos12))

  !  convert atomic coordinates into Cartesian type in unit of Angstrom 
     do i = 1, nions 
         if ( isfc ) then
            pos_c_ang(i,:) = MATMUL(pos(i,:),prim_lat)
         else 
            pos_c_ang(i,:) = ascale*pos(i,:)
         end if 
     end do

  !  convert atomic coordinates into fractional type 
     do i = 1, nions 
         if ( isfc ) then
            pos_f(i,:) = pos(i,:)
         else 
            pos_f(i,:) = MATMUL(pos(i,:)*ascale, inv_prim_lat)
         end if 
     end do

      
      nrepeat =3
      
      y0=1.0d0
      rdh = (halfcell - rmin)/(npoints - 1)
      do i = 1, npoints 
          rad(i) = rmin + rd * (i-1)
          radh(i) = rmin + rd * (i-1)
      end do

      iatom = 0
      ipairs = 0

      do itype= 1, ntype
        iatom = number_atom(itype) + iatom
        jatom = 0
        do jtype = 1, ntype
            pcf = 0.0d0
            pcfh = 0.0d0
            jatom = number_atom(jtype) + jatom
             do i = iatom  - number_atom(itype) + 1, iatom
                 do j = jatom - number_atom(jtype) + 1, jatom

                    do nrepeat = nrepeat_min,  nrepeat_max 
                     k = 0
                     k2= 0
                    do  ii = -nrepeat, nrepeat
                        do jj = -nrepeat, nrepeat
                            do  kk = -nrepeat, nrepeat
                                if (nrepeat > nrepeat_min .and. abs(ii) <nrepeat .and.  &
                                abs(jj) <nrepeat .and. abs(kk) <nrepeat ) cycle
                                !write(*,*) 'try ii,jj,kk',ii,jj,kk
                                vect_pos12(1)= pos_f(i,1) - pos_f(j,1) + ii 
                                vect_pos12(2)= pos_f(i,2) - pos_f(j,2) + jj 
                                vect_pos12(3)= pos_f(i,3) - pos_f(j,3) + kk 
                                dtemp= &
 &                                 sqrt(DOT_PRODUCT(MATMUL(vect_pos12, prim_lat), MATMUL(vect_pos12, prim_lat)))
                                if ( (dtemp .GE. rmin) .AND. (dtemp .LE. rmax) ) then
                                   k = k + 1
                                   y0 =volume /(4.0d0*pi*dtemp**2 *  &
 &                                     number_atom(itype) * number_atom(jtype) * rd) 
!                                  y0 =1.0d0/(number_atom(itype)*invsqrt2pi /sigma)
                                   call gaussian_smearing(dtemp, y0, rad, pcf_temp, npoints, sigma)
                                    pcf = pcf_temp + pcf
                                 else 
                                    y0 = 0.0d0
                                 end if
!!! calculate the atom pair correlation function with half of maximual 
!!!     cell length 
                                if ( (dtemp .GE. rmin) .AND. (dtemp .LE. halfcell) ) then
                                   k2=k2+1
                                   y0 =volume /(4.0d0*pi*dtemp**2 *  &
 &                                     number_atom(itype) * number_atom(jtype) * rdh) 
!                                  y0 =1.0d0/(number_atom(itype)*invsqrt2pi /sigma)
                                   call gaussian_smearing(dtemp, y0, radh, pcfh_temp, npoints, sigma)
                                    pcfh = pcfh_temp + pcfh
                                 else 
                                    y0 = 0.0d0
                                 end if
!!!!!
                            end do       ! kk
                        end do           ! jj
                    end do               ! ii

                    write(*,*)'nrepeat=',nrepeat,'k=',k 
                    if (k==0 .and. k2==0) exit 
                    if (k>0 .and. k2>0 .and. nrepeat == nrepeat_max ) then 
                         write(*,*)'fatal error : increase nrepeat_max '
                         stop 300
                    endif 
                    enddo                !   nrepeat

                 end do                  ! j
             end do                      ! i
             if ( itype .LE. jtype ) then
                ipairs = ipairs + 1
                pair_pcf(ipairs,:) = pcf(:)
                pair_pcfh(ipairs,:) = pcfh(:)
             end if
          end do       ! jtype
      end do           ! itype

      open(29,file='feature_ffpf.dat')
!     save F-fingerprint function of 30 atom pairs at most for plot       
!!    calculated within value of rmax 
      write(29, '(6000F15.8)')  ((pair_pcf(ipairs, ipoints)-1.0d0, ipoints =1, npoints), & 
 &                   ipairs = 1, npairs)     
      close(29)

      open(30,file='ffpf4plot.dat')
!     save F-fingerprint function of 30 atom pairs at most for plot       
!!    calculated within value of rmax 
      do ipoints= 1, npoints
         write(30,'(31F12.3)') rad(ipoints), (pair_pcf(ipairs, ipoints) -1.0d0, ipairs = 1,npairs)
      end do
      close(30)

      open(31,file='pcf4plot.dat')
!     save pair correlation function of 30 atom pairs at most for plot       
!!    calculated within value of rmax 
      do ipoints= 1, npoints
         write(31,'(31F12.3)') rad(ipoints), (pair_pcf(ipairs, ipoints), ipairs = 1,npairs)
      end do
      close(31)

      open(32,file='half_pcf4plot.dat')
!     save pair correlation function of 30 atom pairs at most for plot       
!!    calculated within half of maximual cell length
      do ipoints= 1, npoints
         write(32,'(31F12.3)') radh(ipoints), (pair_pcfh(ipairs, ipoints), ipairs = 1,npairs)
      end do
      close(32)

      open(33,file='feature_geo.dat')
!!! volume, metric tensor, fraction positon of each atom, pair pair correlation
!!! function within hal of maximual cell length for each atom pair
      write(33,'(6007F15.8)') volume, metric_tensor(1,1), metric_tensor(1,2), &
 &        metric_tensor(1,3), metric_tensor(2,1), metric_tensor(2,2), &
 &        metric_tensor(2,3), metric_tensor(3,3),                     &
 &        ((pos_f(i,j), j =1, 3), i=1, nions),          &
 &      ((pair_pcfh(ipairs, ipoints),ipoints =1, npoints), ipairs = 1,npairs) 
      close(33)
    
      DEALLOCATE(element_name, number_atom, pos, pos_c_ang, pos_f, dist_at, &
&          rad, pcf, pcf_temp, pair_pcf,   &
&          radh, pcfh, pcfh_temp, pair_pcfh,   &
&          STAT = DeAllocateStatus)
      IF (DeAllocateStatus /= 0) then
           write(*,*) "*** Trouble deallocating ***"
           stop 100
      endif 

      stop 0  ! normal exit
  END program cal_fingerprint 
