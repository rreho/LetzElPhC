!/******** Must be compiled with same one used for Quantum ESPRESSO *********/
program yambo_rw
  USE HDF5
  implicit none
  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  INTEGER, PARAMETER :: SP = KIND(1.0) ! Size of floating point number
  character(300) :: filein != '../elph_dir/s.dbph_000001'

  real, parameter :: Qe2Yam = 1 !(1./2.)**(3./2.)
  
  integer :: ik, inu
  
  integer :: nbranch_qe, nk_qe, nbands_qe ! Number of branches, k-points, bands
  real(KIND = DP) :: alat ! Lattice constant
  real(KIND = DP), dimension(:, :), allocatable :: klist ! List of k points
  real(KIND = DP), dimension(3) :: q ! q point
  real(KIND = DP), dimension(:), allocatable ::  omega_ph_sq ! omega^2
  
  ! Order of matrix dimensions: kpt, final band, initial band, branch
  complex(KIND = DP), dimension(:,:,:,:), allocatable ::  g ! El-ph matrix elements
  complex(KIND = DP), dimension(:,:,:), allocatable ::  v ! Ph. eigenvectors
  complex(KIND = DP), dimension(:,:,:,:,:), allocatable ::  dV_dR ! MEs of Gradient of V_scf (in Ha^(3/2) (?))
  real(KIND = DP), dimension(:,:), allocatable :: e_k ! Band energies at k
  real(KIND = DP), dimension(:,:), allocatable :: e_kpq ! Band energies at k+q

  ! SP 
  real(KIND = SP), dimension(:, :), allocatable :: klist_s ! List of k points
  real(KIND = SP), dimension(3) :: q_s ! q point
  real(KIND = SP), dimension(:), allocatable ::  omega_ph_sq_s ! omega^2
  
  ! Order of matrix dimensions: kpt, final band, initial band, branch
  real(KIND = SP), dimension(:,:,:,:,:), allocatable ::  g_s ! El-ph matrix elements
  real(KIND = SP), dimension(:,:,:,:), allocatable ::  v_s ! Ph. eigenvectors
  real(KIND = SP), dimension(:,:,:,:,:,:), allocatable ::  dV_dR_s ! MEs of Gradient of V_scf (in Ha^(3/2) (?))
  real(KIND = SP), dimension(:,:), allocatable :: e_k_s ! Band energies at k
  real(KIND = SP), dimension(:,:), allocatable :: e_kpq_s ! Band energies at k+q
  
  logical :: lgamma ! Are we at q=0?
  
  !#######################################################################
  ! datasets identifiers declaration
  integer(HID_T) :: file_id    ! file identifier
  integer(HID_T) :: dset_id_elph_r    ! dataset identifier for elph_real part
  integer(HID_T) :: dset_id_w2    ! dataset identifier for omega^2
  integer(HID_T) :: dset_id_dyn_r    ! dataset identifier for real part of eigenvecs
  integer(HID_T) :: dset_id_eqk    ! dataset identifier for energies at k+q
  !#######################################################################
  !#######################################################################
    ! dataspace identifiers declaration
  integer(HID_T) :: dataspace_id_elph  ! dataspace identifier for elph
  integer(HID_T) :: dataspace_id_w2  ! dataspace identifier for omega^2
  integer(HID_T) :: dataspace_id_dyn  ! dataspace identifier for eigenvecs
  integer(HID_T) :: dataspace_id_eqk  ! dataspace identifier for energies at k+q
    !#######################################################################
  !#######################################################################
  integer(HSIZE_T), dimension(5) :: elph_dims  !dataspace dims
  integer(HSIZE_T), dimension(1) :: w2_dims  !dataspace dims
  integer(HSIZE_T), dimension(4) :: dyn_dims  !dataspace dims
  integer(HSIZE_T), dimension(2) :: eqk_dims  !dataspace dims
  !### space ranks 
  INTEGER :: error, space_rank_elph, space_rank_w2,&
    space_rank_q, space_rank_k, space_rank_ek, space_rank_eqk, &
    space_rank_dyn, arglen, stat

  ! Get the first command-line argument
  if (command_argument_count() < 1) then
    print *, "Error: Please provide a filename as argument"
    stop
  end if

  call get_command_argument(1, filein, arglen, stat)
  if (stat /= 0) then
    print *, "Error reading command-line argument"
    stop
  end if

  print *, "Reading from file: ", trim(filein)
  
  open(unit=99, file=filein, form='unformatted')
  
  read (99) nbranch_qe, nk_qe, nbands_qe
  !Allocation of arrays
  allocate(klist(3, nk_qe))
  allocate(omega_ph_sq(nbranch_qe))
  allocate(g( nbands_qe, nbands_qe, nbranch_qe,nk_qe))
  allocate(v(nbranch_qe, nbranch_qe/3, 3))
  allocate(dV_dR(nk_qe, nbands_qe, nbands_qe, nbranch_qe/3, 3))
  allocate(e_k( nbands_qe,nk_qe))
  allocate(e_kpq(nbands_qe,nk_qe))
  !
  elph_dims = (/2,nbranch_qe,nbands_qe,nbands_qe,nk_qe/)
  !defining the dimensions
  w2_dims = (/nbranch_qe/)
  dyn_dims = (/2, nbranch_qe, nbranch_qe/3, 3/)
  eqk_dims = (/nbands_qe,nk_qe/)


  read (99) alat, q, klist
  read (99) omega_ph_sq
  
  IF (q(1)==0 .AND. q(2)==0  .AND. q(3)==0 ) THEN 
    lgamma = .true.
  ELSE 
    lgamma = .false.
  ENDIF
  
  do ik = 1, nk_qe
      read (99) g(:,:,:,ik) !(final,initial,nu,k)  !Phonon branch, final band, initial band, k-point
      read (99) v(:,:,:)
      if (lgamma) then
          read (99) dV_dR(ik,:,:,:,:)
      end if
      read (99) e_k(:,ik)
      read (99) e_kpq(:,ik)
  end do
  
  close(99)
  deallocate(dV_dR)
  !!!!!! Uncomment the below line to convert to Yambo Units.
  g(:,:,:,:) = Qe2Yam*g(:,:,:,:)    
  allocate(g_s(2,nbranch_qe,nbands_qe,nbands_qe,nk_qe))
  do inu = 1, nbranch_qe
    g_s(1,inu,:,:,:) = real(g(:,:,inu,:) , kind=SP)
    g_s(2,inu,:,:,:) = real(aimag(g(:,:,inu,:)), kind=SP)
  end do
  deallocate(g)
  allocate(klist_s(3, nk_qe))
  allocate(omega_ph_sq_s(nbranch_qe))
  allocate(e_k_s(nbands_qe,nk_qe))
  allocate(e_kpq_s(nbands_qe,nk_qe ))
  !!
  omega_ph_sq_s = real(omega_ph_sq*0.25,kind=SP)  !!0.25 for Ry^2->Ha^2
  klist_s = real(klist,kind=SP)
  e_k_s = real(e_k*0.5,kind=SP)                   !!0.5 for Ry->Ha
  e_kpq_s = real(e_kpq*0.5,kind=SP)               !!0.5 for Ry->Ha
  !!
  deallocate(omega_ph_sq)
  deallocate(klist)
  deallocate(e_k)
  deallocate(e_kpq)
  !!
  allocate(v_s(2, nbranch_qe, nbranch_qe/3, 3))
  v_s(1,:,:,:) = real(v,kind=SP)
  v_s(2,:,:,:) = real(aimag(v),kind=SP)
  deallocate(v) 
  !allocate(dV_dR_s(nk_qe, nbands_qe, nbands_qe, nbranch_qe/3, 3,2))
  !dV_dR_s(:,:,:,:,:,2) = real(dV_dR,kind=SP)
  !dV_dR_s(:,:,:,:,:,2) = real(aimag(dV_dR),kind=SP)
  !!
  !!deallocate(dV_dR)
  !!
  
  !******** Writing the Data in HDF5 file **************
space_rank_elph = 5
space_rank_w2 = 1
space_rank_dyn = 4
space_rank_eqk = 2
call h5open_f(error) 
call h5fcreate_f('elph_yambo_'//trim(filein)//'.h5',H5F_ACC_TRUNC_F, file_id, error) !create file
! ***** create_dataspaces *******************
call h5screate_simple_f(space_rank_elph,elph_dims,dataspace_id_elph,error)  !Dataspace
call h5screate_simple_f(space_rank_w2,w2_dims,dataspace_id_w2,error)  !Dataspace
call h5screate_simple_f(space_rank_dyn,dyn_dims,dataspace_id_dyn,error)  !Dataspace
call h5screate_simple_f(space_rank_eqk,eqk_dims,dataspace_id_eqk,error)  !Dataspace

! ************ Creating datasets *************
call h5dcreate_f(file_id, "ELPH_GKKP_Q1", H5T_NATIVE_REAL,dataspace_id_elph,dset_id_elph_r,error )
call h5dcreate_f(file_id, "PH_FREQS1", H5T_NATIVE_REAL,dataspace_id_w2,dset_id_w2,error )
call h5dcreate_f(file_id, "POLARIZATION_VECTORS", H5T_NATIVE_REAL,dataspace_id_dyn,dset_id_dyn_r,error )
call h5dcreate_f(file_id, "E_K_PLUS_Q1", H5T_NATIVE_REAL,dataspace_id_eqk,dset_id_eqk,error )
!************ Writing data into datasets *****************
call h5dwrite_f(dset_id_elph_r, H5T_NATIVE_REAL,g_s,elph_dims,error )
call h5dwrite_f(dset_id_w2, H5T_NATIVE_REAL,omega_ph_sq_s,w2_dims,error )
call h5dwrite_f(dset_id_dyn_r, H5T_NATIVE_REAL,v_s,dyn_dims,error )
call h5dwrite_f(dset_id_eqk, H5T_NATIVE_REAL,e_kpq_s,eqk_dims,error )
!******* Closing datasets,dataspaces and files
!closing Datasets
call h5dclose_f(dset_id_elph_r, error)
call h5dclose_f(dset_id_w2, error)
call h5dclose_f(dset_id_dyn_r, error)
call h5dclose_f(dset_id_eqk, error)
!closing Dataspace
call h5sclose_f(dataspace_id_elph, error)
call h5sclose_f(dataspace_id_w2, error)
call h5sclose_f(dataspace_id_dyn, error)
call h5sclose_f(dataspace_id_eqk, error)
!closing file
call h5fclose_f(file_id, error)
call h5close_f(error)
!***** Deallocation *********************
  
  deallocate(klist_s)
  deallocate(omega_ph_sq_s)
  deallocate(g_s)
  deallocate(v_s)
  !deallocate(dV_dR_s)
  deallocate(e_k_s)
  deallocate(e_kpq_s)
  
  end program
  
