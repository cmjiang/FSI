program runFSI
  implicit none
  
#ifdef RUN_MPI
  include "mpif.h"
#endif

  integer, parameter :: my_root = 0
  integer :: myid, numnodes, ierr, t_total_start, t_total_end, rate
  integer :: N, L, b, q, mType, bN, LN, m, ns, t, sigma, U, nnz
  integer :: nproc, OMP_GET_NUM_PROCS, inputUnit, m_per_task, iter, nArg
  integer, allocatable :: h(:), myh(:), ci(:), rp(:)

  double precision :: beta, quantity, lquantity, gquantity, eps
  double precision, allocatable :: matsB(:,:,:), matsBh(:,:,:)
  double precision, allocatable :: matMh(:,:)
  double precision, allocatable :: matsG(:,:,:,:), work(:,:,:)
  double precision, allocatable :: vals(:) 
  character(len = 200) :: inputFile

#ifdef RUN_MPI
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, numnodes, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
#else
  myid = 0
  numnodes = 1
#endif

  ! Read parameters from input file
  nArg = iargc()
  if (nArg < 1) then
    print *, 'syntax: runFSI <input file>'
    stop
  endif
  call getarg(1, inputFile)
  inputUnit = 10
  open(inputUnit, file = inputFile)
  read(inputUnit, *)
  read(inputUnit, *)
  read(inputUnit, *) mType
  read(inputUnit, *) m
  read(inputUnit, *) ns
  read(inputUnit, *) L
  read(inputUnit, *) b
  read(inputUnit, *) q
  read(inputUnit, *) beta
  read(inputUnit, *) t
  read(inputUnit, *) sigma
  read(inputUnit, *) U
  close(inputUnit)

  m_per_task = m/numnodes
  N = ns*ns
  bN = b*N
  LN = L*N
  nnz = LN*(N+1)
  ! Initialize OpenMP environment
  call ompinit()
  
  ! Allocate memory
  allocate(matsB(N,N,L))
  allocate(matsBh(N,N,b))
  allocate(work(N,N,b))
  allocate(matMh(bN,bN))
  allocate(matsG(N,N,L,b))
  allocate(myh(L*N*m_per_task))
  allocate(vals(nnz))
  allocate(ci(nnz))
  allocate(rp(LN+1))
  if(myid == my_root) then
    call system_clock(t_total_start)
    allocate(h(L*N*m))
    call randomh(h,N,L,m)
  endif

#ifdef MPI_RUN
  call MPI_Scatter(h,   L*N*m_per_task, MPI_INTEGER, &
                   myh, L*N*m_per_task, MPI_INTEGER, &
                   my_root,                          &
                   MPI_COMM_WORLD, ierr)
#else
  myh = h
#endif

  ! Main iterations
  do iter = 1, m_per_task
    if (mType == 0) then
      call init(matsB,N,L)
    else
      call hubbard(matsB,ns,N,L,beta,t,sigma,U,myh((iter-1)*L*N+1:iter*L*N)) 
    endif

    ! Convert the dense format of M with Bi's and I into CSR format.
    ! Reverse this format for other algorithms and applications.
    call dense2csr(vals,ci,rp,matsB,N,L,nnz)

    ! Run the 3-step FSI:
    ! Step(1). block cyclic reduction
    call reduction(matsB,matsBh,work,N,L,b,q)
    
    ! Step(2). bsofi
    call fullpcyc(matMh,matsBh,N,b)
    call bsofi(N,b,matMh)
    ! Alternatively we can call inversion as:
    ! call inversion(matMh,bN)
    ! but bsofi is more efficient.

    ! Step(3). wrapping
    call seedscopy(matsG,matMh,N,L,b,q)
    call wrapping(matsG,matsB,N,L,b,q,quantity)
    
    lquantity = lquantity + quantity
  enddo

#ifdef MPI_RUN
   call MPI_Reduce(lquantity, gquantity, 1,      &
                  MPI_DOUBLE_PRECISION, MPI_SUM, &
                  my_root, MPI_COMM_WORLD,ierr )
   call MPI_FINALIZE(ierr)
#endif

  if(myid == my_root ) then
    call system_clock(t_total_end,rate)
    write (*,'(i4,A,i4,A)') b, ' out of ', L, ' block columns of G have been computed.'
    write (*,*) 'total time', (t_total_end - t_total_start) /real(rate)
  endif

  deallocate(matsB,matsBh,work,matMh,matsG,myh)
end program runFSI
