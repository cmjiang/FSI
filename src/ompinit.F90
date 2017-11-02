subroutine ompinit()
  implicit none
  integer                 :: nproc    
  character(len=32)       :: argv
  integer                 :: nprocmax, OMP_GET_NUM_PROCS, argc

  argc = iargc()
  nprocmax = OMP_GET_NUM_PROCS()
  call getarg(2, argv)
  if(argv == '-p') then
     if(argc == 2) then
        nproc =  nprocmax
        call OMP_SET_NUM_THREADS(nproc)
        print *, 'number of processors:', nproc
     else
        call getarg(3,argv)
        read(argv,'(I10)') nproc
        if(nproc > 0 .AND. nproc .le. nprocmax) then
           call OMP_SET_NUM_THREADS(nproc)
           print *, nproc
        else
           write(*, "('Invalid number of threads, must be an integer &
              between 1 and the max number of threads (', i2, ').')") nprocmax
           stop
        endif
     endif
  else
    call OMP_SET_NUM_THREADS(nprocmax)
  endif
end subroutine ompinit
