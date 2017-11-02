subroutine init(matA, N, L)
  implicit none
  integer, intent(in)  :: N, L
  double precision, intent(inout), target :: matA(N,N,L) 
  integer :: LN, k
  integer :: i, j
  double precision, pointer :: value(:,:)
  LN = N*L
  do k = 1, L
    value => matA(:,:,k)
    call random_number(value(1:N,1:N))
  enddo
end subroutine init

subroutine fullpcyc(matM, matB, N, L)
  implicit none
  integer, intent(in) :: N, L
  double precision, intent(inout) :: matM(N*L,N*L)
  double precision, intent(in) :: matB(N,N,L)
  integer :: LN, i
  LN = N*L
  matM = 0
  forall(i = 1:LN) matM(i,i) = 1
  matM(1:N,(L-1)*N+1:LN) = matB(:,:,1)
  do i = 1, L-1
    matM(1+i*N:N+i*N,1-N+i*N:i*N) = - matB(:,:,i+1)
  enddo
end subroutine fullpcyc

subroutine formsketch(matS,matsG,N,L,b,q)
  implicit none
  integer, intent(in) :: N,L,b,q
  double precision, intent(in) :: matsG(N,N,L,b)
  double precision, intent(out) :: matS(N*L,N*b)
  integer :: i,j,c,k
  c = L/b
  do i = 1,L
    do j = 1,b
      k = c*j - q
      matS(1+(i-1)*N:i*N, 1+(k-1)*N:k*N) = matsG(1:N,1:N,i,j)
    enddo
  enddo
end subroutine formsketch

subroutine seedscopy(matsG,matMh,N,L,b,q)
  implicit none
  integer, intent(in) :: N, L, b, q
  double precision, intent(in) :: matMh(b*N,b*N)
  double precision, intent(out) :: matsG(N,N,L,b)
  integer k0, l0, kk, ll, c
  c = L/b
  do k0 = 1,b
    do l0 = 1,b
      kk = c*k0 - q
      ll = l0
      matsG(:,:,kk,ll) = matMh(1+(k0-1)*N:k0*N, 1+(l0-1)*N:l0*N)
    enddo
  enddo
end subroutine seedscopy

subroutine normdiff(eps,matM,matsG,N,L,b,q)
  implicit none
  integer, intent(in) :: N,L,b,q
  double precision, intent(out) :: eps
  double precision, intent(in) :: matM(N*L,N*L), matsG(N,N,L,b)
  integer :: i,j,c,k
  double precision :: dlange,work,matD(N,N)
  c = L/b
  eps = 0
  do i = 1, L
    do j = 1, b
      k = c*j - q
      matD = matsG(1:N,1:N,i,j) - matM(1+(i-1)*N:i*N, 1+(k-1)*N:k*N)
      eps = eps + dlange('F',N,N,matD,N,work)/ &
            dlange('F',N,N,matsG(1:N,1:N,i,j),N,work)
      print *, eps
    enddo
  enddo
  eps = eps/L/b
end subroutine normdiff

subroutine inversion(M,n)
  implicit none
  double precision, intent(inout) :: M(n,n)
  integer, intent(in) :: n
  integer :: ipiv(n), info, lwork
  double precision, allocatable :: work(:)
  double precision :: optwork(1)
  call dgetrf(n,n,M,n,ipiv,info)
  call dgetri(n,M,n,ipiv,optwork(1),-1,info)
  lwork = optwork(1)
  allocate(work(lwork))
  call dgetri(n,M,n,ipiv,work,lwork,info)
  deallocate(work)
end subroutine inversion

subroutine mmultiply(Mo,Ml,Mr,n)
  implicit none
  double precision, intent(inout) :: Mo(n,n)
  double precision, intent(in) :: Ml(n,n), Mr(n,n)
  integer, intent(in) :: n
  call dgemm('N','N',n,n,n,1.0d0,Ml,n,Mr,n,0.0d0,Mo,n)
end subroutine mmultiply

subroutine printmatarr(matA, N, L)
  implicit none
  integer, intent(in) :: N, L
  double precision, intent(inout) :: matA(N,N,L) 
  integer :: LN, i, j, k
  LN = N*L
  do k = 1, L
    print *, '-----matrix begin-----'
    do i = 1, N
      do j = 1, N
        write(*,'(f6.2)',advance='no') matA(i,j,k)
      enddo
      write(*,*)
    enddo 
    print *, '------ matrix end ------'
  enddo
end subroutine printmatarr

subroutine printmat(matM, m, n)
  implicit none
  integer, intent(in) :: m, n
  double precision, intent(in) :: matM(m,n)
  integer :: i, j
  print *, '----- matrix begin -----'
  do i = 1, m
    do j = 1, n
      write(*,'(f6.2)',advance='no') matM(i,j)
    enddo
    write(*,*)
  enddo
  print *, '------ matrix end ------'
end subroutine printmat
