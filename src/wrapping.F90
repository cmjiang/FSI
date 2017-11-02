subroutine wrapping(matsG,matsB,N,L,b,q,quantity)
  implicit none
  integer, intent(in) :: N, L, b, q
  double precision, intent(in) :: matsB(N,N,L)
  double precision, intent(out) :: quantity
  double precision, intent(inout) :: matsG(N,N,L,b)
  integer :: thn, OMP_GET_THREAD_NUM 
  double precision :: Gnorm(L,b)
  integer :: c, k0, l0, kk, kn, ll, i, j ,  ii, jj
  c = L/b
  Gnorm = 0
!$omp parallel private(k0, l0, kk, ll, kn, thn)
!$omp do
  do j = 1, b*b
    k0 = mod(j,b) + 1
    l0 = (j-1)/b + 1
    kk = c*k0 - q
    ll = c*l0 - q
    do i = 1, (c-1)/2
      kn = kk - 1
      if ( kn == 0 ) then
        kn = L
      endif
      call wrapup(matsG(:,:,kn,l0), matsG(:,:,kk,l0), matsB, N, L, kk, ll)
      kk = kn
      do ii = 1, N
        do jj = 1, N
          Gnorm(kn,l0) = Gnorm(kn,l0) + matsG(ii,jj,kn,l0)*matsG(ii,jj,kn,l0)
        enddo
      enddo
    enddo
    kk = c*k0 - q
    ll = c*l0 - q
    do i = 1, c/2
      kn = kk + 1
      if ( kn == L+1 ) then
        kn = 1
      endif
      call wrapdown(matsG(:,:,kn,l0), matsG(:,:,kk,l0), matsB, N, L, kk, ll)
      kk = kn
      do ii = 1, N
        do jj = 1, N
          Gnorm(kn,l0) = Gnorm(kn,l0) + matsG(ii,jj,kn,l0)*matsG(ii,jj,kn,l0)
        enddo
      enddo
    enddo
  enddo
!$omp end do
!$omp end parallel
  do i = 1, L
    do j = 1, b
      quantity = quantity + Gnorm(i,j)
    enddo
  enddo
  quantity = quantity / L / b
end subroutine wrapping

subroutine wrapup(Go, Gi, B, N, L, i0 ,j0)
  implicit none
  double precision, intent(out) :: Go(N,N)
  double precision, intent(in) :: Gi(N,N), B(N,N,L)
  integer, intent(in) :: N, L, i0, j0
  integer :: i, k
  double precision, allocatable :: Bi(:,:), Gt(:,:)
  allocate(Bi(N,N),Gt(N,N))
  i = i0
  Go = Gi
  if ( i == j0 ) then
    forall( k = 1:N ) Go(k,k) = Go(k,k) - 1
  endif
  Bi = B(:,:,i)
  call inversion(Bi,N)
  call mmultiply(Gt,Bi,Go,N)
  Go = Gt
 
  if ( i == 1 ) then
    Go = -Go
  endif
  deallocate(Bi,Gt)
end subroutine wrapup

subroutine wrapleft(Go, Gi, B, N, L, i0 ,j0)
  implicit none
  double precision, intent(out) :: Go(N,N)
  double precision, intent(in) :: Gi(N,N), B(N,N,L)
  integer, intent(in) :: N, L, i0, j0
  integer :: j, k
  j = j0
  Go = Gi
  call mmultiply(Go,Gi,B(:,:,j),N)
  if( j == 1 ) then
    Go = -Go
  endif
  if ( i0 == j-1 ) then
    forall( k = 1:N ) Go(k,k) = Go(k,k) + 1
  endif
end subroutine wrapleft

subroutine wrapdown(Go, Gi, B, N, L, i0 ,j0)
  implicit none
  double precision, intent(out) :: Go(N,N)
  double precision, intent(in) :: Gi(N,N), B(N,N,L)
  integer, intent(in) :: N, L, i0, j0
  integer :: i, k
  i = i0 + 1
  if ( i > L ) then
    i = 1
  endif
  call mmultiply(Go,B(:,:,i),Gi,N)
  if ( i == 1 ) then
    Go = -Go
  endif
  if ( i == j0 ) then
    forall( k = 1:N ) Go(k,k) = Go(k,k) + 1
  endif
end subroutine wrapdown

subroutine wrapright(Go, Gi, B, N, L, i0 ,j0)
  implicit none
  double precision, intent(out) :: Go(N,N)
  double precision, intent(in) :: Gi(N,N), B(N,N,L)
  integer, intent(in) :: N, L, i0, j0
  integer :: j, k
  double precision, allocatable :: Bi(:,:), Gt(:,:)
  allocate(Bi(N,N), Gt(N,N))
  Go = Gi
  j = j0 + 1
  if ( j > L ) then
    j = 1
  endif
  if ( j0 == i0 ) then
    forall( k = 1:N ) Go(k,k) = Go(k,k) - 1
  endif
  Bi = B(:,:,j)
  call inversion(Bi,N)
  call mmultiply(Gt,Go,Bi,N)
  Go = Gt
  if( j ==1 ) then
    Go = -Go
  endif
  deallocate(Bi,Gt)
end subroutine wrapright
