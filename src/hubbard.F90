subroutine hubbard(Bl, ns, N, L, beta, t, sigma, U, h)
  implicit none
  integer, intent(in)  :: ns, N, L, t, sigma, U, h(L*N)
  double precision, intent(out) :: Bl(N,N,L)
  double precision, intent(in) :: beta
  double precision :: deltau, v, optwork(1)
  double precision, allocatable :: Vl(:,:,:), Dl(:,:,:), Kx(:,:), B(:,:)
  double precision, allocatable :: K1(:,:), K2(:,:), eye(:,:), W(:), work(:)
  integer :: i, j, lwork, info
  allocate(Vl(N,N,L),Dl(N,N,L),Kx(ns,ns),eye(ns,ns),K1(N,N),K2(N,N))
  allocate(W(N),B(N,N))
  B = 0
  Bl = 0
  Vl = 0
  K1 = 0
  K2 = 0
  eye = 0
  deltau = beta/L
  v = sqrt(deltau*U) + ((deltau*U)**(3/2))/12
  do i = 1, L
    forall(j = 1:N) Vl(j,j,i) = h((i-1)*N+j)
  enddo
  do i = 1, L
    do j = 1, N
      Dl(j,j,i) = exp(sigma*v*Vl(j,j,i))  
    enddo
  enddo
  do i = 1, ns-1
    Kx(i,i+1) = 1
    Kx(i+1,i) = 1
  enddo
  Kx(1,ns) = 1
  Kx(ns,1) = 1
  forall(i = 1:ns) eye(i,i) = 1
  call nkronecker(K1,eye,Kx,ns)
  call nkronecker(K2,Kx,eye,ns)
  K1 = K1 + K2
  K1 = t*deltau*K1
  call dsyev('V','L',N,K1,N,W,optwork,-1,info)
  lwork = optwork(1)
  allocate(work(lwork))
  call dsyev('V','L',N,K1,N,W,work,lwork,info)
  forall(i = 1:N) B(i,i) = exp(W(i))
  call dgemm('N','N',N,N,N,1.0d0,K1,N,B,N,0.0d0,K2,N)
  call dgemm('N','T',N,N,N,1.0d0,K2,N,K1,N,0.0d0,B,N)
  do i = 1,L
    call dgemm('N','N',N,N,N,1.0d0,B,N,Vl(:,:,i),N,0.0d0,Bl(:,:,i),N)
  enddo
  deallocate(Vl,Dl,Kx,eye,K1,K2,W,B)
end subroutine hubbard

subroutine randomh(h, N, L, m)
  implicit none
  integer, intent(in) :: N, L, m
  integer, intent(out) :: h(L*N*m)
  double precision :: r(L*N*m)
  integer :: i, j, k
  call random_number(r)
  do i = 1, L*N*m
    if( r(i) > 0.5 ) then
      h(i) = 1
    else
      h(i) = -1
    endif
  enddo
end subroutine randomh

subroutine nkronecker(A,B,C,n)
  implicit none
  integer, intent(in) :: n
  double precision, intent(in) :: B(n,n),C(n,n)
  double precision, intent(out) :: A(n*n,n*n)
  integer :: i, j
  do i = 1, n
    do j = 1, n
      A((i-1)*n+1:i*n,(j-1)*n+1:j*n)=B(i,j)*C
    enddo
  enddo
end subroutine nkronecker
