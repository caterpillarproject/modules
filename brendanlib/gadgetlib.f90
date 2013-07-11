subroutine cic(pos,mpart,boxsize,dim,npart,mesh_out)
use omp_lib
implicit none 

real*8, intent(in) :: mpart
!f2py intent(in) :: mpart

real*8, intent(in) :: boxsize
!f2py intent(in) :: boxsize

integer, intent(in) :: dim
!f2py intent(in) :: dim

integer, intent(in) :: npart
!f2py intent(hide) :: npart

real*8, intent(in), dimension(0:npart-1,0:2) :: pos
!f2py intent(in) :: pos

real*8, intent(out),  dimension(0:dim-1,0:dim-1,0:dim-1):: mesh_out
!f2py intent(out) :: mesh_out

integer*4   i
integer*8   Ntot
integer   :: xi,yi,zi,xii,yii,zii
real*4,     dimension(:),   allocatable    :: dx1, dy1, dz1, dx2, dy2, dz2 
real*4,     dimension(:),   allocatable    :: x0, y0, z0
real*4,     dimension(:),   allocatable    :: x1, y1, z1, x2, y2, z2
real*8     ndim, seconds
integer :: id
integer :: numthreads

numthreads = 8

call omp_set_num_threads(numthreads)  

ndim = REAL(dim)
Ntot = npart - 1

!###########################################
!     USING POSITIONS - CREATE mesh_out XYZ
!###########################################
seconds = omp_get_wtime()
allocate(x0(Ntot),y0(Ntot),z0(Ntot))
x0(:) = pos(0:Ntot,0)
y0(:) = pos(0:Ntot,1)
z0(:) = pos(0:Ntot,2)

allocate(x1(Ntot))
x1(:) = floor((x0(:))*((ndim-2)/boxsize) - 0.5) + 1
allocate(dx1(Ntot))
dx1(:) = 1 + x1 - 0.5 + (((ndim-2)/boxsize)*(-x0))
deallocate(x0)
allocate(dx2(Ntot))
dx2(:) = 1 - dx1(:)

allocate(y1(Ntot))
y1(:) = floor((y0(:))*((ndim-2)/boxsize) - 0.5) + 1
allocate(dy1(Ntot))
dy1(:) = 1 + y1 - 0.5 + (((ndim-2)/boxsize)*(-y0))
deallocate(y0)
allocate(dy2(Ntot))
dy2(:) = 1 - dy1(:)

allocate(z1(Ntot))
z1(:) = floor((z0(:))*((ndim-2)/boxsize) - 0.5) + 1
allocate(dz1(Ntot))
dz1(:) = 1 + z1 - 0.5 + (((ndim-2)/boxsize)*(-z0))
deallocate(z0)
allocate(dz2(Ntot))
dz2(:) = 1 - dz1(:)

allocate(x2(Ntot), y2(Ntot), z2(Ntot))
x2(:) = x1(:) + 1
y2(:) = y1(:) + 1
z2(:) = z1(:) + 1

seconds = omp_get_wtime() - seconds
write(*,"(a,F8.3)")" Time [s] to allocate vecs: ",seconds

!###########################################
!  LOOP OVER ALL PARTICLES & POPULATE MESH
!###########################################
seconds = omp_get_wtime()
mesh_out(:,:,:) = 0.0
!$OMP PARALLEL DO PRIVATE(xi,yi,zi,xii,yii,zii)
do i = 1,Ntot

  xi = int(x1(i))
  yi = int(y1(i))
  zi = int(z1(i))
  if (xi.ne.0.0) then
      xii = int(x2(i))
      yii = int(y2(i))
      zii = int(z2(i))
      mesh_out(xi,yi,zi)   = mesh_out(xi,yi,zi)   + dx1(i) * dy1(i) * dz1(i) * mpart
      mesh_out(xii,yi,zi)  = mesh_out(xii,yi,zi)  + dx2(i) * dy1(i) * dz1(i) * mpart
      mesh_out(xi,yii,zi)  = mesh_out(xi,yii,zi)  + dx1(i) * dy2(i) * dz1(i) * mpart
      mesh_out(xii,yii,zi) = mesh_out(xii,yii,zi) + dx2(i) * dy2(i) * dz1(i) * mpart
      mesh_out(xi,yi,zii)  = mesh_out(xi,yi,zii)  + dx1(i) * dy1(i) * dz2(i) * mpart
      mesh_out(xii,yi,zii) = mesh_out(xii,yi,zii) + dx2(i) * dy1(i) * dz2(i) * mpart
      mesh_out(xi,yii,zii) = mesh_out(xi,yii,zii) + dx1(i) * dy2(i) * dz2(i) * mpart
      mesh_out(xii,yii,zii)= mesh_out(xii,yii,zii)+ dx2(i) * dy2(i) * dz2(i) * mpart
  end if
end do
!$OMP END PARALLEL DO
        
!c  !$OMP PARALLEL DO PRIVATE(xi,yi,zi,xii,yii,zii)

!c  !$OMP END PARALLEL DO


seconds = omp_get_wtime() - seconds
write(*,"(a,F8.3)")" Time [s] to populate mesh: ",seconds
deallocate(dx1, dy1, dz1, dx2, dy2, dz2, x1, y1, z1, x2, y2, z2)
end subroutine