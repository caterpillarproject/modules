subroutine cic(pos,mpart,boxsize,dim,redshift,npart)
use omp_lib
implicit none 

real*8, intent(in) :: mpart
!f2py intent(in) :: mpart

real*8, intent(in) :: boxsize
!f2py intent(in) :: boxsize

integer, intent(in) :: dim
!f2py intent(in) :: dim

real*4, intent(in) :: redshift
!f2py intent(in) :: redshift

integer, intent(in) :: npart
!f2py intent(hide) :: npart

real*8, intent(in), dimension(0:npart-1,0:2) :: pos
!f2py intent(in) :: pos

integer   i
integer   Ntot
integer   :: xi,yi,zi,xii,yii,zii
real*8, dimension(255,255,255) :: mesh_out
real*8,     dimension(:),   allocatable    :: dx1, dy1, dz1, dx2, dy2, dz2 
real*8,     dimension(:),   allocatable    :: x0, y0, z0
real*8,     dimension(:),   allocatable    :: x1, y1, z1, x2, y2, z2

real*8     ndim,seconds
integer :: numthreads
real*8             :: vol_cell, vol_box, co_ph, mu
character(len=500) :: fname

!real*8,     dimension(:,:,:),   allocatable :: mesh_out


numthreads = 12

call omp_set_num_threads(numthreads)  

ndim = REAL(dim)
Ntot = npart - 1
write(*,*)"Allocating for npart,Ntot",npart,Ntot,huge(npart)
write(*,*)"Dim:",dim,ndim
write(*,*)"Boxsize:",boxsize
write(*,*)"Redshift:",redshift

do i = 1,5
    write(*,*)pos(i,0),pos(i,1),pos(i,2)
end do

!###########################################
!     USING POSITIONS - CREATE mesh_out XYZ
!###########################################
seconds = omp_get_wtime()
allocate(x0(Ntot),y0(Ntot),z0(Ntot))
x0(:) = pos(0:Ntot,0)
y0(:) = pos(0:Ntot,1)
z0(:) = pos(0:Ntot,2)

allocate(x1(Ntot))
x1(:) = floor((x0(:))*((ndim)/boxsize) - 0.5) + 1
allocate(dx1(Ntot))
dx1(:) = 1 + x1 - 0.5 + (((ndim)/boxsize)*(-x0))
deallocate(x0)
allocate(dx2(Ntot))
dx2(:) = 1 - dx1(:)

allocate(y1(Ntot))
y1(:) = floor((y0(:))*((ndim)/boxsize) - 0.5) + 1
allocate(dy1(Ntot))
dy1(:) = 1 + y1 - 0.5 + (((ndim)/boxsize)*(-y0))
deallocate(y0)
allocate(dy2(Ntot))
dy2(:) = 1 - dy1(:)

allocate(z1(Ntot))
z1(:) = floor((z0(:))*((ndim)/boxsize) - 0.5) + 1
allocate(dz1(Ntot))
dz1(:) = 1 + z1 - 0.5 + (((ndim)/boxsize)*(-z0))
deallocate(z0)
allocate(dz2(Ntot))
dz2(:) = 1 - dz1(:)

allocate(x2(Ntot), y2(Ntot), z2(Ntot))
x2(:) = x1(:) + 1
y2(:) = y1(:) + 1
z2(:) = z1(:) + 1

deallocate(dx1, dy1, dz1, dx2, dy2, dz2, x1, y1, z1, x2, y2, z2)

seconds = omp_get_wtime() - seconds
write(*,"(a,F8.3)")" Time [s] to allocate vecs: ",seconds

!###########################################
!  LOOP OVER ALL PARTICLES & POPULATE MESH
!###########################################
seconds = omp_get_wtime()
!allocate(mesh_out(0:dim-1,0:dim-1,0:dim-1))
!allocate(mesh_out(Ntot,Ntot,Ntot))
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
      mesh_out(xi,yi,zi)   = mesh_out(xi,yi,zi)   + dx1(i) * dy1(i) * dz1(i) * mpart*1.98892D33 *0.18
      mesh_out(xii,yi,zi)  = mesh_out(xii,yi,zi)  + dx2(i) * dy1(i) * dz1(i) * mpart*1.98892D33 *0.18
      mesh_out(xi,yii,zi)  = mesh_out(xi,yii,zi)  + dx1(i) * dy2(i) * dz1(i) * mpart*1.98892D33 *0.18
      mesh_out(xii,yii,zi) = mesh_out(xii,yii,zi) + dx2(i) * dy2(i) * dz1(i) * mpart*1.98892D33 *0.18
      mesh_out(xi,yi,zii)  = mesh_out(xi,yi,zii)  + dx1(i) * dy1(i) * dz2(i) * mpart*1.98892D33 *0.18
      mesh_out(xii,yi,zii) = mesh_out(xii,yi,zii) + dx2(i) * dy1(i) * dz2(i) * mpart*1.98892D33 *0.18
      mesh_out(xi,yii,zii) = mesh_out(xi,yii,zii) + dx1(i) * dy2(i) * dz2(i) * mpart*1.98892D33 *0.18
      mesh_out(xii,yii,zii)= mesh_out(xii,yii,zii)+ dx2(i) * dy2(i) * dz2(i) * mpart*1.98892D33 *0.18
  end if
end do
!$OMP END PARALLEL DO
        
vol_cell = ((boxsize*3.08568025D24)/dim/0.6711)**3
vol_box = (boxsize*3.08568025D24/0.6711)**3

co_ph = 1 + redshift
mu = 1/(1.22*1.67262D-24)

write(*,*)"cell volume [cm^3]:",vol_cell
write(*,*)"box volume [cm^3]:",vol_box

write(*,*)"average mass density [g cm^-3]: ",SUM(mesh_out(:,:,:)*co_ph**3)/vol_box
write(*,*)"average number density [cm^-3]: ",SUM(mesh_out(:,:,:)*mu*co_ph**3)/vol_box

!write(9,*)redshift, SUM(mesh_out(:,:,:)*co_ph**3)/vol_box, SUM(mesh_out(:,:,:)*mu*co_ph**3)/vol_box

!write(fname,'(a,i3,a,f5.3,a)')"/bigbang/data/bgriffen/c2ray/cicfiles/parent/",dim,"/density/z",redshift,"_new.dat"
!write(*,*)"writing output to file: ",trim(fname)

!mesh_out = (mesh_out(:,:,:))/vol_cell

!open(unit=8,file=fname,form="unformatted")
!write(8)mesh_out
!close(8)

seconds = omp_get_wtime() - seconds
write(*,"(a,F8.3)")" Time [s] to populate mesh: ",seconds
deallocate(dx1, dy1, dz1, dx2, dy2, dz2, x1, y1, z1, x2, y2, z2)

end subroutine
