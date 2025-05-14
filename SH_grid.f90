program SH2grid

    use SHTOOLS

    implicit none
    integer,parameter :: lm = 20
    integer,parameter :: lmin = 1
    real*8,dimension(:,:),allocatable :: griddh
    real*8,dimension(:,:,:),allocatable :: cilm1,cilm2,cilm3
    integer*4 n,lmax,norm,sampling,csphase,lmax_calc,extend,exitstatus
    integer*4 i,j,l,m,skip
    character*100 header,filename1,filename2,filename3

    real*8,dimension(:,:),allocatable :: rad,theta,phi,total,pot_grid


    real*8,dimension(3) :: value1,value2,value3
    real*8 a,r,lat,lon,f,r0

!-----------------------read the SH file    
    filename1 = 'SH_220.txt'
    filename2 = 'SH_221.txt'
    filename3 = 'SH_219.txt'
    skip = 0
    allocate(cilm1(2,lm+1,lm+1))
    allocate(cilm2(2,lm+1,lm+1))
    allocate(cilm3(2,lm+1,lm+1))
    call SHRead(filename1,cilm1,lmax,skip)
    call SHRead(filename2,cilm2,lmax,skip)
    call SHRead(filename3,cilm3,lmax,skip)


!------------------------------------------------

!-----------------------SH2Grid     time:2019
!-----------------------parameters
    a = 6371200
    r = 3480000
    r0 = 6371200
    n = 2*lmax+2
    allocate(rad(n,2*n))
    allocate(theta(n,2*n))
    allocate(phi(n,2*n))
    allocate(total(n,2*n))
    allocate(pot_grid(n,2*n))
    sampling = 2

!------------------------Br_CMB

    open(400,file='Br_CMB.txt')
    do lon = -180,180,0.5
        do lat = 89.5,-89.5,-0.5
            lmax = 13
            value1 = MakeMagGridPoint(cilm1,lmax,a,r,lat,lon)
            write(400,*) lon,lat,value1(1)
            !print *, lon,lat,value(1)
        end do
    end do
    close(400)

!--------------------------Br_CMB_SV

    open(500,file='Br_CMB_SV.txt')
    do lon = -180,180,0.5
        do lat = 89.5,-89.5,-0.5
            lmax = 17
            value1 = MakeMagGridPoint(cilm1,lmax,a,r,lat,lon)
            value2 = MakeMagGridPoint(cilm2,lmax,a,r,lat,lon)
            write(500,*) lon,lat,(value2(1)-value1(1))/0.09993
            !print *, lon,lat,value(1)
        end do
    end do
    close(500)


 !--------------------------Br_CMB_SA
    open(900,file='Br_CMB_SA.txt')
    do lon = -180,180,0.5
        do lat = 89.5,-89.5,-0.5
            lmax = 15
            value1 = MakeMagGridPoint(cilm1,lmax,a,r,lat,lon)
            value2 = MakeMagGridPoint(cilm2,lmax,a,r,lat,lon)
            value3 = MakeMagGridPoint(cilm3,lmax,a,r,lat,lon)
            write(900,*) lon,lat,((value2(1)-value1(1))/0.09993-(value1(1)-value3(1))/0.09993)/0.09993
            !print *, lon,lat,value(1)
        end do
    end do
    close(900)    
   
!--------------------------Br_surface_SA
    open(600,file='Br_surface_SA.txt')
    do lon = -180,180,0.5
        do lat = 89.5,-89.5,-0.5
            lmax = 20
            value1 = MakeMagGridPoint(cilm1,lmax,a,r0,lat,lon)
            value2 = MakeMagGridPoint(cilm2,lmax,a,r0,lat,lon)
            value3 = MakeMagGridPoint(cilm3,lmax,a,r0,lat,lon)
            write(600,*) lon,lat,((value2(1)-value1(1))/0.09993-(value1(1)-value3(1))/0.09993)/0.09993
            !print *, lon,lat,value(1)
        end do
    end do
    close(600)    

!---------------------------Br_surface_SV
    open(700,file='Br_surface_SV.txt')
    do lon = -180,180,0.5
        do lat = 89.5,-89.5,-0.5
            lmax = 20
            value1 = MakeMagGridPoint(cilm1,lmax,a,r0,lat,lon)
            value2 = MakeMagGridPoint(cilm2,lmax,a,r0,lat,lon)
            write(700,*) lon,lat,(value2(1)-value1(1))/0.09993
            !print *, lon,lat,value(1)
        end do
    end do
    close(700)   




!---------------------------Br_surface

    open(800,file='Br_surface.txt')
    do lon = -180,180,0.5
        do lat = 89.5,-89.5,-0.5
            lmax = 13
            value1 = MakeMagGridPoint(cilm1,lmax,a,r0,lat,lon)
            write(800,*) lon,lat,value1(1)
            !print *, lon,lat,value(1)
        end do
    end do
    close(800)


    

    deallocate(rad)
    deallocate(theta)
    deallocate(phi)
    deallocate(total)
    deallocate(pot_grid)

  

    
end program SH2grid