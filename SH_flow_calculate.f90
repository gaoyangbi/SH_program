program SH2grid

    use SHTOOLS

    implicit none
    integer,parameter :: lm = 20
    integer,parameter :: lmin = 1
    real*8,parameter :: pi = 4*atan(1.0)
    real*8,dimension(:,:,:),allocatable :: cilm1,cilm2
    integer*4 lmax,norm,sampling,csphase,lmax_calc,extend,exitstatus
    integer*4 i,j,skip
    character*100 header,filename1,filename2


    real*8,dimension(3) :: value1,value2,value3,value4
    real*8 a,r,lat,lon,f,r0,Br_t,Br_x
    real*8 m,n

!-----------------------read the SH file    
    filename1 = 'SH_30.txt'
    filename2 = 'SH_31.txt'
    skip = 0
    allocate(cilm1(2,lm+1,lm+1))
    allocate(cilm2(2,lm+1,lm+1))
    call SHRead(filename1,cilm1,lmax,skip)
    call SHRead(filename2,cilm2,lmax,skip)


!------------------------------------------------

!-----------------------SH2Grid
!-----------------------parameters
    a = 6371200
    r = 3480000
    r0 = 6371200
    n = 2*lmax+2
    sampling = 2

!-------------------------calculate flow field\  赤道线上的流场速度
    do lon = -179.75,179.75,0.5
        lmax = 20
        lat = 0.0
        value1 = MakeMagGridPoint(cilm1,lmax,a,r,lat,lon) 
        value2 = MakeMagGridPoint(cilm1,lmax,a,r,lat,lon+0.5)
        value3 = MakeMagGridPoint(cilm1,lmax,a,r,lat,lon+0.25)
        value4 = MakeMagGridPoint(cilm2,lmax,a,r,lat,lon+0.25)
        Br_t = (value4(1)-value3(1))/0.10047
        Br_x = (value2(1)-value1(1))/0.5
        m = (-Br_t) / (Br_x)
        print *, m*pi/180.0*r
    end do

    


    deallocate(cilm1)
    deallocate(cilm2) 

    
end program SH2grid

