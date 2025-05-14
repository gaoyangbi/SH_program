program SHread2grid

    use SHTOOLS

    implicit none
    integer,parameter :: lmax = 18
    integer,parameter :: lmin = 1

    integer,parameter :: epoch = 84    !168

    
    integer,parameter :: skip = 0   ! skip lines of txt
    real*8,dimension(:,:),allocatable :: Clm_array,time_epoch,Slm_array,read_array
    real*8,dimension(:,:),allocatable :: griddh
    real*8,dimension(:,:,:),allocatable :: cilm1,cilm2,cilm3
    integer*4 lsum,epoch_index,lm
    integer*4 i,j,l,m
    character*100 header
    character*100 coef_path,grid_path,txt_name
    character*100 filename,filenum

    real*8,dimension(3) :: value1
    real*8,parameter :: a = 6371200   ! 参考半径
    real*8,parameter :: r0 = 6371200
    real*8 lat,lon


    txt_name  = 'coefficients_2014_2020_3.txt'
    coef_path = '/data2/ztwang/project/zyg/SH_program/2014_2020_pacific_result_coef/'
    grid_path = '/data2/ztwang/project/zyg/SH_program/2014_2020_pacific_result_grid/'


    lsum = (lmax+1)*(lmax+1)-1
    allocate(Clm_array(1,epoch))
    allocate(Slm_array(1,epoch))
    allocate(read_array(1,epoch))
    allocate(time_epoch(1,epoch))

    allocate(cilm1(2,lmax+1,lmax+1))



!----------------------------------------------------- read SH2coef

    do epoch_index=1,epoch,1
        write(filenum,'(I3.3)') epoch_index
        filename = trim(coef_path)//'SH_coef_'//trim(filenum)//'.txt'
        open(200, file = trim(filename))

        open(100, file = trim(txt_name))   
        if (skip .gt. 0) then
            do i=1,skip,1
                read(100,*) header
            end do
        end if
        read(100,*) (time_epoch(1,j) , j = 1 , epoch)

        do i=1,lsum,1
            read(100,*) l,m,(read_array(1,j) , j = 1 , epoch)
            if (m == 0) then
                write(200,*) l,m,read_array(1,epoch_index)
            else if (m > 0) then
                Clm_array = read_array
            else if (m < 0) then
                Slm_array = read_array
                write(200,*) l,-m,Clm_array(1,epoch_index),Slm_array(1,epoch_index)
            end if
            
            !print *, Clm_array(i,2),Clm_array(i,epoch+2)
        end do                
        close(100)
        close(200)
    end do

    
!---------------------------------------------------------   

    deallocate(Clm_array)
    deallocate(Slm_array)
    deallocate(time_epoch)
    deallocate(read_array)
!----------------------------------------------SH2grid
    do epoch_index=1,epoch,1

    !-----------------------read the SH file  
        write(filenum,'(I3.3)') epoch_index
        filename = trim(coef_path)//trim('SH_coef_')//trim(filenum)//'.txt' 

        call SHRead(filename,cilm1,lm,0)           
    
    !---------------------------Br_surface
        
        filename = trim(grid_path)//trim('Br_surface_')//trim(filenum)//'.txt' 
        open(800,file=trim(filename))
        do lon = -180,180,0.5
            do lat = 89.5,-89.5,-0.5
                value1 = MakeMagGridPoint(cilm1,lmax,a,r0,lat,lon)
                write(800,*) lon,lat,value1(1)
                !print *, lon,lat,value(1)
            end do
        end do
        close(800)
     
    
    
    end do
    
    deallocate(cilm1)


    
end program SHread2grid