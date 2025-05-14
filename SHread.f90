    program SHreadchaos

        use SHTOOLS

        implicit none
        integer,parameter :: lmax = 20
        integer,parameter :: lmin = 1
        integer,parameter :: epoch = 168
        integer,parameter :: epoch_index = 219
        real*8,dimension(:,:),allocatable :: Clm_array,time_epoch,Slm_array,read_array
        integer*4 lsum
        integer*4 i,j,l,m
        character*100 header
        character*100 filename


        lsum = (lmax+1)*(lmax+1)-1
        allocate(Clm_array(1,epoch))
        allocate(Slm_array(1,epoch))
        allocate(read_array(1,epoch))
        allocate(time_epoch(1,epoch))



!----------------------------------------------------- 
        open(100,file = 'CHAOS-7.12_core.shc')
        open(200, file = 'SH_219.txt')

        do i=1,4,1
            read(100,*) header
        end do

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
!---------------------------------------------------------   
        

        deallocate(Clm_array)
        deallocate(Slm_array)
        deallocate(time_epoch)
        deallocate(read_array)
!-------------------------------------------------test
!        filename = 'SH.txt'
!       call SHRead(filename,clim,lm,skip)
!        open(300,file='see.txt')
!        do i=1,21
!            do j=1,21
!                write(300,*) i-1,j-1,clim(1,i,j),clim(2,i,j)
!            end do
!        end do
!        close(300)
!        print *,clim(1,8,2)

        
    end program SHreadchaos