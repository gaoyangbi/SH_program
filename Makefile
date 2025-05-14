FC = ifort
OMPFC = mpif90
MPIFC = mpiifort
LIBPATH = -L/data2/ztwang/project/zyg/SHtools/SHTOOLS-master/lib 
INCLUDE = -I/data2/ztwang/project/zyg/SHtools/SHTOOLS-master/include
SCALAPACK = -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64
LIBS = -lSHTOOLS -lfftw3 -lm -qmkl -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread  -lmkl_core -lpthread -liomp5
#COPTS =   -qopenmp -Wl,-stack_size,0x40000000,-stack_addr,0xf0000000
COPTS =  -fopenmp
#all: shale01 shale02 shale03 shale04 shale05 shale06 shalecg02 shalecg03 shalecg04 shalecg05 shalecg06 shalecg07 shalecg08
all: SHread2grid readspline_zyg
#ifeq ($(FC),gfortran)
#	COPTS = -fopenmp
#	LIBS = -llapack -lblas
#endif
SHread2grid: SHread2grid.f90
	$(FC) -o SHread2grid SHread2grid.f90  $(LIBPATH) $(INCLUDE) $(LIBS) -m64 -fpp -free -O3
chaosreadspline: readspline_zyg.f   # 进行6次B样条插值，参数设置在.f文件中进行设置。输出随时间变化的球鞋系数矩阵
	$(FC) -o chaosreadspline readspline_zyg.f -static-intel 
SH_grid: SH_grid.f90
	$(FC) -o SH_grid SH_grid.f90 $(LIBPATH) $(INCLUDE) $(LIBS) -m64 -fpp -free -O3
clean:
	rm -f *.o SHread2grid  chaosreadspline SH_grid

