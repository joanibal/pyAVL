      PROGRAM Parallel_Hello_World
        USE OMP_LIB
        
        !$OMP PARALLEL
        
            PRINT *, "Hello from process: ", OMP_GET_THREAD_NUM()
        
        !$OMP END PARALLEL
        
        END
        ! gfortran hello_omp_world.f -o hello_world -fopenmp