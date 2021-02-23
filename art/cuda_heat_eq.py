import pycuda.autoinit
import pycuda.driver as drv
import numpy
from pylab import *
from pycuda.compiler import SourceModule


class CUDADiffusion(object):
    def __init__(self,array):
        self.ROWS,self.COLS = shape(array)

        mod = SourceModule("""
        const int COLS = """+str(self.COLS)+""";
        const int ROWS = """+str(self.ROWS)+""";

        __global__ void Laplace(float *T_old, float *T_new, float *bc)
        {
            // compute the "i" and "j" location of the node point
            // handled by this thread

            int i = blockIdx.x * blockDim.x + threadIdx.x ;
            int j = blockIdx.y * blockDim.y + threadIdx.y ;

            // get the natural index values of node (i,j) and its neighboring nodes
                                          //                         N
            int P = i + j*COLS;           // node (i,j)              |
            int N = i + (j+1)*COLS;       // node (i,j+1)            |
            int S = i + (j-1)*COLS;       // node (i,j-1)     W ---- P ---- E
            int E = (i+1) + j*COLS;       // node (i+1,j)            |
            int W = (i-1) + j*COLS;       // node (i-1,j)            |
                                          //                         S

            // only update "interior" node points
            if(i>0 && i<COLS-1 && j>0 && j<ROWS-1) {
                T_new[P] = 0.25*( T_old[E] + T_old[W] + T_old[N] + T_old[S] );
            }
            T_new[P] = (bc[P]==-1)*T_new[P] + (bc[P]!=-1)*bc[P];
        }
        """,options=['-std=c++11'])


        ## initialize memory on GPU
        self.u  = drv.mem_alloc(array.nbytes)
        self.next_u = drv.mem_alloc(array.nbytes)
        self.boundary_conditions = drv.mem_alloc(array.nbytes)

        ##
        self.cuda_laplace = mod.get_function("Laplace")
        self.cuda_laplace.prepare("P P P")

    def equilibriate(self,init_con,boundary_conditions,tol,n_its=100):
        """ takes the initial condition and boundary conditions. The latter are specified 
        using values other than -1, and approximates equilibrium."""

        u_h      = np.zeros_like(init_con,dtype=np.float32)
        next_u_h = np.zeros_like(init_con)

        drv.memcpy_htod(self.u, init_con)
        drv.memcpy_htod(self.boundary_conditions, boundary_conditions)
        
        block_w = 32
        block_h = 32
        grid_w  = int(ceil(float(self.COLS)/block_w))
        grid_h  = int(ceil(float(self.ROWS)/block_h))

        
        n_calls = 5
        it = 0
        while True:
            it += 1
            self.cuda_laplace.prepared_call((grid_w,grid_h),
                                            (block_w,block_h,1), # grid
                                            self.u,self.next_u,
                                            self.boundary_conditions) # block
            if (it % n_calls) == 0:
                """ check magnitude of change """
                drv.memcpy_dtoh(u_h,self.u)
                drv.memcpy_dtoh(next_u_h,self.next_u)
                max_change = abs(u_h-next_u_h).max()
                if max_change < tol:
                    # print('it: %d \t max_change: %f'%(it,max_change))
                    break
                n_calls *= 2
            drv.memcpy_dtod(self.u,self.next_u,init_con.nbytes)

        drv.memcpy_dtoh(init_con,self.next_u)
        
        
        # u_host = numpy.zeros_like(init_con)
        # drv.memcpy_dtoh(u_host, self.u)
        # du_host = numpy.zeros_like(init_con)
        # drv.memcpy_dtoh(du_host, self.next_u)
        # bc_host = numpy.zeros_like(init_con)
        # drv.memcpy_dtoh(bc_host, self.boundary_conditions)
        # figure()
        # subplot2grid((1,3),(0,0))
        # imshow(u_host)
        # subplot2grid((1,3),(0,1))
        # imshow(du_host)
        # subplot2grid((1,3),(0,2))
        # imshow(bc_host)
        # show()



# COLS,ROWS = 256,128
# init_con = np.zeros((COLS,ROWS)).astype(numpy.float32)

# test = cuda_heat_eq(init_con)

# ## setup initial conditions
# r = 20
# init_con[COLS/2-r:COLS/2+r,ROWS/2-r:ROWS/2+r] = 1.0

# ## setup boundary conditions
# boundary_conditions = -1*np.ones_like(init_con)
# boundary_conditions[1:5,:] = 1.0
# boundary_conditions[50:55,:] = 1.0
# boundary_conditions[:,60:65] = 0.0

# test.equilibriate(init_con,boundary_conditions)

# COLS,ROWS = 256,128
# init_con = np.zeros((COLS,ROWS),dtype=np.float32)

# test = CUDAHeatEquilibriator(init_con)

# ## setup initial conditions
# r = 20
# init_con[COLS/2-r:COLS/2+r,ROWS/2-r:ROWS/2+r] = 1.0

# ## setup boundary conditions
# boundary_conditions = -1*np.ones_like(init_con,dtype=np.float32)
# boundary_conditions[1:5,:] = 1.0
# boundary_conditions[50:55,:] = 1.0
# boundary_conditions[:,60:65] = 0.0

# test.equilibriate(init_con,boundary_conditions,1E-4)
