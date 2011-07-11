"""
treatment planning model 
"""
# import needed modules
import petsc4py, numpy, sys
PetscOptions =  sys.argv
PetscOptions.append("-ksp_monitor")
#PetscOptions.append("-help")
petsc4py.init(PetscOptions)

from petsc4py import PETSc
from mpi4py import MPI

# break processors into separate communicators
petscRank = PETSc.COMM_WORLD.getRank()
petscSize = PETSc.COMM_WORLD.Get_size()
sys.stdout.write("petsc rank %d petsc nproc %d\n" % (petscRank, petscSize))

# break up processors into communicators
NumProcsPerSubComm = 4
color = petscRank/NumProcsPerSubComm 
NumSubCommunicators = petscSize / NumProcsPerSubComm + 1
subcomm = PETSc.Comm(MPI.COMM_WORLD.Split(color))
subRank = subcomm.Get_rank()
subSize = subcomm.Get_size()
sys.stdout.write("number of sub communictors %d \n" % ( NumSubCommunicators ))
sys.stdout.write("subcomm rank %d subcomm nproc %d\n" % (subRank, subSize))

# set shell context
import femLibrary
fem = femLibrary.PyFEMInterface(4.0) 
fem.SetuplibMesh( subcomm ) # initialize libMesh data structures

# inifile = None ==> setup from command line
inifile=None
fem.SetupIni( inifile ) 
# set any additional parameters needed
fem.SetIniValue( "optical/mu_a_healthy","4.0" ) 
fem.SetIniValue( "optical/mu_a_tumor","400.0" ) 

# initialize FEM Mesh
# must setup Ini File first
fem.SetupUnStructuredGrid( "phantomMesh.e",0 ) 
#fem.SetupStructuredGrid( (10,10,4) ,[0.0,1.0],[0.0,1.0],[0.0,1.0]) 

# add the data structures for the Background System Solve
# set deltat, number of time steps, power profile, and add system
deltat = 5.0
ntime  = 7 
fem.AddPennesSystem(2,deltat,ntime,[[1,28,46,78,119],[0.0,4.0,0.0,9.0,0.0]]) 

# initialize libMesh data structures
fem.InitializeEquationSystems( ) 

# print info
fem.printSelf() 

# loop over time steps and solve
for timeID in range(1,ntime):
   fem.UpdateTransientSystemTimeStep("StateSystem",timeID ) 
   fem.SystemSolve( "StateSystem" ) 
   #fem.StoreTransientSystemTimeStep("StateSystem",timeID ) 
   fem.WriteTimeStep("fem_data.e" , timeID, timeID )  











