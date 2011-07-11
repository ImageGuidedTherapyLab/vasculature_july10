"""
treatment planning model for RF ablation
"""
# import needed modules
import petsc4py, numpy, sys
PetscOptions =  sys.argv
PetscOptions.append("-ksp_monitor")
PetscOptions.append("-ksp_rtol")
PetscOptions.append("1.0e-15")
#PetscOptions.append("-help")
petsc4py.init(PetscOptions)

from petsc4py import PETSc
from mpi4py import MPI

# break processors into separate communicators
petscRank = PETSc.COMM_WORLD.getRank()
petscSize = PETSc.COMM_WORLD.Get_size()
sys.stdout.write("petsc rank %d petsc nproc %d\n" % (petscRank, petscSize))

# set shell context
# TODO import vtk should be called after femLibrary ???? 
# FIXME WHY IS THIS????
import femLibrary
# initialize libMesh data structures
libMeshInit = femLibrary.PyLibMeshInit(PetscOptions,PETSc.COMM_WORLD) 
  

# inifile = None ==> setup from command line
inifile=None
fem.SetupIni( inifile ) 
# set any additional parameters needed
fem.SetIniValue( "optical/mu_a_healthy","4.0" ) 
fem.SetIniValue( "optical/mu_a_tumor","400.0" ) 

# store control variables
getpot = femLibrary.PylibMeshGetPot(PetscOptions) 
# from Duck table 2.15
getpot.SetIniValue( "material/specific_heat","3840.0" ) 
# set ambient temperature 
getpot.SetIniValue( "initial_condition/u_init","21.0" ) 
# from Duck
getpot.SetIniValue( "electric_conductivity/s_0_healthy", "6.0") 

# initialize FEM Mesh
femMesh = femLibrary.PylibMeshMesh()
#Setup Affine Transformation for registration
RotationMatrix = [[1.,0.,0.],
                  [0.,1.,0.],
                  [0.,0.,1.]]
Translation =     [0.,0.,0.]
# TODO input full path to FEM mesh here
femMesh.SetupUnStructuredGrid( "clusterVessel.e",0,RotationMatrix, Translation  ) 
MeshOutputFile = "fem_data.e" 

# add the data structures for the Background System Solve
# set deltat, number of time steps, power profile, and add system
nsubstep = 1
acquisitionTime = 5.00
deltat = acquisitionTime / nsubstep
ntime  = 60 
eqnSystems =  femLibrary.PylibMeshEquationSystems(femMesh,getpot)
getpot.SetIniPower(nsubstep,  [ [1,5,41,ntime],[1.0,0.0,1.0,0.0] ])
eqnSystems.AddPennesRFSystem("StateSystem",deltat,ntime) 

# initialize libMesh data structures
eqnSystems.init( ) 

# print info
eqnSystems.PrintSelf() 

# write IC
exodusII_IO = femLibrary.PylibMeshExodusII_IO(femMesh)
exodusII_IO.WriteTimeStep(MeshOutputFile,eqnSystems, 1, 0.0 )  

# loop over time steps and solve
for timeID in range(1,ntime):
   print "time step = " ,timeID
   eqnSystems.UpdateTransientSystemTimeStep("StateSystem",timeID ) 
   eqnSystems.SystemSolve( "StateSystem" ) 
   #eqnSystems.StoreTransientSystemTimeStep("StateSystem",timeID ) 
   exodusII_IO.WriteTimeStep(MeshOutputFile,eqnSystems, timeID+1, timeID*deltat )  
