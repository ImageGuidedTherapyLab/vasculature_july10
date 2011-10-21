# Read DAKOTA parameters file (aprepro or standard format) and call a
# Python module for fem analysis.

# DAKOTA will execute this script as
#   linearWFS.py params.in results.out
# so sys.argv[1] will be the parameters file and
#    sys.argv[2] will be the results file to return to DAKOTA

# necessary python modules
import sys
import re
import os

# FIXME global vars are prob bad idea...
ntime  = 10 
nsubstep = 1
acquisitionTime = 6.00
deltat = acquisitionTime / nsubstep
SolnOutputTemplate = "soln.%04d.out"

# set imaging dimensions of final projection
imageDimensions = (10,10,20)

def WriteCubitJouFile(diameter,distance):
  print """
  reset
  #create tissue
  create brick x 74 y 74 z 84
  volume 1 move 0 0 12
  
  #create probe
  create cylinder radius 1.36 height 25
  volume 2 move -2.5  0 0
  create cylinder radius 1.36 height 25
  volume 3 move  2.5  0 0
  create cylinder radius 1.36 height 25
  # 5 * cos(pi/6)
  volume 4 move  0  4.3301270 0
  create cylinder radius 8 height 25
  
  # create vessel
  # create cylinder radius 4.0 height 120
  create cylinder radius %s height 120
  volume 6 move %s 0  0
  
  # cut up domain to mesh
  webcut volume 1 with tool volume 6
  webcut volume 1 with Plane Surface 18
  webcut volume 8 with Plane Surface 17
  webcut volume 9 with tool volume 5
  webcut volume 10 with tool volume 2
  webcut volume 10 with tool volume 3
  webcut volume 10 with tool volume 4
  webcut volume 7  with Plane Surface 18
  webcut volume 14 with Plane Surface 17
  
  # clean up
  delete volume 2 3 4 5 6
  
  # only mesh half
  webcut volume 1 8 9 10 13 with plane xplane 
  delete volume 11 16 17 18 19 20 
  
  # merge
  imprint volume  1 7 8 9 10 12 13 14 15
  merge volume    1 7 8 9 10 12 13 14 15
  
  # set size
  volume 1 7 8 9 14 15 size {_coarseSize}
  volume  10 12 13     size {_fineSize}
  
  # mesh tets
  mesh volume 10 12 13
  mesh volume 9 1 7 8 14 15 
  
  #
  # export in pieces
  reset genesis
  block 1 volume  10 9 1 8 
  block 2 volume  7 14 15 
  block 3 volume  12 13
  # add BC
  sideset 2 surface 84 130 118 119 104 105 94
  sideset 2 name "neumann" 
  
  export mesh "clustertmp.e" overwrite
  reset
  import mesh geometry "clustertmp.e" feature_angle 0
  ##
  # add BC
  sideset 2 surface 7 10 
  sideset 2 name "neumann" 
  sideset 3 surface 6 3 4
  sideset 3 name "cauchy" 
  #skin volume 3 4 make group 3
  nodeset 1 volume 3 4
  nodeset 1 name "dirichlet" 
  ##
  # write in pieces 
  block 1 volume 1
  block 1 name "kidney" 
  block 2 volume 2
  block 2 name "vessel" 
  block 3 volume 3 4 
  block 3 name "probe" 
  #
  # scale from [mm] to [m] and write
  volume all scale 0.001
  ${_meshWriteCmd = "export mesh 'clusterVessel"//tostring(_idRes)//".e' overwrite" }
  {rescan(_meshWriteCmd)}
  comment "res id is" {_idRes} 
  comment "fine size" {_fineSize} "coarse size" {_coarseSize} 
  comment "radius" {_vesselRadius} "distance" {_vesselDistance}
  """  %  (diameter,distance) 

# Write Template Image
def WriteVTKTemplateImage( TemplateFilename ):
  import vtk
  import vtk.util.numpy_support as vtkNumPy 
  import numpy
  # set Image Template Dimensions
  femBounds  = (0.0,0.04,-0.04, 0.04, -0.03,0.06)
  origin = (femBounds[0], femBounds[2], femBounds[4])
  spacing = ( (femBounds[1]-femBounds[0])/ imageDimensions[0] ,
              (femBounds[3]-femBounds[2])/ imageDimensions[1] ,
              (femBounds[5]-femBounds[4])/ imageDimensions[2]  
            )
  print femBounds, origin, spacing
  # imports raw data and stores it.
  dataImporter = vtk.vtkImageImport()
  # array is converted to a string of chars and imported.
  # numpy array stored as ROW MAJOR
  # MUST write out in COLUMN MAJOR format to be the same as VTK
  data_string = numpy.zeros(imageDimensions,dtype=numpy.float,order='F').tostring(order='F')
  dataImporter.CopyImportVoidPointer(data_string, len(data_string))
  # The type of the newly imported data is set to unsigned char (uint8)
  dataImporter.SetDataScalarTypeToDouble()
  # Because the data that is imported only contains an intensity value (it isnt RGB-coded or someting similar), the importer
  # must be told this is the case.
  dataImporter.SetNumberOfScalarComponents(1)
  # The following two functions describe how the data is stored and the dimensions of the array it is stored in. For this
  # simple case, all axes are of length 75 and begins with the first element. For other data, this is probably not the case.
  # I have to admit however, that I honestly dont know the difference between SetDataExtent() and SetWholeExtent() although
  # VTK complains if not both are used.
  dataImporter.SetDataExtent( 0, imageDimensions[0]-1, 0, imageDimensions[1]-1, 0, imageDimensions[2]-1)
  dataImporter.SetWholeExtent(0, imageDimensions[0]-1, 0, imageDimensions[1]-1, 0, imageDimensions[2]-1)
  dataImporter.SetDataSpacing( spacing )
  dataImporter.SetDataOrigin(  origin  )
  dataImporter.SetScalarArrayName(  "scalars" )
  dataImporter.Update()
  vtkTemplateWriter = vtk.vtkDataSetWriter()
  vtkTemplateWriter.SetFileName( TemplateFilename )
  vtkTemplateWriter.SetInput( dataImporter.GetOutput() )
  vtkTemplateWriter.Update()
  return 

####################################################################
def RunSubProcessQueue(CodeExeCmds,ErrorFile):
  import subprocess
  import time
  # max number of local jobs
  numlocalJob = 8 
  process = []; 
  # create an error log file
  errorFile=open(ErrorFile ,"w")
  while ( len(CodeExeCmds) > 0 or len(process) > 0 ):
      # only run numlocalJob at a time
      if (len(process) > numlocalJob):
        raise RuntimeError("\n\n running too many jobs at a time??")
      elif (len(process) == numlocalJob):
        print len(CodeExeCmds), " jobs remaining..."
        time.sleep(30) # pause wait for jobs to finish
      elif( len(CodeExeCmds) > 0 ):
        cmd = CodeExeCmds.pop(0)
        print "running " , cmd
        process.append( [subprocess.Popen(cmd,shell=True),cmd] )
      if( len(process) > 0 ):
        runningJob = process.pop(0)
        if ( runningJob[0].poll() == None ):
          # job not done put it back in the list
          # not that we pop from the front of the list and pushback at the 
          # end to cycle through
          print " pid ",runningJob[0].pid, " still running"
          process.append( runningJob )
          time.sleep(2) # pause 
        elif ( runningJob[0].poll() == 0 ):
          pass # job is done 
        else:
          print "job exiting with ", runningJob[0].poll() 
          errorFile.write("error in  %s   \n" % runningJob[1] )
          #raise RuntimeError("\n\n unknown exit code ")
  errorFile.close; errorFile.flush() 
####################################################################
def rfModeling(**kwargs):
  """
  treatment planning model 
  """
  ## # import petsc and numpy
  ## import petsc4py, numpy
  ## # init petsc
  ## PetscOptions =  sys.argv
  ## PetscOptions.append("-ksp_monitor")
  ## PetscOptions.append("-ksp_rtol")
  ## PetscOptions.append("1.0e-15")
  ## #PetscOptions.append("-help")
  ## petsc4py.init(PetscOptions)
  ## 
  ## # break processors into separate communicators
  ## from petsc4py import PETSc
  ## petscRank = PETSc.COMM_WORLD.getRank()
  ## petscSize = PETSc.COMM_WORLD.Get_size()
  ## sys.stdout.write("petsc rank %d petsc nproc %d\n" % (petscRank, petscSize))

  ## # set shell context
  ## # TODO import vtk should be called after femLibrary ???? 
  ## # FIXME WHY IS THIS????
  ## import femLibrary
  ## # initialize libMesh data structures
  ## libMeshInit = femLibrary.PyLibMeshInit(PetscOptions,PETSc.COMM_WORLD) 
  ## 
  ## # store control variables
  ## getpot = femLibrary.PylibMeshGetPot(PetscOptions) 
  ## # set dirichlet bc on vessel
  ## getpot.SetIniValue( "bc/u_dirichletid","1" ) 
  ## # from Duck table 2.15
  ## getpot.SetIniValue( "material/specific_heat","3840.0" ) 
  ## # set ambient temperature 
  ## getpot.SetIniValue( "initial_condition/u_init","21.0" ) 
  ## # from Duck
  ## getpot.SetIniValue( "electric_conductivity/s_0_healthy", "6.0") 

  ## # from Duck table 2.15
  ## getpot.SetIniValue( "material/specific_heat","3840.0" ) 
  ## # set ambient temperature 
  ## u_init = 34.3
  ## getpot.SetIniValue( "initial_condition/u_init","%f" % u_init ) 
  ## getpot.SetIniValue( "initial_condition/probe_init","%f" % u_init ) 
  ## # thermal conductivity from Duck/CRC Handbook
  ## getpot.SetIniValue( "thermal_conductivity/k_0_healthy",
  ##                             kwargs['cv']['k_0_healthy'] ) 
  ## getpot.SetIniValue( "thermal_conductivity/k_0_tumor",
  ##                             kwargs['cv']['k_0_tumor'] ) 
  ## # perfusion from Duck/CRC Handbook
  ## getpot.SetIniValue( "perfusion/w_0_healthy",
  ##                  kwargs['cv']['w_0_healthy'] ) 
  ## getpot.SetIniValue( "perfusion/w_0_tumor",
  ##                  kwargs['cv']['w_0_tumor'] ) 
  ## #
  ## #  given the original orientation as two points along the centerline z = x2 -x1
  ## #     the transformed orienteation would be \hat{z} = A x2 + b - A x1 - b = A z
  ## #  ie transformation w/o translation which is exactly w/ vtk has implemented w/ TransformVector
  ## #  TransformVector = TransformPoint - the transation
  ## #Setup Affine Transformation for registration
  ## RotationMatrix = [[1.,0.,0.],
  ##                   [0.,1.,0.],
  ##                   [0.,0.,1.]]
  ## Translation =     [0.,0.,0.] 
  ## # original coordinate system laser input
  ## laserTip    = Translation 
  ## 
  ## # set laser orientation values
  ## getpot.SetIniValue( "probe/domain","2") 
  ## getpot.SetIniValue( "probe/x_0","%f" % laserTip[0]) 
  ## getpot.SetIniValue( "probe/y_0","%f" % laserTip[1]) 
  ## getpot.SetIniValue( "probe/z_0","%f" % laserTip[2]) 
  ## 
  ## # initialize FEM Mesh
  ## femMesh = femLibrary.PylibMeshMesh()
  #femMesh.SetupUnStructuredGrid(kwargs['mesh_file'],0,RotationMatrix, Translation  ) 
  # TODO input full path to FEM mesh here
  print "reading vessel mesh with diameter %s distance %s" % (kwargs['cv']['vessel_diameter'],kwargs['cv']['vessel_distance'])
  WriteCubitJouFile(kwargs['cv']['vessel_diameter'],kwargs['cv']['vessel_distance'])
  FEMMeshFileName="/data/fuentes/utsa/vasculature_july10/clusterVessel.e"
  ## femMesh.ReadFile(FEMMeshFileName)
  ## # http://en.wikipedia.org/wiki/Tmpfs
  ## # tmpfs /dev/shm file system should write to RAM = fast!
  ## MeshOutputFile = "/dev/shm/fem_data.%04d.e" % kwargs['fileID'] 
  ## #fem.SetupStructuredGrid( (10,10,4) ,[0.0,1.0],[0.0,1.0],[0.0,1.0]) 
  ## 
  ## # add the data structures for the Background System Solve
  ## # set deltat, number of time steps, power profile, and add system
  ## eqnSystems =  femLibrary.PylibMeshEquationSystems(femMesh,getpot)
  ## getpot.SetIniPower(nsubstep,[ [1,2,4,7,ntime],[0.0,4.0,0.0,9.0,0.0] ])
  ## rfSystem = eqnSystems.AddPennesRFSystem("StateSystem",deltat) 
  ## rfSystem.AddStorageVectors(ntime)
  ## 
  ## # initialize libMesh data structures
  ## eqnSystems.init( ) 
  ## 
  ## # quick error check
  ## #errCheckSoln = fem.GetSolutionVector( "StateSystem" )[...]
  ## #if (  errCheckSoln.size + 1 != kwargs['functions'] ):
  ## #  print "ERROR!! number of response functions incorrect!!"
  ## #  raise RuntimeError("soln vec + 1 = %d .NE. num_response = %d"%(errCheckSoln.size+1,kwargs['functions']) )

  ## # print info
  ## eqnSystems.PrintSelf() 
  ## 
  ## # write IC
  ## exodusII_IO = femLibrary.PylibMeshExodusII_IO(femMesh)
  ## exodusII_IO.WriteTimeStep(MeshOutputFile,eqnSystems, 1, 0.0 )  
  ## 
  ## # loop over time steps and solve
  ## ObjectiveFunction = 0.0
  ## for timeID in range(1,ntime*nsubstep):
  ## #for timeID in range(1,3):
  ##    print "time step = " ,timeID
  ##    eqnSystems.UpdatePetscFEMSystemTimeStep("StateSystem",timeID ) 
  ##    rfSystem.SystemSolve( ) 
  ##    # write soln to disk for processing
  ##    if ( timeID%nsubstep == 0 ):
  ##      exodusII_IO.WriteTimeStep(MeshOutputFile ,eqnSystems, timeID+1, timeID*deltat )  
  ##      # project to imaging data
  ##      if ( petscRank == 0 ):
  ##        import vtk
  ##        import vtk.util.numpy_support as vtkNumPy 
  ##        # read exodus file
  ##        vtkExodusIIReader = vtk.vtkExodusIIReader()
  ##        vtkExodusIIReader.SetFileName(MeshOutputFile)
  ##        vtkExodusIIReader.SetPointResultArrayStatus("u0",1)
  ##        vtkExodusIIReader.SetTimeStep(timeID-1) 
  ##        vtkExodusIIReader.Update()
  ##        # read template image
  ##        vtkTemplateReader = vtk.vtkDataSetReader() 
  ##        vtkTemplateReader.SetFileName( "../imageTemplate.vtk" ) 
  ##        vtkTemplateReader.Update() 
  ##        # project to imaging
  ##        vtkResample = vtk.vtkCompositeDataProbeFilter()
  ##        vtkResample.SetInput(  vtkTemplateReader.GetOutput() )
  ##        vtkResample.SetSource( vtkExodusIIReader.GetOutput() ) 
  ##        vtkResample.Update()
  ##        # write numpy to disk
  ##        imageData = vtkResample.GetOutput().GetPointData().GetArray('u0')
  ##        soln      = vtkNumPy.vtk_to_numpy(imageData) 
  ##        numpy.savetxt(SolnOutputTemplate % timeID,soln)

  ## # clean up
  ## os.remove(MeshOutputFile)

  # return values
  retval = dict([])
  retval['fns'] = [0.0]
  ##retval['fns'] = [ObjectiveFunction]
  ##retval['rank'] = petscRank 
  return(retval)
# end def rfModeling(**kwargs):
####################################################################
def WriteDakotaInputFile(filename,analysis_driver,
                         parameters_file,results_file,dakota_work_directory,
                         allow_existing_results,num_response):
  dakotafile=open( filename,"w")
  # TODO better to pass None and capture ?
  if(allow_existing_results == None):
    allow_existing_results  = ""
  dakotafile.write(
"""
# see reference manual
#   http://dakota.sandia.gov/licensing/stable/html-ref/index.html

# run single method and output to a file
#  http://dakota.sandia.gov/licensing/stable/html-ref/StratCommands.html
strategy,
	single_method tabular_graphics_data

method,
	polynomial_chaos
	  #quadrature_order   = 4 4 4 4 4 4 1 1
	  quadrature_order   = 1 1 1 1 1 1 4 4
	  samples = 10000		
	  seed = 12347 rng rnum2	
          # vector response input 
          # http://dakota.sandia.gov/licensing/stable/html-ref/IntroCommands.html#IntroCmdsInpSpec
	  response_levels =		
	  %d*10.:10.:100.
	  variance_based_decomp #univariate_effects
	  #output verbose
	  output silent

variables,
	uniform_uncertain = 8			
      	  lower_bounds      =  4.  .59   3.   4.  .59   3. .0001  .0004
	  upper_bounds      =  8.  .61   12.  8.  .61  12. .003   .004
	  descriptors       = 's_0_healthy' 'k_0_healthy' 'w_0_healthy' 's_0_tumor' 'k_0_tumor' 'w_0_tumor' 'vessel_distance' 'vessel_diameter'

interface,
	system 
          asynchronous			#0,#p0
          evaluation_concurrency 10	#3,#8,#18,#19
	  analysis_driver = '%s'
	  #analysis_driver = './ibrun_par_driver'
	  # this will guarantee that evaluations are replaced with
          # evaluations modulo the evaluation concurrency
 	  local_evaluation_static_scheduling
	  parameters_file = '%s'
	  results_file = '%s'
          work_directory named = "%s"
	  file_save file_tag
	  directory_save directory_tag
          %s

responses,
	num_response_functions = %d
	no_gradients
	no_hessians
""" % (num_response,analysis_driver,parameters_file,results_file,dakota_work_directory,allow_existing_results,num_response) )
  dakotafile.close; dakotafile.flush() 
# end def WriteDakotaInputFile
##################################################################
def ParseAndRun(param_file):
  # ----------------------------
  # Parse DAKOTA parameters file
  # ----------------------------
  
  # setup regular expressions for parameter/label matching
  e = '-?(?:\\d+\\.?\\d*|\\.\\d+)[eEdD](?:\\+|-)?\\d+' # exponential notation
  f = '-?\\d+\\.\\d*|-?\\.\\d+'                        # floating point
  i = '-?\\d+'                                         # integer
  value = e+'|'+f+'|'+i                                # numeric field
  tag = '\\w+(?::\\w+)*'                               # text tag field
  
  # regular expression for aprepro parameters format
  aprepro_regex = re.compile('^\s*\{\s*(' + tag + ')\s*=\s*(' + value +')\s*\}$')
  # regular expression for standard parameters format
  standard_regex = re.compile('^\s*(' + value +')\s+(' + tag + ')$')
  
  # open DAKOTA parameters file for reading
  paramsfile = open(param_file, 'r')
  #fileID = int(param_file.split(".").pop())
  fileID = int(os.getcwd().split(".").pop())
  
  # extract the parameters from the file and store in a dictionary
  paramsdict = {}
  for line in paramsfile:
      m = aprepro_regex.match(line)
      if m:
          paramsdict[m.group(1)] = m.group(2)
      else:
          m = standard_regex.match(line)
          if m:
              paramsdict[m.group(2)] = m.group(1)
  
  paramsfile.close()
  
  # crude error checking; handle both standard and aprepro cases
  num_vars = 0
  if ('variables' in paramsdict):
      num_vars = int(paramsdict['variables'])
  elif ('DAKOTA_VARS' in paramsdict):
      num_vars = int(paramsdict['DAKOTA_VARS'])
  
  num_fns = 0
  if ('functions' in paramsdict):
      num_fns = int(paramsdict['functions'])
  elif ('DAKOTA_FNS' in paramsdict):
      num_fns = int(paramsdict['DAKOTA_FNS'])
  
  # -------------------------------
  # Convert and send to application
  # -------------------------------
  
  # set up the data structures the rosenbrock analysis code expects
  # for this simple example, put all the variables into a single hardwired array
  continuous_vars = { 
                      'k_0_healthy'     :paramsdict['k_0_healthy'],
                      'k_0_tumor'       :paramsdict['k_0_tumor'  ],
                      's_0_healthy'     :paramsdict['s_0_healthy'],
                      's_0_tumor'       :paramsdict['s_0_tumor'  ],
                      'vessel_distance' :paramsdict['vessel_distance'],
                      'vessel_diameter' :paramsdict['vessel_diameter'],
                    }
  try:
     continuous_vars['w_0_healthy'] = paramsdict['w_0_healthy' ]  
     continuous_vars['w_0_tumor'  ] = paramsdict['w_0_tumor'   ] 
  except KeyError:
     continuous_vars['w_0_healthy'] = "0.0"
     continuous_vars['w_0_tumor'  ] = "0.0"
  
  try:
    active_set_vector = [ int(paramsdict['ASV_%d:response_fn_%d' % (i,i) ]) for i in range(1,num_fns+1)  ] 
  except KeyError:
    active_set_vector = [ int(paramsdict['ASV_%d:obj_fn' % (i) ]) for i in range(1,num_fns+1)  ] 
  
  # set a dictionary for passing to rosenbrock via Python kwargs
  fem_params              = {}
  fem_params['cv']        = continuous_vars
  fem_params['asv']       = active_set_vector
  fem_params['functions'] = num_fns
  fem_params['fileID']    = fileID 
  
  # execute the rosenbrock analysis as a separate Python module
  print "Running RF model..."
  fem_results = rfModeling(**fem_params)
  print "RF complete."
  
  
  ## ----------------------------
  ## Return the results to DAKOTA
  ## ----------------------------
  #
  #if (fem_results['rank'] == 0 ):
  #  # write the results.out file for return to DAKOTA
  #  # this example only has a single function, so make some assumptions;
  #  # not processing DVV
  #  outfile = open('results.out.tmp.%d' % fileID, 'w')
  #  
  #  # write functions
  #  for func_ind in range(0, num_fns):
  #      if (active_set_vector[func_ind] & 1):
  #          functions = fem_results['fns']    
  #          outfile.write(str(functions[func_ind]) + ' f' + str(func_ind) + '\n')
  #  
  #  ## write gradients
  #  #for func_ind in range(0, num_fns):
  #  #    if (active_set_vector[func_ind] & 2):
  #  #        grad = rosen_results['fnGrads'][func_ind]
  #  #        outfile.write('[ ')
  #  #        for deriv in grad: 
  #  #            outfile.write(str(deriv) + ' ')
  #  #        outfile.write(']\n')
  #  #
  #  ## write Hessians
  #  #for func_ind in range(0, num_fns):
  #  #    if (active_set_vector[func_ind] & 4):
  #  #        hessian = rosen_results['fnHessians'][func_ind]
  #  #        outfile.write('[[ ')
  #  #        for hessrow in hessian:
  #  #            for hesscol in hessrow:
  #  #                outfile.write(str(hesscol) + ' ')
  #  #            outfile.write('\n')
  #  #        outfile.write(']]')
  #  #
  #  outfile.close();outfile.flush
  #  #
  #  ## move the temporary results file to the one DAKOTA expects
  #  #import shutil
  #  #shutil.move('results.out.tmp.%d' % fileID, sys.argv[2])
# end def ParseAndRun:
##################################################################

# setup command line parser to control execution
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--dakota_exe", 
                  action="store", dest="dakota_exe", default="dakota",
                  help="full path to dakota EXE", metavar="EXE")
parser.add_option("--work_dir", action="store", dest="work_dir", 
                  default="vessel",
                  help="path to work DIR", metavar="DIR")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
parser.add_option( "--pre_run", "--setup",
                  action="store", dest="pre_run", default=None,
                  help="path to driver setup FILE=[mpich2_setup_driver,ranger_setup_driver]", metavar="FILE")
parser.add_option( "--run_fem","--param_file", 
                  action="store", dest="param_file", default=None,
                  help="run code with parameter FILE", metavar="FILE")
parser.add_option( "--run_queue",
                  action="store", dest="run_queue", default=None,
                  help="run jobs in DIR", metavar="DIR")
parser.add_option( "--execution",
                  action="store", dest="execution", default="sh",
                  help="run jobs with EXE to run exe.qsub", metavar="EXE")
parser.add_option( "--post_run",
                  action="store", dest="post_run", default=None,
                  help="post process stats in DIR ", metavar = "DIR")
(options, args) = parser.parse_args()

# write initial dakota input scripts to setup directories and files
if (options.pre_run != None):
  # create work dir
  if(os.path.isdir(options.work_dir)):
    raise RuntimeError("work dir %s ALREADY EXISTS!!!!" % options.work_dir)
  os.mkdir(options.work_dir)
  # all in/out files should be different to avoid potential race condition
  WriteDakotaInputFile("%s/pce_setup.in" % options.work_dir,
                       "%s/%s" % (os.getcwd(),options.pre_run) , "pce.in", "pce.out" ,
                       "realization", None,1)
  # get # nodes
  num_func = imageDimensions[0] * imageDimensions[1] * imageDimensions[2] 
  POSTEXEC=[]
  for idtime in range(1,ntime):
      timeDir = "%s/time.%04d" % (options.work_dir,idtime)
      os.mkdir( timeDir )
      # all in/out files should be different to avoid potential race condition
      WriteDakotaInputFile("%s/pce_time.%04d.in" %  (timeDir,idtime) ,
                           "echo " , "pce.%04d.in" % idtime ,SolnOutputTemplate % idtime ,
                           "../realization" ,"allow_existing_results",num_func)
      fullpathJob =  "%s/%s/time.%04d" %(os.getcwd(),options.work_dir,idtime)
      POSTEXEC.append("cd %s; %s pce_time.%04d.in > /dev/null" % ( fullpathJob,options.dakota_exe,idtime ) )
  # write script for post processing stats
  postRunStatFile=open( "%s/paramlist" % options.work_dir ,"w")
  for cmdEXE in POSTEXEC:
    postRunStatFile.write(cmdEXE )
    postRunStatFile.write("\n" )
  postRunStatFile.close; postRunStatFile.flush() 
  # set it up...
  os.system("cd %s; %s pce_setup.in" % (options.work_dir,options.dakota_exe) )
  # create image template and write to disk to visualize WRT mesh
  WriteVTKTemplateImage( "%s/imageTemplate.vtk" % options.work_dir )

#run all jobs in the queue
elif (options.run_queue != None):
  CODEEXEC=[]
  for job in os.listdir(options.run_queue):
    fullpathJob =  "%s/%s/%s" %(os.getcwd(),options.run_queue,job)
    # only consider directories
    if(os.path.isdir(fullpathJob)):
      # filter out directories not of interest
      if( job.find("realization") != -1 ) :
        CODEEXEC.append("cd %s; %s exe.qsub" % (fullpathJob,options.execution) )
  RunSubProcessQueue(CODEEXEC,"%s/error_run.log" % options.run_queue)
elif (options.param_file != None):
  ParseAndRun(options.param_file)
elif (options.post_run):
  POSTEXEC=[]
  for idtime in range(1,ntime):
    fullpathJob =  "%s/%s/time.%04d" %(os.getcwd(),options.post_run,idtime)
    POSTEXEC.append("cd %s; %s pce_time.%04d.in > /dev/null" % ( fullpathJob,options.dakota_exe,idtime ) )
  RunSubProcessQueue(POSTEXEC,"%s/error_post.log" % options.post_run)
else:
  parser.print_help()
  print options
