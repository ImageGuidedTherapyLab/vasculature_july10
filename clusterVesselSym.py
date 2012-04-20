# clarox -nographics -nojournal -noecho -batch journalfilename.py
# from cubit cmd script editor
#   f = file("file.py"); cmdtext = f.read() ; exec(cmdtext)

def GenerateRFMesh(VesselDiameter ,VesselDistance):
  """
  Vessel Radius and Distance input in mm
  """
  cubit.cmd('set developer on')
  cubit.cmd('   reset')
  
  #from	Yusheng Feng Yusheng.Feng@utsa.edu
  #to	Carlos Acosta <carlos.acosta.berlin@gmail.com>
  #subject	Blood vessel sizes
  #
  #The following are blood vessel sizes (diameters) that you need to incorporate
  #them for mesh generation for each model.
  #
  #Human:
  #
  #Major Arteries:                 3.0 mm
  #First order arterial branch:    1.0 mm
  #Second order arterial branch:   0.6 mm
  #Arterioles:                     0.02 mm
  #Capillaries:                    0.008 mm
  #
  #Dogs:
  #
  #Major Arteries:                 4.0 mm
  #First order arterial branch:    1.3 mm
  #Second order arterial branch:   0.45 mm
  #Third order arterial branch:    0.15 mm
  #Arterioles:                     0.05 mm
  #Capillaries:                    0.008 mm
  #
  #Thus, you have 10 sizes for diameter, 
  #and coupled with 4 distances: 2.5 mm, 5 mm, 10 mm, and 20 mm.
  # A total of 40 cases.
  #
  #YSF
  
  #create tissue block
  cubit.cmd('create brick x 74 y 74 z 84')
  cubit.cmd('volume 1 move 0 0 12')
  cubit.cmd('volume 1 name "healthy" ')
  idhealthy = cubit.get_id_from_name('healthy')
  print "id" ,idhealthy
  
  #create probe
  ApplicatorRadius   = 1.36 # mm
  ApplicatorDistance = 2.5 # mm
  cubit.cmd('create cylinder radius %12.5e height 120' % ApplicatorRadius)
  cubit.cmd('create cylinder radius %12.5e height 120' % ApplicatorRadius)
  cubit.cmd('volume 2 move  %12.5e  0 0' % ( ApplicatorDistance) )
  cubit.cmd('volume 3 move  %12.5e  0 0' % (-ApplicatorDistance) )
  cubit.cmd('create cylinder radius %12.5e height 120' % ApplicatorRadius)
  # 5 * cos(pi/6)
  ApplicatorHeight = 4.3301270
  cubit.cmd('volume 4 move  0 %12.5e 0' % ApplicatorHeight )
  
  # start distance just outside applicator
  VesselRadius = VesselDiameter / 2.0
  cubit.cmd('create cylinder radius %12.5e height 120' %  VesselRadius )
  DistanceStart = ApplicatorDistance + ApplicatorRadius + VesselRadius #mm
  cubit.cmd('volume 5 move %12.5e 0  0' % (DistanceStart + VesselDistance) )
  
  # cut up domain to mesh
  cubit.cmd('webcut volume 1  tool volume 2 ')
  cubit.cmd('webcut volume 1  tool volume 3 ')
  cubit.cmd('webcut volume 1  tool volume 4 ')
  cubit.cmd('webcut volume 1  tool volume 5 ')
  
  cubit.cmd('delete volume 2 3 4 5')
  # create second vessel
  #cubit.cmd('create cylinder radius %12.5e height 120' %  VesselRadius )
  #cubit.cmd('volume 10 move 0 %12.5e 0' % (ApplicatorHeight+ApplicatorRadius+VesselRadius + ApplicatorDistance ) )
  #cubit.cmd('volume 10 rotate 90 about y ' )
  #cubit.cmd('webcut volume 1  tool volume 10 ')
  #
  # create box to group
  cubit.cmd('create brick x %12.5e y %12.5e z 25' % (2.*(ApplicatorHeight+ApplicatorRadius),2.*(ApplicatorHeight+ApplicatorRadius)) )
  #
  ## cut with box
  cubit.cmd('webcut volume  1 6 7 8 9 with plane surface 44')
  cubit.cmd('webcut volume  1 6 7 8 9 with plane surface 43')
  cubit.cmd('delete volume 10 ')
  # 
  # cubit.cmd('webcut volume 14 Plane Surface 20')
  # cubit.cmd('webcut volume 16 Plane Surface 21')
  # 
  # merge
  cubit.cmd('imprint volume   1 6 7 8 9 11 12 13 14 15 16 17 18 19 20')
  cubit.cmd('merge volume     1 6 7 8 9 11 12 13 14 15 16 17 18 19 20')
  
  # set size
  meshResolutions = [(1.5 , 3.8 ), (0.5 , 3.0 ), (0.25, 2.5 )]
  resolutionID = 0
  (fineSize,coarseSize) = meshResolutions[resolutionID]
  
  cubit.cmd('volume  17 18 19 20 size %f' %   fineSize)
  cubit.cmd('volume  12 13 14 15 6 7 8 9 1 16 11 size %f' % coarseSize)
  ## # 
  ## # mesh 
  cubit.cmd('mesh volume  17 18 19 20')
  cubit.cmd('mesh volume  12 13 14 15')
  cubit.cmd('mesh volume   6  7  8  9')
  cubit.cmd('mesh volume   1 16 11 ')
  ## cubit.cmd('mesh volume    1 ')
  ## cubit.cmd('surface in volume 10 scheme trimesh  ')
  ## cubit.cmd('mesh surface in volume 10  ')
  ## #
  ## 
  ## cubit.cmd('hexset hex in volume  6 7 8  separate tri in surface in volume 10')
  
  # export in pieces
  cubit.cmd('reset genesis')
  cubit.cmd('block 1 volume  12 13 14 1 16 11 6 7 8')
  cubit.cmd('block 1 name "liver"  ')
  cubit.cmd('block 2 volume  9 15 20 ')
  cubit.cmd('block 2 name "vessel"  ')
  cubit.cmd('block 3 volume  17 18 19')
  cubit.cmd('block 3 name "applicator"  ')
  # add BC
  cubit.cmd('skin volume all make sideset 2')
  #cubit.cmd('sideset 2 surface 84 130 118 119 104 105 94')
  cubit.cmd('sideset 2 name "neumann" ')
  cubit.cmd('nodeset 1 volume 17 18 19')
  cubit.cmd('nodeset 1 name "dirichletApplicator"')
  cubit.cmd('nodeset 2 volume  9 15 20')
  cubit.cmd('nodeset 2 name "dirichletVessel"')
  #
  # scale from [mm] to [m] and write'
  cubit.cmd('volume all scale 0.001')
  cubit.cmd('export mesh "clusterVesselDistance%fDiameter%f.e" overwrite' % (VesselDistance,VesselDiameter) )
# end def GenerateRFMesh
##################################################################
def ParseDakotaFile(param_file):
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
                      'vessel_distance' :paramsdict['vessel_distance'],
                      'vessel_diameter' :paramsdict['vessel_diameter'],
                    }
  
  try:
    active_set_vector = [ int(paramsdict['ASV_%d:response_fn_%d' % (i,i) ]) for i in range(1,num_fns+1)  ] 
  except KeyError:
    active_set_vector = [ int(paramsdict['ASV_%d:obj_fn' % (i) ]) for i in range(1,num_fns+1)  ] 
  
  # set a dictionary for passing to rosenbrock via Python kwargs
  fem_params              = {}
  fem_params['cv']        = continuous_vars
  fem_params['asv']       = active_set_vector
  fem_params['functions'] = num_fns

  return fem_params              
# end def ParseDakotaFile:
##################################################################
#params = ParseDakotaFile("/data/fuentes/utsa/vasculature_july10/vessel/realization.1/pce.in")
#vessel_distance = float(params['cv']['vessel_distance'])
#vessel_diameter = float(params['cv']['vessel_diameter'])
diameterList = [0.15,0.45,0.6,1.0,1.3,3.0,4.0] #mm
distanceList = [2.5,5.0,10.0,20.0] #mm

#diameterList = [4.0] #mm
#distanceList = [20.0] #mm
for vessel_diameter in diameterList:
  for vessel_distance in distanceList:
    GenerateRFMesh(vessel_diameter,vessel_distance) 

