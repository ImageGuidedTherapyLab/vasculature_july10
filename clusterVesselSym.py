# clarox -nographics -nojournal -noecho -batch journalfilename.py
# from cubit cmd script editor
#   f = file("file.py"); cmdtext = f.read() ; exec(cmdtext)
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

# create vessel
VesselRadius   = 1.0 # mm
VesselDistance = 0.0 # mm
VesselDistance = 2.5 # mm

# start distance just outside applicator
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
meshResolutions = [(1.0 , 3.8 ), (0.5 , 3.0 ), (0.25, 2.5 )]
resolutionID = 0
(fineSize,coarseSize) = meshResolutions[resolutionID]

cubit.cmd('volume  17 18 19 20 size %f' %   fineSize)
## cubit.cmd('volume  1 10 13    size %f' % coarseSize)
## # 
## # mesh 
cubit.cmd('mesh volume  17 18 19 20 ')
## cubit.cmd('mesh volume    1 ')
## cubit.cmd('surface in volume 10 scheme trimesh  ')
## cubit.cmd('mesh surface in volume 10  ')
## #
## 
## cubit.cmd('hexset hex in volume  6 7 8  separate tri in surface in volume 10')

#cubit.cmd('#')
#cubit.cmd('# export in pieces')
#cubit.cmd('reset genesis')
#cubit.cmd('block 1 volume  10 9 1 8 ')
#cubit.cmd('block 2 volume  7 14 15 ')
#cubit.cmd('block 3 volume  12 13')
#cubit.cmd('# add BC')
#cubit.cmd('sideset 2 surface 84 130 118 119 104 105 94')
#cubit.cmd('sideset 2 name "neumann" ')
#
#cubit.cmd('export mesh "clustertmp.e" overwrite')
#cubit.cmd('reset')
#cubit.cmd('import mesh geometry "clustertmp.e" feature_angle 0')
#cubit.cmd('##')
#cubit.cmd('# add BC')
#cubit.cmd('sideset 2 surface 7 10 ')
#cubit.cmd('sideset 2 name "neumann" ')
#cubit.cmd('sideset 3 surface 6 3 4')
#cubit.cmd('sideset 3 name "cauchy" ')
#cubit.cmd('#skin volume 3 4 make group 3')
#cubit.cmd('nodeset 1 volume 3 4')
#cubit.cmd('nodeset 1 name "dirichlet" ')
#cubit.cmd('##')
#cubit.cmd('# write in pieces ')
#cubit.cmd('block 1 volume 1')
#cubit.cmd('block 1 name "kidney" ')
#cubit.cmd('block 2 volume 2')
#cubit.cmd('block 2 name "vessel" ')
#cubit.cmd('block 3 volume 3 4 ')
#cubit.cmd('block 3 name "probe" ')
#cubit.cmd('#')
#cubit.cmd('# scale from [mm] to [m] and write')
#cubit.cmd('volume all scale 0.001')
#cubit.cmd('${_meshWriteCmd = "export mesh \'clusterVessel"//tostring(_idRes)//".e\' overwrite" }')
#cubit.cmd('{rescan(_meshWriteCmd)}')
#cubit.cmd('comment "res id is" {_idRes} ')
#cubit.cmd('comment "fine size" {_fineSize} "coarse size" {_coarseSize} ')
#cubit.cmd('comment "radius" {_vesselRadius} "distance" {_vesselDistance}')
#
