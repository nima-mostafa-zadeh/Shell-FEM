#=========================================================
## INPUTS
#=========================================================

MODEL_NAME = 'Model-2'
JOB_NAME = 'Stretch-EvansSkalak-Test-1'

PART_NAME = 'Part-Cylinder'
CYLINDER_RADIUS = 1.0
CYLINDER_HEIGHT = 1.0
PART_ELEMENTS_NUM = 10

MATERIAL_NAME = 'Material-Cylinder'
MATERIAL_DENSITY = 10.0
MATERIAL_DAMPING = 100.0
MATERIAL_SHEAR_MODULUS = 500.0
MATERIAL_BULK_MODULUS = 10000.0
MATERIAL_POISSON_RATIO = (3 * MATERIAL_BULK_MODULUS - 2 * MATERIAL_SHEAR_MODULUS) / (6 * MATERIAL_BULK_MODULUS + 2 * MATERIAL_SHEAR_MODULUS)
MATERIAL_YOUNGS_MODULUS = 2 * MATERIAL_SHEAR_MODULUS * (1 + MATERIAL_POISSON_RATIO);

SECTION_NAME = 'Section-Cylinder'
SECTION_THICKNESS = 0.01

STEP_IS_AUTOMATIC = False
STEP_TOTAL_TIME = 100
STEP_DT = 5.0e-4
STEP_AMP_MAX_TIME = 10

STRETCHING_CONCENTRATED_FORCE = 0.01

#=========================================================
## LIBRARIES
#=========================================================

from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *

#=========================================================
## MODEL DEFINITION
#=========================================================
MODEL = mdb.Model(modelType=STANDARD_EXPLICIT, name=MODEL_NAME)

#=========================================================
## PART DEFINITION
#=========================================================
MODEL.ConstrainedSketch(name='__profile__', sheetSize=10.0)
MODEL.sketches['__profile__'].sketchOptions.setValues(
    viewStyle=AXISYM)
MODEL.sketches['__profile__'].ConstructionLine(point1=(0.0, 
    -5.0), point2=(0.0, 5.0))
MODEL.sketches['__profile__'].FixedConstraint(entity=
    MODEL.sketches['__profile__'].geometry[2])
MODEL.sketches['__profile__'].Line(point1=(CYLINDER_RADIUS, 0.0), point2=(
    CYLINDER_RADIUS, CYLINDER_HEIGHT))
MODEL.sketches['__profile__'].VerticalConstraint(addUndoState=
    False, entity=MODEL.sketches['__profile__'].geometry[3])
MODEL.Part(dimensionality=AXISYMMETRIC, name=PART_NAME, 
    type=DEFORMABLE_BODY)
MODEL.parts[PART_NAME].BaseWire(sketch=
    MODEL.sketches['__profile__'])
del MODEL.sketches['__profile__']

#=========================================================
## MESHING
#=========================================================
MODEL.parts[PART_NAME].seedPart(deviationFactor=0.1, 
    minSizeFactor=0.1, size=CYLINDER_HEIGHT / PART_ELEMENTS_NUM)
MODEL.parts[PART_NAME].generateMesh()

#=========================================================
## MATERIAL DEFINITION
#=========================================================
MODEL.Material(name=MATERIAL_NAME)
MODEL.materials[MATERIAL_NAME].Density(table=((MATERIAL_DENSITY, ), ))
MODEL.materials[MATERIAL_NAME].Damping(alpha=MATERIAL_DAMPING)
MODEL.materials[MATERIAL_NAME].Depvar(n=4)
MODEL.materials[MATERIAL_NAME].UserMaterial(mechanicalConstants=
    (MATERIAL_SHEAR_MODULUS, MATERIAL_BULK_MODULUS))

#=========================================================
## SECTION DEFINITION
#=========================================================
MODEL.HomogeneousShellSection(idealization=NO_IDEALIZATION, 
    integrationRule=GAUSS, material=MATERIAL_NAME, name=
    SECTION_NAME, nodalThicknessField='', numIntPts=3, poissonDefinition=
    DEFAULT, preIntegrate=OFF, temperature=GRADIENT, thickness=SECTION_THICKNESS, 
    thicknessField='', thicknessModulus=None, thicknessType=UNIFORM, 
    useDensity=OFF)

MODEL.parts[PART_NAME].Set(edges=
    MODEL.parts[PART_NAME].edges.getSequenceFromMask((
    '[#1 ]', ), ), name='Set-1')
MODEL.parts[PART_NAME].SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=
    MODEL.parts[PART_NAME].sets['Set-1'], sectionName=
    SECTION_NAME, thicknessAssignment=FROM_SECTION)

#=========================================================
## ASSEMBLY DEFINITION
#=========================================================
MODEL.rootAssembly.DatumCsysByThreePoints(coordSysType=
    CYLINDRICAL, origin=(0.0, 0.0, 0.0), point1=(1.0, 0.0, 0.0), point2=(0.0, 
    0.0, -1.0))
MODEL.rootAssembly.Instance(dependent=ON, name=
    PART_NAME + '-1', part=MODEL.parts[PART_NAME])

MODEL.rootAssembly.Set(name='Set-Top', vertices=
    MODEL.rootAssembly.instances[PART_NAME + '-1'].vertices.getSequenceFromMask(
    ('[#2 ]', ), ))
MODEL.rootAssembly.Set(name='Set-Anchor', vertices=
    MODEL.rootAssembly.instances[PART_NAME + '-1'].vertices.getSequenceFromMask(
    ('[#1 ]', ), ))

MODEL.rootAssembly.regenerate()

#=========================================================
## STEP DEFINITION
#=========================================================
if STEP_IS_AUTOMATIC:
    MODEL.ExplicitDynamicsStep(improvedDtMethod=ON, name='Step-1', 
        maxIncrement=None, previous='Initial', scaleFactor=1.0, timeIncrementationMethod=
        AUTOMATIC_GLOBAL, timePeriod=STEP_TOTAL_TIME, nlgeom=ON)
else:
    MODEL.ExplicitDynamicsStep(improvedDtMethod=ON, name='Step-1', 
        nlgeom=ON, previous='Initial', timeIncrementationMethod=
        FIXED_USER_DEFINED_INC, timePeriod=STEP_TOTAL_TIME, userDefinedInc=STEP_DT)
        

#=========================================================
## LOAD DEFINITION
#=========================================================
mdb.models['Model-2'].SmoothStepAmplitude(data=((0.0, 0.0), (STEP_AMP_MAX_TIME, 1.0), (STEP_TOTAL_TIME, 1.0)), name=
    'Amp-1', timeSpan=STEP)
    
MODEL.ConcentratedForce(cf2=STRETCHING_CONCENTRATED_FORCE, createStepName='Step-1', 
    distributionType=UNIFORM, field='', localCsys=None, name='Load-1', region=
    MODEL.rootAssembly.sets['Set-Top'], amplitude='Amp-1')

#=========================================================
## BOUNDARY CONDITION DEFINITION
#=========================================================
MODEL.YsymmBC(createStepName='Step-1', localCsys=None, name=
    'BC-1', region=MODEL.rootAssembly.sets['Set-Anchor'])

#=========================================================
## JOB DEFINITION
#=========================================================
mdb.Job(activateLoadBalancing=False, atTime=None, contactPrint=OFF, 
    description='', echoPrint=OFF, explicitPrecision=DOUBLE_PLUS_PACK, 
    historyPrint=OFF, memory=90, memoryUnits=PERCENTAGE, model=MODEL_NAME, 
    modelPrint=OFF, multiprocessingMode=DEFAULT, name=JOB_NAME, 
    nodalOutputPrecision=FULL, numCpus=1, numDomains=1, queue=None, resultsFormat=ODB, 
    scratch='', type=ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
    
    
#mdb.jobs['Stretch-test1'].submit(consistencyChecking=OFF)