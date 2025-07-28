import numpy as np
from stl import mesh
import os

# CONFIGURATION
AIRFOIL_FILE = "airfoil.txt"
STL_OUTPUT = "constant/triSurface/airfoil.stl"
DOMAIN_LENGTH = 5.0
DOMAIN_HEIGHT = 4.0
CHORD_LENGTH = 1.0
EXTRUDE_DEPTH = 0.01
MESH_CELLS_X = 400
MESH_CELLS_Y = 320

# Prepare directories
def prepare_directories():
    os.makedirs("system", exist_ok=True)
    os.makedirs("constant/triSurface", exist_ok=True)
    os.makedirs("0", exist_ok=True)

# Read and scale airfoil coordinates
def read_airfoil(filename):
    coords = np.loadtxt(filename)
    coords = coords / 1000.0
    max_x = np.max(coords[:, 0])
    coords /= max_x
    # Close the curve if not already closed
    if not np.array_equal(coords[0], coords[-1]):
        coords = np.vstack([coords, coords[0]])
    return coords

# Create STL from scaled coordinates
def create_stl(coords, filename):
    vertices = []
    faces = []
    
    n = len(coords)
    
    # Create vertices for front and back faces
    for x, y in coords:
        vertices.append([x, y, 0.0])           # Front face
        vertices.append([x, y, EXTRUDE_DEPTH]) # Back face
    
    # Create triangular faces for the airfoil surface
    for i in range(n - 1):  # n-1 because last point = first point
        v0_front = i * 2
        v0_back = i * 2 + 1
        v1_front = (i + 1) * 2
        v1_back = (i + 1) * 2 + 1
        
        # Triangle 1: front -> next_front -> next_back
        faces.append([v0_front, v1_front, v1_back])
        # Triangle 2: front -> next_back -> back
        faces.append([v0_front, v1_back, v0_back])
    
    # Convert to STL mesh
    airfoil_mesh = mesh.Mesh(np.zeros(len(faces), dtype=mesh.Mesh.dtype))
    for i, face in enumerate(faces):
        for j in range(3):
            airfoil_mesh.vectors[i][j] = vertices[face[j]]
    
    airfoil_mesh.save(filename)
    print(f"STL saved to {filename}")

# Write controlDict
def write_controlDict():
    with open("system/controlDict", "w") as f:
        f.write("""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version: 2506                                   |
|   \\  /    A nd           |                                                 |
|    \\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      controlDict;
}
application     simpleFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         1000;
deltaT          1;
writeControl    timeStep;
writeInterval   100;
purgeWrite      0;
writeFormat     ascii;
writePrecision 6;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable true;
""")
    print("controlDict written.")

# Write blockMeshDict
def write_blockMeshDict():
    with open("system/blockMeshDict", "w") as f:
        f.write(f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}}
convertToMeters 1.0;
vertices
(
    ({-DOMAIN_LENGTH} {-DOMAIN_HEIGHT} 0)
    ({ DOMAIN_LENGTH} {-DOMAIN_HEIGHT} 0)
    ({ DOMAIN_LENGTH} { DOMAIN_HEIGHT} 0)
    ({-DOMAIN_LENGTH} { DOMAIN_HEIGHT} 0)
    ({-DOMAIN_LENGTH} {-DOMAIN_HEIGHT} {EXTRUDE_DEPTH})
    ({ DOMAIN_LENGTH} {-DOMAIN_HEIGHT} {EXTRUDE_DEPTH})
    ({ DOMAIN_LENGTH} { DOMAIN_HEIGHT} {EXTRUDE_DEPTH})
    ({-DOMAIN_LENGTH} { DOMAIN_HEIGHT} {EXTRUDE_DEPTH})
);
blocks
(
    hex (0 1 2 3 4 5 6 7) ({MESH_CELLS_X} {MESH_CELLS_Y} 1) simpleGrading (1 1 1)
);
boundary
(
    inlet
    {{ type patch; faces ((0 4 7 3)); }}
    outlet
    {{ type patch; faces ((1 2 6 5)); }}
    top
    {{ type patch; faces ((3 7 6 2)); }}
    bottom
    {{ type patch; faces ((0 1 5 4)); }}
    frontAndBack
    {{ type empty; faces ((0 3 2 1)(4 5 6 7)); }}
);
mergePatchPairs ();""")
    print("blockMeshDict written.")

# Write surfaceFeatureExtractDict
def write_surfaceFeatureExtractDict():
    with open("system/surfaceFeatureExtractDict", "w") as f:
        f.write("""FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      surfaceFeatureExtractDict;
}
airfoil.stl
{
    extractionMethod extractFromSurface;
    extractFromSurfaceCoeffs
    {
        includedAngle 150;
    }
    writeObj yes;
}
""")
    print("surfaceFeatureExtractDict written.")

# Write snappyHexMeshDict
def write_snappyHexMeshDict():
    with open("system/snappyHexMeshDict", "w") as f:
        f.write("""FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      snappyHexMeshDict;
}
castellatedMesh true;
snap            true;
addLayers       false;

geometry
{
    airfoil.stl
    {
        type triSurfaceMesh;
        name airfoil;
    }
}

castellatedMeshControls
{
    maxLocalCells       200000;
    maxGlobalCells      400000;
    minRefinementCells  10;
    nCellsBetweenLevels 2;

    features
    (
        {
            file "airfoil.eMesh";
            level 2;
        }
    );

    refinementSurfaces
    {
        airfoil
        {
            level (2 2);
            patchInfo
            {
                type wall;
            }
        }
    }

    refinementRegions
    {
    }

    resolveFeatureAngle 30;
    locationInMesh      (0 0.5 0.005);
    allowFreeStandingZoneFaces true;
}

snapControls
{
    nSmoothPatch 3;
    tolerance    2.0;
    nSolveIter   30;
    nRelaxIter   5;
    nFeatureSnapIter 10;
    explicitFeatureSnap    false;
    implicitFeatureSnap    true;
}

addLayersControls
{
    relativeSizes true;
    layers
    {
    }
}

meshQualityControls
{
    maxNonOrtho          65;
    maxBoundarySkewness  20;
    maxInternalSkewness  4;
    maxConcave           80;
    minVol               1e-13;
    minTetQuality        1e-15;
    minArea              -1;
    minTwist             0.02;
    minDeterminant       0.001;
    minFaceWeight        0.05;
    minVolRatio          0.01;
    minTriangleTwist     -1;
    nSmoothScale         4;
    errorReduction       0.75;
}

debug          0;
mergeTolerance 1e-6;
""")
    print("snappyHexMeshDict written.")

# Write fvSchemes
def write_fvSchemes():
    with open("system/fvSchemes", "w") as f:
        f.write("""FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      fvSchemes;
}

ddtSchemes
{
    default         steadyState;
}

gradSchemes
{
    default         Gauss linear;
}

divSchemes
{
    default         none;
    div(phi,U)      bounded Gauss linearUpwind grad(U);
    div(phi,k)      bounded Gauss upwind;
    div(phi,epsilon) bounded Gauss upwind;
    div(phi,omega) bounded Gauss upwind;
    div(phi,v2)    bounded Gauss upwind;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}

laplacianSchemes
{
    default         Gauss linear corrected;
}

interpolationSchemes
{
    default         linear;
}

snGradSchemes
{
    default         corrected;
}
""")
    print("fvSchemes written.")

# Write fvSolution
def write_fvSolution():
    with open("system/fvSolution", "w") as f:
        f.write("""FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      fvSolution;
}

solvers
{
    p
    {
        solver          GAMG;
        tolerance       1e-06;
        relTol          0.1;
        smoother        GaussSeidel;
    }

    U
    {
        solver          smoothSolver;
        smoother        GaussSeidel;
        tolerance       1e-05;
        relTol          0.1;
    }

    k
    {
        solver          smoothSolver;
        smoother        GaussSeidel;
        tolerance       1e-05;
        relTol          0.1;
    }

    epsilon
    {
        solver          smoothSolver;
        smoother        GaussSeidel;
        tolerance       1e-05;
        relTol          0.1;
    }

    omega
    {
        solver          smoothSolver;
        smoother        GaussSeidel;
        tolerance       1e-05;
        relTol          0.1;
    }
}

SIMPLE
{
    nNonOrthogonalCorrectors 0;
    consistent yes;
}

relaxationFactors
{
    fields
    {
        p       0.3;
    }
    equations
    {
        U       0.7;
        k       0.7;
        epsilon 0.7;
        omega   0.7;
    }
}
""")
    print("fvSolution written.")

if __name__ == '__main__':
    prepare_directories()
    coords = read_airfoil(AIRFOIL_FILE)
    create_stl(coords, STL_OUTPUT)
    write_controlDict()
    write_blockMeshDict()
    write_surfaceFeatureExtractDict()
    write_snappyHexMeshDict()
    write_fvSchemes()
    write_fvSolution()
    print("System files made")
    print("\nTo generate the mesh with airfoil boundary, run these commands in sequence:")
    print("1. blockMesh")
    print("2. surfaceFeatureExtract") 
    print("3. snappyHexMesh -overwrite")