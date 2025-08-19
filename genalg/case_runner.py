import numpy as np
from stl import mesh
import os
import subprocess
from dotenv import load_dotenv

# --- Load environment variables from .env file ---
load_dotenv()

def get_config():
    """
    Dynamically read configuration from environment variables.
    This ensures we always get the most up-to-date values.
    """
    return {
        'AIRFOIL_FILE': os.getenv("AIRFOIL_FILE"),
        'STL_OUTPUT': os.getenv("STL_OUTPUT"),
        'DOMAIN_LENGTH': float(os.getenv("DOMAIN_LENGTH")),
        'DOMAIN_HEIGHT': float(os.getenv("DOMAIN_HEIGHT")),
        'CHORD_LENGTH': float(os.getenv("CHORD_LENGTH")),
        'EXTRUDE_DEPTH': float(os.getenv("EXTRUDE_DEPTH")),
        'MESH_CELLS_X': int(os.getenv("MESH_CELLS_X")),
        'MESH_CELLS_Y': int(os.getenv("MESH_CELLS_Y")),
        'AIRFOIL_SCALE': float(os.getenv("AIRFOIL_SCALE")),
        'INLET_VELOCITY': float(os.getenv("INLET_VELOCITY")),
        'ANGLE_OF_ATTACK': float(os.getenv("ANGLE_OF_ATTACK")),
        'KINEMATIC_VISCOSITY': float(os.getenv("KINEMATIC_VISCOSITY")),
        'DENSITY': float(os.getenv("DENSITY")),
        'TURBULENCE_MODEL': os.getenv("TURBULENCE_MODEL"),
    }

# --- LEGACY GLOBAL VARIABLES (for backward compatibility) ---
# These will be updated dynamically when needed
AIRFOIL_FILE = os.getenv("AIRFOIL_FILE")
STL_OUTPUT = os.getenv("STL_OUTPUT")
DOMAIN_LENGTH = float(os.getenv("DOMAIN_LENGTH"))
DOMAIN_HEIGHT = float(os.getenv("DOMAIN_HEIGHT"))
CHORD_LENGTH = float(os.getenv("CHORD_LENGTH"))
EXTRUDE_DEPTH = float(os.getenv("EXTRUDE_DEPTH"))
MESH_CELLS_X = int(os.getenv("MESH_CELLS_X"))
MESH_CELLS_Y = int(os.getenv("MESH_CELLS_Y"))
AIRFOIL_SCALE = float(os.getenv("AIRFOIL_SCALE"))
INLET_VELOCITY = float(os.getenv("INLET_VELOCITY"))
ANGLE_OF_ATTACK = float(os.getenv("ANGLE_OF_ATTACK"))
KINEMATIC_VISCOSITY = float(os.getenv("KINEMATIC_VISCOSITY"))
DENSITY = float(os.getenv("DENSITY"))
TURBULENCE_MODEL = os.getenv("TURBULENCE_MODEL")

def prepare_directories():
    """Creates the necessary OpenFOAM directory structure."""
    os.makedirs("system", exist_ok=True)
    os.makedirs("constant/triSurface", exist_ok=True)
    os.makedirs("0", exist_ok=True)

def read_and_scale_airfoil(filename, scale_factor):
    """
    Reads airfoil coordinates, normalizes, applies a scale factor,
    and ensures the geometry is closed.
    """
    config = get_config()  # Get fresh config values
    
    coords = np.loadtxt(filename)
    
    # Normalize by chord length
    max_x = np.max(coords[:, 0])
    min_x = np.min(coords[:, 0])
    chord = max_x - min_x
    coords[:, 0] = (coords[:, 0] - min_x) / chord
    coords[:, 1] = coords[:, 1] / chord
    
    # Apply the scale factor from the genetic algorithm
    coords *= scale_factor

    print("Scale applied is: ", scale_factor)  # Use the actual parameter, not env variable
    
    # Ensure airfoil is closed
    if not np.allclose(coords[0], coords[-1], atol=1e-6):
        coords = np.vstack([coords, coords[0]])
    return coords

def create_stl(coords, filename):
    """Creates a watertight 3D STL file for the airfoil."""
    config = get_config()
    n_points = len(coords) - 1
    
    # Create vertices for both Z planes
    vertices_front = np.hstack([coords[:-1], np.zeros((n_points, 1))])
    vertices_back = np.hstack([coords[:-1], np.full((n_points, 1), config['EXTRUDE_DEPTH'])])
    vertices = np.vstack([vertices_front, vertices_back])

    faces = []
    # Side faces
    for i in range(n_points):
        next_i = (i + 1) % n_points
        faces.append([i, next_i, i + n_points])
        faces.append([next_i, next_i + n_points, i + n_points])

    # Front cap (fan triangulation)
    for i in range(1, n_points - 1):
        faces.append([0, i + 1, i])
        
    # Back cap (fan triangulation)
    offset = n_points
    for i in range(1, n_points - 1):
        faces.append([offset, offset + i, offset + i + 1])

    # Create the mesh
    airfoil_mesh = mesh.Mesh(np.zeros(len(faces), dtype=mesh.Mesh.dtype))
    for i, f in enumerate(faces):
        for j in range(3):
            airfoil_mesh.vectors[i][j] = vertices[f[j]]
    
    airfoil_mesh.save(filename)
    print(f"STL saved to {filename}")

def write_decomposeParDict():
    """Writes the decomposeParDict for parallel processing."""
    num_cores = int(os.getenv("NUM_CORES"))
    with open("system/decomposeParDict", "w") as f:
        f.write(f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2312                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      decomposeParDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

numberOfSubdomains {num_cores};

method          scotch;

// ************************************************************************* //
""")
    print("decomposeParDict written.")

def write_controlDict():
    config = get_config()
    Aref = config['CHORD_LENGTH'] * config['EXTRUDE_DEPTH'] * config['AIRFOIL_SCALE']

    with open("system/controlDict", "w") as f:
        f.write(f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version: 2506                                   |
|   \\  /    A nd           |                                                 |
|    \\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      controlDict;
}}

application     simpleFoam;

startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         5000;

deltaT          1;
writeControl    timeStep;
writeInterval   500;

purgeWrite      0;
writeFormat     ascii;
writePrecision  6;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable true;

libs            ("libforces.so" "libfieldFunctionObjects.so");

functions
{{
    forces
    {{
        type            forces;
        
        writeControl    timeStep;
        writeInterval   1;
        
        patches         (airfoil);
        
        rho             rhoInf;
        rhoInf          {config['DENSITY']};
        
        CofR            (0 0 0);
        liftDir         (0 1 0);
        dragDir         (1 0 0);
        pitchAxis       (0 0 1);
        
        writeFields     false;
    }}

    forceCoeffs
    {{
        type            forceCoeffs;
        
        writeControl    timeStep;
        writeInterval   1;
        
        patches         (airfoil);
        
        rho             rhoInf;
        rhoInf          {config['DENSITY']};
        
        CofR            (0.25 0 0);
        liftDir         (0 1 0);
        dragDir         (1 0 0);
        pitchAxis       (0 0 1);
        
        magUInf         {config['INLET_VELOCITY']};
        lRef            {config['CHORD_LENGTH']};
        Aref            {Aref};
        
        writeFields     false;
    }}
}}
""")
    print("controlDict written with correct syntax.")
    print(f"   Reference area: {Aref:.6f} m²")
    print(f"   Reference velocity: {config['INLET_VELOCITY']:.2f} m/s")
    print(f"   Reference length: {config['CHORD_LENGTH']:.2f} m")


def write_blockMeshDict():
    config = get_config()
    print("Block mesh is:", MESH_CELLS_X)
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
    ({-config['DOMAIN_LENGTH']} {-config['DOMAIN_HEIGHT']} 0)
    ({ config['DOMAIN_LENGTH']} {-config['DOMAIN_HEIGHT']} 0)
    ({ config['DOMAIN_LENGTH']} { config['DOMAIN_HEIGHT']} 0)
    ({-config['DOMAIN_LENGTH']} { config['DOMAIN_HEIGHT']} 0)
    ({-config['DOMAIN_LENGTH']} {-config['DOMAIN_HEIGHT']} {config['EXTRUDE_DEPTH']})
    ({ config['DOMAIN_LENGTH']} {-config['DOMAIN_HEIGHT']} {config['EXTRUDE_DEPTH']})
    ({ config['DOMAIN_LENGTH']} { config['DOMAIN_HEIGHT']} {config['EXTRUDE_DEPTH']})
    ({-config['DOMAIN_LENGTH']} { config['DOMAIN_HEIGHT']} {config['EXTRUDE_DEPTH']})
);
blocks
(
    hex (0 1 2 3 4 5 6 7) ({config['MESH_CELLS_X']} {config['MESH_CELLS_Y']} 1) simpleGrading (1 1 1)
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
addLayers       true;

geometry
{
    airfoil.stl
    {
        type triSurfaceMesh;
        name airfoil;
    }
    
    nearfield
    {
        type searchableBox;
        min (-2.0 -1.0 0.0);
        max ( 4.0  1.0 0.01);
    }
}

castellatedMeshControls
{
    // Limiting the total cell count for faster meshing
    maxLocalCells    100000;
    maxGlobalCells   200000;
    minRefinementCells 10;
    nCellsBetweenLevels 2;

    features
    (
        {
            file "airfoil.eMesh";
            level 2; // Lower feature refinement level for speed
        }
    );

    refinementSurfaces
    {
        airfoil
        {
            level (2 2); // Lower surface refinement level for a coarser mesh
            patchInfo
            {
                type wall;
            }
        }
    }

    refinementRegions
    {
        nearfield
        {
            mode inside;
            levels ((1E15 1)); // Lower region refinement for a faster mesh
        }
    }
    
    resolveFeatureAngle 30;
    locationInMesh      (0 0.5 0.005);
    allowFreeStandingZoneFaces true;
}

snapControls
{
    nSmoothPatch      3;     // Changed from 1
    tolerance         2.0;
    nSolveIter        50;    // Changed from 5 (Crucial for proper snapping)
    nRelaxIter        5;     // Changed from 1 (Helps the solver converge)
    nFeatureSnapIter  10;    // Changed from 2
}

addLayersControls
{
    relativeSizes true;
    
    layers
    {
        airfoil
        {
            nSurfaceLayers 5; // Fewer layers for a faster mesh
        }
    }
    
    // Simplifies layer addition for faster processing
    expansionRatio 1.3;
    finalLayerThickness 0.3;
    minThickness 0.05;
    
    nGrow 0;
    featureAngle 60;
    slipFeatureAngle 30;
    nRelaxIter 1;
    nSmoothSurfaceNormals 1;
    nSmoothNormals 3;
    nSmoothThickness 10;
    maxFaceThicknessRatio 0.5;
    maxThicknessToMedialRatio 0.3;
    minMedialAxisAngle 90;
    nBufferCellsNoExtrude 0;
    nLayerIter 10;
    nRelaxedIter 10;
}

meshQualityControls
{
    maxNonOrtho         65;  // Changed from 75
    maxBoundarySkewness 20;  // Changed from 25
    maxInternalSkewness 4;   // Changed from 8 (This is a very important change)
    maxConcave 80;
    minVol 1e-15;
    minTetQuality 1e-30;
    minArea -1;
    minTwist 0.02;
    minDeterminant 0.001;
    minFaceWeight 0.02;
    minVolRatio 0.01;
    minTriangleTwist -1;
    nSmoothScale 4;
    errorReduction 0.75;
    
    relaxed
    {
        maxNonOrtho 80;
        maxBoundarySkewness 30;
        maxInternalSkewness 10;
        maxConcave 85;
        minVol 1e-20;
        minTetQuality 1e-50;
        minArea -1;
        minTwist 0.001;
        minDeterminant 1e-6;
        minFaceWeight 1e-6;
        minVolRatio 1e-6;
        minTriangleTwist -1;
    }
}

debug 0;
mergeTolerance 1e-6;
""")
    print("snappyHexMeshDict written.")

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
wallDist
{
    method meshWave;
}
""")
    print("fvSchemes written.")

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
        relTol          0.01;
        smoother        GaussSeidel;
        nPreSweeps      0;
        nPostSweeps     2;
        cacheAgglomeration on;
        agglomerator    faceAreaPair;
        nCellsInCoarsestLevel 50;
        mergeLevels     1;
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
    nNonOrthogonalCorrectors 2;
    consistent yes;
    
    residualControl
    {
        p               1e-2;       // Relaxed from 1e-3
        U               1e-3;       // Relaxed from 1e-4
    }
}

relaxationFactors
{
    fields
    {
        p               0.2;        // More aggressive relaxation from 0.3
    }
    equations
    {
        U               0.5;        // More aggressive relaxation from 0.7
        k               0.5;        // More aggressive relaxation from 0.7
        epsilon         0.5;        // More aggressive relaxation from 0.7
        omega           0.5;        // More aggressive relaxation from 0.7
    }
}

""")
    print("fvSolution written.")

def write_transportProperties():
    """Writes the transportProperties file with kinematic viscosity."""
    with open("constant/transportProperties", "w") as f:
        f.write(f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      transportProperties;
}}

transportModel  Newtonian;
nu              {KINEMATIC_VISCOSITY};
""")
    print("transportProperties written.")

def write_turbulenceProperties():
    """
    Writes turbulenceProperties based on the selected model.
    Includes enhanced coefficients for k-omega-SST for low-Re treatment.
    """
    with open("constant/turbulenceProperties", "w") as f:
        if TURBULENCE_MODEL == "laminar":
            f.write("""FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      turbulenceProperties;
}

simulationType laminar;
""")
        else:
            model_map = {
                "kEpsilon": "kEpsilon",
                "kOmega": "kOmega", 
                "kOmegaSST": "kOmegaSST"
            }
            openfoam_model = model_map.get(TURBULENCE_MODEL, "kOmegaSST")
            
            f.write(f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      turbulenceProperties;
}}

simulationType RAS;

RAS
{{
    RASModel        {openfoam_model};
    turbulence      on;
    printCoeffs     on;
""")
            
            if TURBULENCE_MODEL == "kOmegaSST":
                f.write("""    
    kOmegaSSTCoeffs
    {
        alphaK1          0.85;
        alphaK2          1.0;
        alphaOmega1      0.5;
        alphaOmega2      0.856;
        gamma1           0.5532;
        gamma2           0.4403;
        beta1            0.075;
        beta2            0.0828;
        betaStar         0.09;
        a1               0.31;
        b1               1.0;
        c1               10.0;
        
        // Disable F3 blending for better low-Re performance
        F3               false;
    }
""")
            
            f.write("}")
            
    print(f"turbulenceProperties written. Using {TURBULENCE_MODEL} model.")

def calculate_inlet_conditions():
    """
    Calculates the inlet turbulence conditions (k, epsilon, omega)
    based on the turbulence intensity (I) and length scale (L).
    """
    config = get_config()

    # Get values from the configuration
    U_inf = config['INLET_VELOCITY']
    # Use a turbulence intensity from the config file if available, otherwise use a default
    I = config.get('TURBULENCE_INTENSITY', 0.05)
    
    # Calculate the turbulence length scale (L).
    # L is typically based on the chord length of the airfoil.
    # The 'scale' variable can be a multiplier for the chord length.
    base_chord_length = config.get('CHORD_LENGTH', 1.0) # Use a default of 1.0 if not found
    
    # Common practice is to set L as 10% of the chord.
    # Multiply by an additional 'scale' factor if provided in the config.
    scale_factor = config.get('SCALE', 1.0)
    L = 0.1 * base_chord_length * scale_factor

    # Ensure L is not zero to avoid division by zero errors
    if L == 0:
        L = 1e-6 # Set a small non-zero value as a safeguard
    
    # Calculate turbulent kinetic energy
    # k = 1.5 * (U_inf * I)^2
    k = 1.5 * (U_inf * I) ** 2
    
    # Calculate dissipation rate and specific dissipation rate
    # Epsilon = Cmu^(3/4) * k^(3/2) / L
    # OpenFOAM uses the simplified relation epsilon = 0.09 * k^(1.5) / L
    epsilon = 0.09 * k**(1.5) / L
    
    # Omega = epsilon / (Cmu * k)
    # OpenFOAM uses the simplified relation omega = epsilon / (0.09 * k)
    omega = epsilon / (0.09 * k)
    
    return k, epsilon, omega


def calculate_nut():
    config = get_config()
    k, epsilon, omega = calculate_inlet_conditions()
    
    if config['TURBULENCE_MODEL'] == "kEpsilon":
        Cmu = 0.09
        nut = Cmu * k**2 / epsilon
    elif config['TURBULENCE_MODEL'] in ["kOmega", "kOmegaSST"]:
        nut = k / omega
    else:
        nut = 0.0
    
    return nut

def write_boundary_conditions():
    config = get_config()
    
    # Calculate inlet velocity components using FRESH config values
    aoa_rad = np.radians(config['ANGLE_OF_ATTACK'])
    U_x = config['INLET_VELOCITY'] * np.cos(aoa_rad)
    U_y = config['INLET_VELOCITY'] * np.sin(aoa_rad)
    
    k, epsilon, omega = calculate_inlet_conditions()
    nut = calculate_nut()
    
    print(f"   Writing boundary conditions with AoA = {config['ANGLE_OF_ATTACK']:.2f}°")
    
    # Pressure field
    with open("0/p", "w") as f:
        f.write(f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      p;
}}

dimensions      [0 2 -2 0 0 0 0];
internalField   uniform 0;

boundaryField
{{
    inlet
    {{
        type            zeroGradient;
    }}
    
    outlet
    {{
        type            fixedValue;
        value           uniform 0;
    }}
    
    top
    {{
        type            zeroGradient;
    }}
    
    bottom
    {{
        type            zeroGradient;
    }}
    
    airfoil
    {{
        type            zeroGradient;
    }}
    
    frontAndBack
    {{
        type            empty;
    }}
}}
""")
    
    # Velocity field
    with open("0/U", "w") as f:
        f.write(f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volVectorField;
    object      U;
}}

dimensions      [0 1 -1 0 0 0 0];
internalField   uniform ({U_x} {U_y} 0);

boundaryField
{{
    inlet
    {{
        type            fixedValue;
        value           uniform ({U_x} {U_y} 0);
    }}
    
    outlet
    {{
        type            zeroGradient;
    }}
    
    top
    {{
        type            slip;
    }}
    
    bottom
    {{
        type            slip;
    }}
    
    airfoil
    {{
        type            noSlip;
    }}
    
    frontAndBack
    {{
        type            empty;
    }}
}}
""")
    
    # Turbulent kinetic energy
    if config['TURBULENCE_MODEL'] != "laminar":
        with open("0/k", "w") as f:
            f.write(f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      k;
}}

dimensions      [0 2 -2 0 0 0 0];
internalField   uniform {k:.6f};

boundaryField
{{
    inlet
    {{
        type            fixedValue;
        value           uniform {k:.6f};
    }}
    
    outlet
    {{
        type            zeroGradient;
    }}
    
    top
    {{
        type            slip;
    }}
    
    bottom
    {{
        type            slip;
    }}
    
    airfoil
    {{
        type            fixedValue;
        value           uniform 1e-12;
    }}
    
    frontAndBack
    {{
        type            empty;
    }}
}}
""")
    
        # Turbulence dissipation rate or specific dissipation rate
        if config['TURBULENCE_MODEL'] in ["kEpsilon"]:
            with open("0/epsilon", "w") as f:
                f.write(f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      epsilon;
}}

dimensions      [0 2 -3 0 0 0 0];
internalField   uniform {epsilon:.6f};

boundaryField
{{
    inlet
    {{
        type            fixedValue;
        value           uniform {epsilon:.6f};
    }}
    
    outlet
    {{
        type            zeroGradient;
    }}
    
    top
    {{
        type            slip;
    }}
    
    bottom
    {{
        type            slip;
    }}
    
    airfoil
    {{
        type            epsilonWallFunction;
        value           uniform {epsilon:.6f};
    }}
    
    frontAndBack
    {{
        type            empty;
    }}
}}
""")
        
        elif config['TURBULENCE_MODEL'] in ["kOmega", "kOmegaSST"]:
            with open("0/omega", "w") as f:
                f.write(f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      omega;
}}

dimensions      [0 0 -1 0 0 0 0];
internalField   uniform {omega:.6f};

boundaryField
{{
    inlet
    {{
        type            fixedValue;
        value           uniform {omega:.6f};
    }}
    
    outlet
    {{
        type            zeroGradient;
    }}
    
    top
    {{
        type            slip;
    }}
    
    bottom
    {{
        type            slip;
    }}
    
    airfoil
    {{
        type            omegaWallFunction;
        value           uniform {omega:.6f};
    }}
    
    frontAndBack
    {{
        type            empty;
    }}
}}
""")
        
        with open("0/nut", "w") as f:
            f.write(f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      nut;
}}

dimensions      [0 2 -1 0 0 0 0];
internalField   uniform {nut:.6e};

boundaryField
{{
    inlet
    {{
        type            calculated;
        value           uniform {nut:.6e};
    }}
    
    outlet
    {{
        type            calculated;
        value           uniform {nut:.6e};
    }}
    
    top
    {{
        type            calculated;
        value           uniform {nut:.6e};
    }}
    
    bottom
    {{
        type            calculated;
        value           uniform {nut:.6e};
    }}
    
    airfoil
    {{
        type            fixedValue;
        value           uniform 0;
    }}
    
    frontAndBack
    {{
        type            empty;
    }}
}}
""")
    
    else:
        # For laminar flow, still need nut file but with zero values
        with open("0/nut", "w") as f:
            f.write(f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      nut;
}}

dimensions      [0 2 -1 0 0 0 0];
internalField   uniform 0;

boundaryField
{{
    inlet
    {{
        type            calculated;
        value           uniform 0;
    }}
    
    outlet
    {{
        type            calculated;
        value           uniform 0;
    }}
    
    top
    {{
        type            calculated;
        value           uniform 0;
    }}
    
    bottom
    {{
        type            calculated;
        value           uniform 0;
    }}
    
    airfoil
    {{
        type            calculated;
        value           uniform 0;
    }}
    
    frontAndBack
    {{
        type            empty;
    }}
}}
""")
    
    print(f"   Inlet velocity: {U_x:.2f}, {U_y:.2f} m/s (AoA: {config['ANGLE_OF_ATTACK']}°)")
    if config['TURBULENCE_MODEL'] != "laminar":
        print(f"   Turbulence: k={k:.6f}, ε={epsilon:.6f}, ω={omega:.6f}")
        print(f"   Turbulent viscosity: nut={nut:.6e}")
        print("   Using LOW-RE boundary conditions for fine mesh")
    else:
        print("   Laminar flow: nut=0")

def handle_subprocess(command, log_filename, cwd=None):
    """
    A helper function to run a command in a specific directory (cwd)
    and handle errors.
    """
    try:
        # Ensure log files are also written to the correct directory
        log_filepath = os.path.join(cwd, log_filename) if cwd else log_filename
        with open(log_filepath, "w") as log_file:
            process = subprocess.run(
                command,
                check=True,
                stdout=log_file,
                stderr=subprocess.PIPE,
                cwd=cwd # <-- Explicitly sets the working directory for the command
            )
        return True
    except subprocess.CalledProcessError as e:
        print(f"  !!! ERROR during '{' '.join(command)}' !!!")
        print(f"  Error output:\n{e.stderr.decode()}")
        return False
    except FileNotFoundError:
        print(f"  !!! ERROR: Command '{command[0]}' not found. Is OpenFOAM sourced?")
        return False

def run_meshing_workflow(cwd=None):
    """Executes only the meshing steps in the specified directory."""
    print("  Starting meshing workflow...")
    num_cores = int(os.getenv("NUM_CORES"))

    meshing_commands = {
        "blockMesh": ["blockMesh"],
        "surfaceFeatureExtract": ["surfaceFeatureExtract"],
        "decomposePar": ["decomposePar", "-force"],
        "snappyHexMesh": ["mpirun", "-np", str(num_cores), "--oversubscribe", "snappyHexMesh", "-parallel", "-overwrite"],
        "reconstructParMesh": ["reconstructParMesh", "-constant"]
    }

    for step, command in meshing_commands.items():
        print(f"    - Running {step}...")
        if not handle_subprocess(command, f"{step}.log", cwd=cwd):
            return False

    print("  Meshing workflow completed successfully.")
    return True

def run_solver_workflow(cwd=None):
    """Executes only the solver step in the specified directory."""
    print("  Starting solver workflow...")

    if not handle_subprocess(["simpleFoam"], "simpleFoam.log", cwd=cwd):
        return False

    print("  Solver workflow completed successfully.")
    return True