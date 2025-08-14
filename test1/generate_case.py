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
MESH_CELLS_Y = 280

# FLOW CONDITIONS
INLET_VELOCITY = 14.6285  # m/s
ANGLE_OF_ATTACK = 0  # degrees
KINEMATIC_VISCOSITY = 1.5e-5  # m²/s (air at 15°C)
DENSITY = 1.225  # kg/m³ (air at 15°C)
TURBULENCE_MODEL = "kOmegaSST"  # kEpsilon, kOmega, kOmegaSST, laminar

# Prepare directories
def prepare_directories():
    os.makedirs("system", exist_ok=True)
    os.makedirs("constant/triSurface", exist_ok=True)
    os.makedirs("0", exist_ok=True)

# 1. Read and scale airfoil coordinates
def read_airfoil(filename):
    """
    Reads airfoil coordinates, scales, normalizes by chord length,
    and ensures the geometry is closed.
    """
    coords = np.loadtxt(filename)
    
    # Scale from mm to m if needed (assuming max > 10 means mm)
    if np.max(coords) > 10:
        coords = coords / 1000.0
    
    # Normalize by chord length
    max_x = np.max(coords[:, 0])
    min_x = np.min(coords[:, 0])
    chord = max_x - min_x
    coords[:, 0] = (coords[:, 0] - min_x) / chord
    coords[:, 1] = coords[:, 1] / chord
    
    # Ensure airfoil is closed by adding the first point at the end
    if not np.allclose(coords[0], coords[-1], atol=1e-6):
        coords = np.vstack([coords, coords[0]])
    
    return coords

def create_stl(coords, filename):
    """
    Create a thin 3D airfoil STL with caps to form a watertight volume.
    """
    n_points = len(coords) - 1 # Remove duplicate closing point for triangulation
    
    # Create vertices for both Z planes
    vertices = []
    for i in range(n_points):
        # Bottom surface (z=0)
        vertices.append([coords[i][0], coords[i][1], 0.0])
    
    for i in range(n_points):
        # Top surface (z=EXTRUDE_DEPTH)
        vertices.append([coords[i][0], coords[i][1], EXTRUDE_DEPTH])
    
    # Create triangular faces
    faces = []
    
    # 1. Side faces connecting bottom and top surfaces (your original code)
    for i in range(n_points):
        next_i = (i + 1) % n_points
        
        # Triangle 1: bottom[i], bottom[next_i], top[i] (outward normal)
        faces.append([i, next_i, i + n_points])
        
        # Triangle 2: bottom[next_i], top[next_i], top[i] (outward normal)
        faces.append([next_i, next_i + n_points, i + n_points])

    # 2. Front cap faces (at z=0)
    # We use a fan triangulation with the first point as the center
    # The normal should point in the -z direction
    for i in range(1, n_points - 1):
        faces.append([0, i + 1, i])
        
    # 3. Back cap faces (at z=EXTRUDE_DEPTH)
    # Normal should point in the +z direction
    offset = n_points # Index offset for the top vertices
    for i in range(1, n_points - 1):
        faces.append([offset, offset + i, offset + i + 1])

    # Create the mesh object
    airfoil_mesh = mesh.Mesh(np.zeros(len(faces), dtype=mesh.Mesh.dtype))
    
    for i, face in enumerate(faces):
        for j in range(3):
            airfoil_mesh.vectors[i][j] = vertices[face[j]]
    
    airfoil_mesh.save(filename)
    print(f"Watertight 3D STL saved to {filename} with {len(faces)} triangular faces")


# Write controlDict
def write_controlDict():
    Aref = CHORD_LENGTH * EXTRUDE_DEPTH  # Assuming unit span

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
writeInterval   100;

purgeWrite      0;
writeFormat     ascii;
writePrecision  6;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable true;

functions
{{
    forces
    {{
        type            forces;
        libs            ("libforces.so");
        
        writeControl    timeStep;
        writeInterval   1;
        
        patches         (airfoil);
        
        rho             rhoInf;
        rhoInf          {DENSITY};
        
        CofR            (0 0 0);
        liftDir         (0 1 0);
        dragDir         (1 0 0);
        pitchAxis       (0 0 1);
    }}

    forceCoeffs
    {{
        type            forceCoeffs;
        libs            ("libforces.so");
        
        writeControl    timeStep;
        writeInterval   1;
        
        patches         (airfoil);
        
        rho             rhoInf;
        rhoInf          {DENSITY};
        
        CofR            (0.25 0 0);
        liftDir         (0 1 0);
        dragDir         (1 0 0);
        pitchAxis       (0 0 1);
        
        magUInf         {INLET_VELOCITY};
        lRef            {CHORD_LENGTH};
        Aref            {Aref};
    }}

    pressureCoeff
    {{
        type            pressure;
        libs            ("libfieldFunctionObjects.so");
        
        writeControl    writeTime;
        
        mode            staticCoeff;
        
        rho             rhoInf;
        rhoInf          {DENSITY};
        
        pInf            0;
        UInf            ({INLET_VELOCITY} 0 0);
    }}

    wallPressure
    {{
        type            surfaces;
        libs            ("libsampling.so");
        
        writeControl    writeTime;
        
        surfaceFormat   raw;
        
        surfaces
        (
            airfoil
            {{
                type        patch;
                patches     (airfoil);
                triangulate false;
            }}
        );
        
        fields (p wallShearStress);
    }}
}}
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

# Write snappyHexMeshDict - ORIGINAL WORKING VERSION
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
addLayers       true;  // Changed to enable boundary layers

geometry
{
    airfoil.stl
    {
        type triSurfaceMesh;
        name airfoil;
    }
    
    // Nearfield refinement box around airfoil
    nearfield
    {
        type searchableBox;
        min (-2.0 -1.0 0.0);     // Adjust coordinates based on your airfoil position
        max ( 4.0  1.0 0.01);    // Extend further downstream for wake capture
    }
    
}

castellatedMeshControls
{
    maxLocalCells       1000000;  // Increased for boundary layers
    maxGlobalCells      2000000; // Increased for boundary layers
    minRefinementCells  10;
    nCellsBetweenLevels 2;

    features
    (
        {
            file "airfoil.eMesh";
            level 3;
        }
    );

    refinementSurfaces
    {
        airfoil
        {
            level (4 4);  // Increased refinement for better boundary layer attachment
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
            levels ((1E15 1));  // Level 2 refinement in nearfield
        }
        
    }

    resolveFeatureAngle 15;
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
        airfoil
        {
            nSurfaceLayers 15;    // INCREASED from 8 - much higher BL resolution
        }
    }
    
    // High-resolution boundary layer parameters
    expansionRatio    1.1;        // REDUCED from 1.2 - tighter spacing
    finalLayerThickness 0.2;      // REDUCED from 0.5 - thinner layers
    minThickness      0.02;       // REDUCED from 0.1 - thinner minimum
    
    // Enhanced layer controls for accuracy
    nGrow             2;          // INCREASED from 0 - buffer cells
    featureAngle      45;         // REDUCED from 60 - better feature handling
    slipFeatureAngle  20;         // REDUCED from 30 - tighter slip control
    nRelaxIter        5;          // INCREASED from 3 - more relaxation
    nSmoothSurfaceNormals 3;      // INCREASED from 1 - smoother normals
    nSmoothNormals    10;         // INCREASED from 3 - much smoother
    nSmoothThickness  20;         // INCREASED from 10 - smoother thickness
    maxFaceThicknessRatio 0.3;    // REDUCED from 0.5 - more uniform layers
    maxThicknessToMedialRatio 0.2; // REDUCED from 0.3 - better quality
    minMedialAxisAngle 70;        // REDUCED from 90 - stricter criteria
    nBufferCellsNoExtrude 1;      // INCREASED from 0 - buffer zones
    nLayerIter        100;        // INCREASED from 50 - more iterations
    nRelaxedIter      50;         // INCREASED from 20 - more relaxed iterations
    
    // Additional advanced controls
    additionalReporting true;     // Better layer reporting
    mergePatchFacesAngle 30;      // Merge patch faces for smoother layers
}

meshQualityControls
{
    // Relaxed quality controls for boundary layer meshing
    maxNonOrtho          70;    // Slightly relaxed for layers
    maxBoundarySkewness  20;
    maxInternalSkewness  4;
    maxConcave           80;
    minVol               1e-13;
    minTetQuality        1e-30; // Much more relaxed for layers
    minArea              -1;
    minTwist             0.02;
    minDeterminant       0.001;
    minFaceWeight        0.02;  // Relaxed for layers
    minVolRatio          0.01;
    minTriangleTwist     -1;
    nSmoothScale         4;
    errorReduction       0.75;
    
    // Additional quality controls for layers
    relaxed
    {
        maxNonOrtho         75;
        maxBoundarySkewness 25;
        maxInternalSkewness 8;
        maxConcave          80;
        minVol              1e-20;
        minTetQuality       1e-50;
        minArea             -1;
        minTwist            0.001;
        minDeterminant      1e-6;
        minFaceWeight       1e-6;
        minVolRatio         1e-6;
        minTriangleTwist    -1;
    }
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
wallDist
{
    method meshWave;
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
        relTol          0.01;  // Less strict
        smoother        GaussSeidel;
        nPreSweeps      0;
        nPostSweeps     2;
        cacheAgglomeration on;
        agglomerator    faceAreaPair;
        nCellsInCoarsestLevel 50;  // Increased
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
    nNonOrthogonalCorrectors 1;  // Add non-orthogonal correctors
    consistent yes;
    
    residualControl
    {
        p               1e-3;  // Less strict
        U               1e-4;
    }
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

# Write transport properties
def write_transportProperties():
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

# Write turbulence properties
def write_turbulenceProperties():
    with open("constant/turbulenceProperties", "w") as f:
        turbulence_models = {
            "laminar": "laminar",
            "kEpsilon": "kEpsilon", 
            "kOmega": "kOmega",
            "kOmegaSST": "kOmegaSST"
        }
        
        model = turbulence_models.get(TURBULENCE_MODEL, "kEpsilon")
        
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
    RASModel        {model};
    turbulence      on;
    printCoeffs     on;
}}
""")
    print(f"turbulenceProperties written ({TURBULENCE_MODEL})")

# Calculate inlet turbulence conditions
def calculate_inlet_conditions():
    U_inf = INLET_VELOCITY
    I = 0.05  # 5% turbulence intensity
    L = 0.1   # turbulence length scale (10% of chord)
    
    k = 1.5 * (U_inf * I) ** 2
    epsilon = 0.09 * k**(1.5) / L
    omega = epsilon / (0.09 * k)
    
    return k, epsilon, omega

# NEW: Calculate turbulent viscosity
def calculate_nut():
    k, epsilon, omega = calculate_inlet_conditions()
    
    if TURBULENCE_MODEL == "kEpsilon":
        # For k-epsilon: nut = Cmu * k^2 / epsilon
        Cmu = 0.09
        nut = Cmu * k**2 / epsilon
    elif TURBULENCE_MODEL in ["kOmega", "kOmegaSST"]:
        # For k-omega: nut = k / omega
        nut = k / omega
    else:
        # For laminar flow
        nut = 0.0
    
    return nut

# Write boundary condition files
def write_boundary_conditions():
    # Calculate inlet velocity components
    aoa_rad = np.radians(ANGLE_OF_ATTACK)
    U_x = INLET_VELOCITY * np.cos(aoa_rad)
    U_y = INLET_VELOCITY * np.sin(aoa_rad)
    
    k, epsilon, omega = calculate_inlet_conditions()
    nut = calculate_nut()
    
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
    if TURBULENCE_MODEL != "laminar":
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
        if TURBULENCE_MODEL in ["kEpsilon"]:
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
        
        elif TURBULENCE_MODEL in ["kOmega", "kOmegaSST"]:
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
    
    print(f"   Inlet velocity: {U_x:.2f}, {U_y:.2f} m/s (AoA: {ANGLE_OF_ATTACK}°)")
    if TURBULENCE_MODEL != "laminar":
        print(f"   Turbulence: k={k:.6f}, ε={epsilon:.6f}, ω={omega:.6f}")
        print(f"   Turbulent viscosity: nut={nut:.6e}")
        print("   Using LOW-RE boundary conditions for fine mesh")
    else:
        print("   Laminar flow: nut=0")


def write_turbulence_properties():
    """Write turbulenceProperties with enhanced k-omega-SST coefficients for low-Re treatment"""
    
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
            # Map turbulence models to OpenFOAM names
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
            
            # Add enhanced coefficients for k-omega-SST low-Re treatment
            if TURBULENCE_MODEL == "kOmegaSST":
                f.write("""    
    kOmegaSSTCoeffs
    {
        alphaK1         0.85;
        alphaK2         1.0;
        alphaOmega1     0.5;
        alphaOmega2     0.856;
        gamma1          0.5532;
        gamma2          0.4403;
        beta1           0.075;
        beta2           0.0828;
        betaStar        0.09;
        a1              0.31;
        b1              1.0;
        c1              10.0;
        
        // Enhanced low-Re treatment
        F3              false;      // Disable F3 blending for better low-Re performance
    }
""")
            
            f.write("}")
            
    print(f"   Turbulence model: {TURBULENCE_MODEL}")
    if TURBULENCE_MODEL == "kOmegaSST":
        print("   Enhanced k-omega-SST coefficients for low-Re treatment")


def calculate_inlet_conditions():
    """
    Calculate inlet turbulence conditions for high Reynolds number external flow
    Using very low turbulence intensity appropriate for clean wind tunnel conditions
    """
    # Very low turbulence intensity for clean external flow
    turbulence_intensity = 0.001  # 0.1% - very clean conditions
    
    # Turbulent kinetic energy
    k = 1.5 * (INLET_VELOCITY * turbulence_intensity) ** 2
    
    # Turbulent length scale (1% of chord length for external flow)
    L = 0.01  # 1% of 1m chord
    
    # Specific dissipation rate (omega)
    C_mu = 0.09
    omega = k**0.5 / (C_mu**0.25 * L)
    
    # Dissipation rate (epsilon) - if needed
    epsilon = C_mu * k**2 / (k / omega)
    
    return k, epsilon, omega


def calculate_nut():
    """
    Calculate turbulent viscosity for external flow
    Using much lower values appropriate for low turbulence intensity
    """
    # For very low turbulence intensity external flow
    # Use turbulent viscosity ratio μt/μ ≈ 0.1-1.0
    kinematic_viscosity = 1.5e-5  # Air at 20°C
    turbulent_viscosity_ratio = 0.143  # Conservative for external flow
    
    nut = turbulent_viscosity_ratio * kinematic_viscosity
    
    return nut


# Main script
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
    write_transportProperties()
    write_turbulenceProperties()
    write_boundary_conditions()
    
    print("All system files generated. Run (Don't forget to source /usr/lib/openfoam/openfoam2506/etc/bashrc)")
    print("1. blockMesh")
    print("2. surfaceFeatureExtract") 
    print("3. decomposePar, then mpirun -np 4 snappyHexMesh -parallel -overwrite (for 4 cores mpi) OR snappyHexMesh -overwrite")
    print("4. reconstructParMesh -constant -overwrite if you used mpi")
    print("Note: If you get an error in Paraview that says the mesh count don't match, you needa remove 0/ file then run the script again once you do snappyHexMesh")
    print("5. simpleFoam")
    print("")
    print(f"   INLET_VELOCITY = {INLET_VELOCITY} m/s")
    print(f"   ANGLE_OF_ATTACK = {ANGLE_OF_ATTACK}°") 
    print(f"   TURBULENCE_MODEL = {TURBULENCE_MODEL}")