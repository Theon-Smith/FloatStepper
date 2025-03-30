//
// Created by theon on 2025/03/26.
//
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    overInterDyMFoam

Group
    grpMultiphaseSolvers grpMovingMeshSolvers

Description
    Solver for two incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach,
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "incompressibleInterPhaseTransportModel.H"
//#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include "cellCellStencilObject.H"
#include "localMin.H"
#include "oversetAdjustPhi.H"
#include "oversetPatchPhiErr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Experimental extension of floatStepper for overset meshes"
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"

    #include "createDyMControls.H"
    #include "createFields.H"

    // might need to set correctphi flag set in fvSolution
    #include "initCorrectPhi.H"
    //#include "createUfIfPresent.H"
    #include "createUf.H"

    // In OverInterOnly:
    #include "createAlphaFluxes.H"  //todo: probably needs a rework. MULES dependent
    #include "createFvOptions.H"
    #include "createControls.H"
    #include "setCellMask.H"    // crucial overset stuff
    #include "setInterpolatedCells.H"

    // turbulence->validate(); hopefully don't need this anymore

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "readOversetDyMControls.H"
        if (LTS)
        {
            #include "setRDeltaT.H"    //is this necessary? not present in floatstepper
        }
        else
        {
            #include "CourantNo.H"
            #include "alphaCourantNo.H"
            #include "setDeltaT.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        //scalar timeBeforeMeshUpdate = runTime.elapsedCpuTime();    //No mention of this. Looks like its diagnostic

        mesh.update();

        //////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////// This stuff should probably be moved to the mesh update section in ??floatermotion???
        /////////////////////////////////////////////////////////////////////////////////////////////////////
        if (mesh.changing())
        {
            // Do not apply previous time-step mesh compression flux
            // if the mesh topology changed
            if (mesh.topoChanging())
            {
                talphaPhi1Corr0.clear();
            }

            // Update cellMask field for blocking out hole cells
            #include "setCellMask.H"
            #include "setInterpolatedCells.H"
            #include "correctPhiFaceMask.H"
        }
        ///////////////////////////////////////////////////////////////////////////////////////////////////////
        // There is also an MRF.update() after these two lines in floatstepper

        gh = (g & mesh.C()) - ghRef;
        ghf = (g & mesh.Cf()) - ghRef;

        // floatstepper stuff. unsure about placement
        MRF.update();
/*
        if (correctPhi)
        {
            // Calculate absolute flux
            // from the mapped surface velocity
            // phi = mesh.Sf() & Uf(); todo: reintroduce this

            #include "correctPhi.H"

            // Make the flux relative to the mesh motion
            fvc::makeRelative(phi, U);

            mixture.correct();
        }*/

        if (checkMeshCourantNo)
        {
            #include "meshCourantNo.H"
        }

        // not sure this is the best place for this
        U.correctBoundaryConditions();

        while (pimple.loop())
        {
            // No clue what this is. Just going to leave it in todo: understand this
            if (adjustFringe)
            {
                oversetAdjustPhi(phi, U, zoneIdMass);
            }

            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            rhoPhi *= faceMask;

            mixture.correct();

	        #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
 	           #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}
