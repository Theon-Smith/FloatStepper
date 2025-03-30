//
// Created by theon on 2025/03/28.
//

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "dynamicOverFloaterMotionSolversFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "motionSolver.H"
#include "pointMesh.H"
#include "volFields.H"
#include "floaterMotionSolver.H"
#include "fvCFD.H"
#include "isoAdvection.H"
#include "myCorrectPhi.H"
#include "pimpleControl.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "incompressibleInterPhaseTransportModel.H"
#include "processorPolyPatch.H" //Used by calculateGhf0.H

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicOverFloaterMotionSolversFvMesh, 0);
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        dynamicOverFloaterMotionSolversFvMesh,
        IOobject
    );
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        dynamicOverFloaterMotionSolversFvMesh,
        doInit
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicOverFloaterMotionSolversFvMesh::dynamicOverFloaterMotionSolversFvMesh
(
    const IOobject& io,
    const bool doInit
)
:
    dynamicFvMesh(io, doInit),
    motionSolvers_()
{
    if (doInit)
    {
        init(false);    // do not initialise lower levels
    }
}


bool Foam::dynamicOverFloaterMotionSolversFvMesh::init
(
    const bool doInit,
    const bool mandatory
)
{
    if (doInit)
    {
        dynamicFvMesh::init(doInit);
    }

    IOobject ioDict
    (
        "dynamicMeshDict",
        time().constant(),
        *this,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    IOdictionary dict(ioDict);

    label i = 0;

    const auto* dictptr = dict.findDict("solvers");
    if (dictptr)
    {
        const dictionary& solverDict = *dictptr;

        motionSolvers_.setSize(solverDict.size());

        for (const entry& dEntry : solverDict)
        {
            if (dEntry.isDict())
            {
                IOobject io(ioDict);
                io.readOpt(IOobject::NO_READ);
                io.writeOpt(IOobject::AUTO_WRITE);
                io.rename(dEntry.dict().dictName());

                IOdictionary IOsolverDict
                (
                    io,
                    dEntry.dict()
                );

                motionSolvers_.set
                (
                    i++,
                    motionSolver::New(*this, IOsolverDict)
                );
            }
        }
        motionSolvers_.setSize(i);
    }
    else if (mandatory)
    {
        motionSolvers_.setSize(1);
        motionSolvers_.set(i++, motionSolver::New(*this, dict));
    }

    // Assume something changed
    return true;
}


bool Foam::dynamicOverFloaterMotionSolversFvMesh::init(const bool doInit)
{
    // Fall-back to always constructing motionSolver
    return init(doInit, true);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicOverFloaterMotionSolversFvMesh::~dynamicOverFloaterMotionSolversFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::dynamicOverFloaterMotionSolversFvMesh::mapFields
(
    const mapPolyMesh& mpm
)
{
    dynamicFvMesh::mapFields(mpm);

    // Update the motionSolvers for any topo change ...
    for (auto& ms : motionSolvers_)
    {
        ms.updateMesh(mpm);
    }
}


bool Foam::dynamicOverFloaterMotionSolversFvMesh::update()
{
  /*//////////////////////////////////////////////////////
    //// Old solver codes
    /////////////////////////////////////////////////////
    if (motionSolvers_.size())
    {
        // Accumulated displacement
        pointField disp(motionSolvers_[0].newPoints() - fvMesh::points());

        for (label i = 1; i < motionSolvers_.size(); i++)
        {
            disp += motionSolvers_[i].newPoints() - fvMesh::points();
        }

        fvMesh::movePoints(points() + disp);

        volVectorField* Uptr = getObjectPtr<volVectorField>("U");

        if (Uptr)
        {
            Uptr->correctBoundaryConditions();
        }
    }

    return true;

   */
    if (motionSolvers_.size())
    {
        // Accumulated displacement
        pointField disp(points().size(), Zero);

        // First accummulate all non-floaterMotion displacements
        forAll(motionSolvers_, i)
        {
            if (!isA<floaterMotionSolver>(motionSolvers_[i]))
            {
                Info << "motionSolver_[" << i << "] is not a floaterMotionSolver" << endl;
                if (isA<displacementMotionSolver>(motionSolvers_[i]))
                {
                    Info << "is a displacementMotionSolver" << endl;
                    displacementMotionSolver& p0Solver =
                       dynamic_cast<displacementMotionSolver&>
                       (
                           motionSolvers_[i]
                       );
                   p0Solver.solve();
                   disp += p0Solver.pointDisplacement().primitiveField();
                }
            }
        }

        // Then add accumulated displacements to floaterMotion displacements
        forAll(motionSolvers_, i)
        {
            if (isA<floaterMotionSolver>(motionSolvers_[i]))
            {
                Info << "motionSolver_[" << i << "] is a floaterMotionSolver" << endl;

                // Grab handle to mesh motion solver
                floaterMotionSolver& bodySolver =
                    static_cast<floaterMotionSolver&>(motionSolvers_[i]);
                Info << "Body state before any update:" << endl;
                bodySolver.motion().status();

                const scalar deltaT = time().deltaTValue();

                if (time().timeIndex() == time().startTimeIndex() + 1)
                {
                    Info << "First time step - not calculating added mass" << endl;
                    vector a = bodySolver.motion().state().a();
                    vector alpha = bodySolver.motion().state().domegadt();
                    scalarField dvwdt(6, 0);
                    dvwdt[0] = a[0], dvwdt[1] = a[1], dvwdt[2] = a[2];
                    dvwdt[3] = alpha[0], dvwdt[4] = alpha[1], dvwdt[5] = alpha[2];
                    bodySolver.motion().updateFloaterState(dvwdt, deltaT);
                    //Info << "Body state after first time step:" << endl;
                    //bodySolver.motion().status();
                }
                else
                {
                    // todo: ensure necessary overset things are saved
                    pointField oldPoints(polyMesh::points());
//                    pointField oldOldPoints(polyMesh::oldPoints());
                    // oldMeshPhi used in alternative way to reset mesh motion (see below)
                    surfaceScalarField oldMeshPhi(mesh.phi());
                    const floaterMotionState oldBodyState(bodySolver.motion().state());
                    const scalarField noAcceleration(6,0);
                    bodySolver.motion().setAcceleration(Zero);
                    bodySolver.motion().setAngularAcceleration(Zero);
                    bodySolver.motion().updateFloaterState(noAcceleration, deltaT);
                    //Info << "Body state after zero acceleration time step:" << endl;
                    //bodySolver.motion().status();

                    // Update mesh points in accordance with 0-acceleration body motion
                    #include "updateOversetMesh0.H"
                    // Calculate fluid motion corresponding to 0-acceleration body motion
                    #include "updateFluid0.H" //alphaEqn, UEqn and Piso loop
                    // Recording fluid force associated with 0-acceleration motion
                    vector F0(Zero);
                    vector tau0(Zero);
                    bodySolver.motion().calcForceAndTorque(rho0, p0, U0, F0, tau0);
                    Info << "F0 = " << F0 << ", tau0 = " << tau0 << endl;

                    // todo: check that oversetting is properly restored
                    // Resetting mesh going first two steps back, then one forward to
                    // regain correct mesh.phi()
                    // fvMesh::movePoints(oldOldPoints);
                    // fvMesh::movePoints(oldPoints);
                    // Alternative approach with storage of old meshPhi
                    mesh.setPhi() = oldMeshPhi;
                    mesh.movePoints(oldPoints);

                    // restore overset mesh
                    oversetFvMeshBase::update();
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
                    // end of overset restore

                    //Reset body state
                    bodySolver.motion().setState(oldBodyState);
                    //Info << "Body state after reset:" << endl;
                    //bodySolver.motion().status();

                    // Calculating the 6 columns of the added mass tensor
                    const labelList DoFs = bodySolver.DoFs();

                    if (bodySolver.motion().MaddUpdateTime())
                    {
                        moving(false);
                        bodySolver.motion().calcAddedMass(*this, DoFs);
                        moving(true);
                    }
                    // Calculating 6 element acceleration vector
                    scalarField dvwdt = bodySolver.motion().calcAcceleration(F0, tau0, DoFs);
                    // Updating body state to new time
                    bodySolver.motion().updateFloaterState(dvwdt, deltaT);
                    Info << "Body state after real time step:" << endl;
                    bodySolver.motion().status();
                    Info << "Centre of mass      : " << bodySolver.motion().centreOfMass() << endl;

                }

                bodySolver.solve();
                disp += bodySolver.pointDisplacement().primitiveField();
                fvMesh::movePoints(bodySolver.points0() + disp);
                // todo: check that this is a good jumping off point for oversetting
            }
        }


        volVectorField* Uptr = getObjectPtr<volVectorField>("U");

        if (Uptr)
        {
            Uptr->correctBoundaryConditions();
        }
    }

    return true;
}


// ************************************************************************* //
