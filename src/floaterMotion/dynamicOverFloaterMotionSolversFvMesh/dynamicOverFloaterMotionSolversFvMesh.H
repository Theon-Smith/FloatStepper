//
// Created by theon on 2025/03/28.
//

#ifndef DYNAMICOVERFLOATERMOTIONSOLVERSFVMESH_H
#define DYNAMICOVERFLOATERMOTIONSOLVERSFVMESH_H


/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

Class
    Foam::dynamicOversetFvMesh

Description
    dynamicFvMesh with support for overset meshes.

SourceFiles
    dynamicOversetFvMesh.C

\*---------------------------------------------------------------------------*/

#include "dynamicFvMesh.H"
#include "oversetFvMeshBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class motionSolver;

/*---------------------------------------------------------------------------*\
                   Class dynamicOverFloaterMotionSolversFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicOverFloaterMotionSolversFvMesh
:
    public dynamicFvMesh,
    public oversetFvMeshBase
{
    // Private data
    PtrList<motionSolver> motionSolvers_;

    // Private Member Functions

    //- No copy construct
    dynamicOverFloaterMotionSolversFvMesh(const dynamicOverFloaterMotionSolversFvMesh&) = delete;

    //- No copy assignment
    void operator=(const dynamicOverFloaterMotionSolversFvMesh&) = delete;


public:

    //- Runtime type information
    TypeName("dynamicOverFloaterMotionSolversFvMesh");


    // Constructors

    //- Construct from IOobject
    dynamicOverFloaterMotionSolversFvMesh(const IOobject& io,  const bool doInit=true);

    //- Destructor
    virtual ~dynamicOverFloaterMotionSolversFvMesh();

    // Member Functions

        //- Override ldu addressing
        virtual const lduAddressing& lduAddr() const
        {
            return oversetFvMeshBase::lduAddr();
        }

        //- Override ldu addressing
        virtual lduInterfacePtrsList interfaces() const
        {
            return oversetFvMeshBase::interfaces();
        }


        // Overset

            // Explicit interpolation

                virtual void interpolate(scalarField& psi) const
                {
                    Stencil::New(*this).interpolate<scalar>(*this, psi);
                }

                virtual void interpolate(vectorField& psi) const
                {
                    Stencil::New(*this).interpolate<vector>(*this, psi);
                }

                virtual void interpolate(sphericalTensorField& psi) const
                {
                    Stencil::New
                    (
                        *this
                    ).interpolate<sphericalTensor>(*this, psi);
                }

                virtual void interpolate(symmTensorField& psi) const
                {
                    Stencil::New(*this).interpolate<symmTensor>(*this, psi);
                }

                virtual void interpolate(tensorField& psi) const
                {
                    Stencil::New(*this).interpolate<tensor>(*this, psi);
                }

                virtual void interpolate(volScalarField& psi) const
                {
                    Stencil::New(*this).interpolate<volScalarField>(psi);
                }

                virtual void interpolate(volVectorField& psi) const
                {
                    Stencil::New(*this).interpolate<volVectorField>(psi);
                }

                virtual void interpolate(volSphericalTensorField& psi) const
                {
                    Stencil::New(*this).interpolate
                    <
                        volSphericalTensorField
                    >(psi);
                }

                virtual void interpolate(volSymmTensorField& psi) const
                {
                    Stencil::New(*this).interpolate<volSymmTensorField>(psi);
                }

                virtual void interpolate(volTensorField& psi) const
                {
                    Stencil::New(*this).interpolate<volTensorField>(psi);
                }


            // Implicit interpolation (matrix manipulation)

                //- Solve returning the solution statistics given convergence
                //- tolerance. Use the given solver controls
                virtual SolverPerformance<scalar> solve
                (
                    fvMatrix<scalar>& m,
                    const dictionary& dict
                ) const
                {
                    return oversetFvMeshBase::solveOverset<scalar>(m, dict);
                }

                //- Solve returning the solution statistics given convergence
                //- tolerance. Use the given solver controls
                virtual SolverPerformance<vector> solve
                (
                    fvMatrix<vector>& m,
                    const dictionary& dict
                ) const
                {
                    return oversetFvMeshBase::solveOverset<vector>(m, dict);
                }

                //- Solve returning the solution statistics given convergence
                //- tolerance. Use the given solver controls
                virtual SolverPerformance<symmTensor> solve
                (
                    fvMatrix<symmTensor>& m,
                    const dictionary& dict
                ) const
                {
                    return oversetFvMeshBase::solveOverset<symmTensor>
                    (m, dict);
                }

                //- Solve returning the solution statistics given convergence
                //- tolerance. Use the given solver controls
                virtual SolverPerformance<tensor> solve
                (
                    fvMatrix<tensor>& m,
                    const dictionary& dict
                ) const
                {
                    return oversetFvMeshBase::solveOverset<tensor>(m, dict);
                }

        //- Initialise all non-demand-driven data
        virtual bool init(const bool doInit);

        //- Initialise all non-demand-driven data. Used for initialising
        //  without default motionSolver (used by dynamicRefineFvMesh)
        virtual bool init(const bool doInit, const bool mandatory);

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Map all fields in time using given map. Triggered by topo change
        virtual void mapFields(const mapPolyMesh& mpm);

        //- Write using given format, version and compression
        virtual bool writeObject
        (
            IOstreamOption streamOpt,
            const bool writeOnProc
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif //DYNAMICOVERFLOATERMOTIONSOLVERSFVMESH_H

// ************************************************************************* //