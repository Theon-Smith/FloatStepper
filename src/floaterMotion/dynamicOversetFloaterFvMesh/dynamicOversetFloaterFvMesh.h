//
// Created by theon on 2025/03/27.
//

#ifndef DYNAMICOVERSETFLOATERFVMESH_H
#define DYNAMICOVERSETFLOATERFVMESH_H

#include "dynamicFloaterMotionSolversFvMesh.H"
#include "oversetFvMeshBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class dynamicOversetFloaterFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicOversetFloaterFvMesh
:
    public dynamicFloaterMotionSolversFvMesh,
    public oversetFvMeshBase
{
    // Private Member Functions

        //- No copy construct
        dynamicOversetFloaterFvMesh(const dynamicOversetFloaterFvMesh&) = delete;

        //- No copy assignment
        void operator=(const dynamicOversetFloaterFvMesh&) = delete;


public:

    //- Runtime type information
    TypeName("dynamicOversetFloaterFvMesh");


    // Constructors

        //- Construct from IOobject
        dynamicOversetFloaterFvMesh(const IOobject& io,  const bool doInit=true);


    //- Destructor
    virtual ~dynamicOversetFloaterFvMesh();


    // Member Functions

        //- Override ldu addressing
        virtual const lduAddressing& lduAddr() const
        {
            return dynamicOversetFloaterFvMesh::lduAddr();
        }

        //- Override ldu addressing
        virtual lduInterfacePtrsList interfaces() const
        {
            return dynamicOversetFloaterFvMesh::interfaces();
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


        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

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

#endif //DYNAMICOVERSETFLOATERFVMESH_H
