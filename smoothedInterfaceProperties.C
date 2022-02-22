/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "smoothedInterfaceProperties.H"
//#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"
#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvc.H"
// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //

const Foam::scalar Foam::smoothedInterfaceProperties::convertToRad =
    Foam::constant::mathematical::pi/180.0;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Correction for the boundary condition on the unit normal nHat on
// walls to produce the correct contact angle.

// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the wall.

void Foam::smoothedInterfaceProperties::correctContactAngle
(
    surfaceVectorField::Boundary& nHatb,
    const surfaceVectorField::Boundary& gradAlphaf
) const
{
    const fvMesh& mesh = alpha1_.mesh();
    const volScalarField::Boundary& abf = alpha1_.boundaryField();

    const fvBoundaryMesh& boundary = mesh.boundary();

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleTwoPhaseFvPatchScalarField>(abf[patchi]))
        {
            alphaContactAngleTwoPhaseFvPatchScalarField& acap =
                const_cast<alphaContactAngleTwoPhaseFvPatchScalarField&>
                (
                    refCast<const alphaContactAngleTwoPhaseFvPatchScalarField>
                    (
                        abf[patchi]
                    )
                );

            fvsPatchVectorField& nHatp = nHatb[patchi];
            const scalarField theta
            (
                convertToRad*acap.theta(U_.boundaryField()[patchi], nHatp)
            );

            const vectorField nf
            (
                boundary[patchi].nf()
            );

            // Reset nHatp to correspond to the contact angle

            const scalarField a12(nHatp & nf);
            const scalarField b1(cos(theta));

            scalarField b2(nHatp.size());
            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            const scalarField det(1.0 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nHatp = a*nf + b*nHatp;
            nHatp /= (mag(nHatp) + deltaN_.value());

            acap.gradient() = (nf & nHatp)*mag(gradAlphaf[patchi]);
            acap.evaluate();
        }
    }
}

void Foam::smoothedInterfaceProperties::smoothen
(
    volScalarField& smooth_func
) 
{
    const fvMesh& mesh = smooth_func.mesh();
    const surfaceVectorField& Sf = mesh.Sf();

    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();

    for(int iter = 0; iter < smoothItr_; iter++)
    {
    	scalarField smooth_cal(mesh.nCells(),scalar(0));

    	scalarField sum_area(mesh.nCells(),scalar(0));

		surfaceScalarField smoothF = fvc::interpolate(smooth_func);

		for(int facei = 0; facei < nei.size(); facei++) //KVA note: should be own???
		{
			smooth_cal[own[facei]] += smoothF[facei]*mag(Sf[facei]);
			sum_area[own[facei]] += mag(Sf[facei]);
		}

		forAll(nei,facei)
		{
			smooth_cal[nei[facei]] += smoothF[facei]*mag(Sf[facei]);
			sum_area[nei[facei]] += mag(Sf[facei]);
		}

		forAll(mesh.boundary(), patchi)
		{
			const UList<label>& pFaceCells = mesh.boundary()[patchi].faceCells();

			const fvsPatchScalarField& pssf = smoothF.boundaryField()[patchi];

			forAll(mesh.boundary()[patchi], facei)
			{
			   smooth_cal[pFaceCells[facei]] += pssf[facei]*mag(Sf[facei]);
			   sum_area[pFaceCells[facei]] += mag(Sf[facei]);
			}
		}

		forAll(mesh.cells(),celli)
		{
			smooth_func[celli] = smooth_cal[celli]/sum_area[celli];
		}

		smooth_func.correctBoundaryConditions();

    }
}

void Foam::smoothedInterfaceProperties::calculateK()
{
    const fvMesh& mesh = alpha1_.mesh();
    const surfaceVectorField& Sf = mesh.Sf();

    smoothAlpha_ = alpha1_;

	// 1)
	// smoothen(smoothAlpha_);

	// 2)
    for (int i = 0; i < smoothItr_; ++i)
	{
    	//Lafaurie smooth function
        smoothAlpha_ = fvc::average(fvc::interpolate(smoothAlpha_));
	}


    // Cell gradient of alpha
    volVectorField gradAlpha(fvc::grad(smoothAlpha_));

    // Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));

    //gradAlphaf -=
    //    (mesh.Sf()/mesh.magSf())
    //   *(fvc::snGrad(alpha1_) - (mesh.Sf() & gradAlphaf)/mesh.magSf());

    // Face unit interface normal
    surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN_));
    // surfaceVectorField nHatfv
    // (
    //     (gradAlphaf + deltaN_*vector(0, 0, 1)
    //    *sign(gradAlphaf.component(vector::Z)))/(mag(gradAlphaf) + deltaN_)
    // );
    correctContactAngle(nHatfv.boundaryFieldRef(), gradAlphaf.boundaryField());

    // Face unit interface normal flux
    nHatf_ = nHatfv & Sf;

    // Simple expression for curvature
    K_ = -fvc::div(nHatf_);

	// Curvature smoothing
	if (kSmoothItr_ > 0)
	{
		volScalarField KsStar = K_;
		volScalarField Ks = K_;
		volScalarField kernel = 2.0*sqrt(mag(smoothAlpha_*(1.0-smoothAlpha_)));
		volScalarField w = sqrt(smoothAlpha_*(1.0 - smoothAlpha_) + 1e-6);
        for (int i = 0; i < kSmoothItr_; i++)
	    {
        	KsStar = fvc::average(fvc::interpolate(Ks*w))/fvc::average(fvc::interpolate(w));
	    	Ks = kernel*K_ + (1.0-kernel)*KsStar;
        }
		K_ = Ks;
	}

    // Complex expression for curvature.
    // Correction is formally zero but numerically non-zero.
    /*
    volVectorField nHat(gradAlpha/(mag(gradAlpha) + deltaN_));
    forAll(nHat.boundaryField(), patchi)
    {
        nHat.boundaryFieldRef()[patchi] = nHatfv.boundaryField()[patchi];
    }

    K_ = -fvc::div(nHatf_) + (nHat & fvc::grad(nHatfv) & nHat);
    */
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::smoothedInterfaceProperties::smoothedInterfaceProperties
(
    const volScalarField& alpha1,
    const volVectorField& U,
    const IOdictionary& dict
)
:
    transportPropertiesDict_(dict),
    cAlpha_
    (
        readScalar
        (
            alpha1.mesh().solverDict(alpha1.name()).lookup("cAlpha")
        )
    ),
    smoothItr_
    (
        readScalar
        (
            alpha1.mesh().solverDict(alpha1.name()).lookup("smoothItr")
        )
    ),
    kSmoothItr_
    (
        readScalar
        (
            alpha1.mesh().solverDict(alpha1.name()).lookup("kSmoothItr")
        )
    ),
    sigma_("sigma", dimensionSet(1, 0, -2, 0, 0), dict),

    deltaN_
    (
        "deltaN",
        1e-8/pow(average(alpha1.mesh().V()), 1.0/3.0)
    ),
    smoothAlpha_(alpha1),
    alpha1_(alpha1),
    U_(U),

    nHatf_
    (
        IOobject
        (
            "nHatf",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar("nHatf", dimArea, 0.0)
    ),

    K_
    (
        IOobject
        (
            "smoothedInterfaceProperties:K",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar("K", dimless/dimLength, 0.0)
    )
{
    calculateK();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField>
Foam::smoothedInterfaceProperties::surfaceTensionForce() const
{
    return fvc::interpolate(sigmaK())*fvc::snGrad(alpha1_);
}


Foam::tmp<Foam::volScalarField>
Foam::smoothedInterfaceProperties::nearInterface() const
{
    return pos(alpha1_ - 0.01)*pos(0.99 - alpha1_);
}


bool Foam::smoothedInterfaceProperties::read()
{
    alpha1_.mesh().solverDict(alpha1_.name()).lookup("cAlpha")     >> cAlpha_;
    alpha1_.mesh().solverDict(alpha1_.name()).lookup("smoothItr")  >> smoothItr_;
    alpha1_.mesh().solverDict(alpha1_.name()).lookup("kSmoothItr") >> kSmoothItr_;
    transportPropertiesDict_.lookup("sigma") >> sigma_;

    return true;
}


// ************************************************************************* //
