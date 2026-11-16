/*---------------------------------------------------------------------------*\
  Kabuto CFD modules for OpenFOAM-14
  Copyright (C) 2026 Kabuto CFD
-------------------------------------------------------------------------------
License
    This file is part of kabuto-foam-modules, add-on libraries for OpenFOAM.

    It is free software: you can redistribute it and/or modify it under the
    terms of the GNU General Public License as published by the Free Software
    Foundation, either version 3 of the License, or (at your option) any later
    version.

    It is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
    details.

    You should have received a copy of the GNU General Public License along
    with it.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "nutSpaldingLESWallFunctionFvPatchScalarField.H"
#include "wallLaw.H"
#include "momentumTransportModel.H"
#include "fieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::nutSpaldingLESWallFunctionFvPatchScalarField::calcSampling() const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const cellList& cells = mesh.cells();
    const labelUList& own = mesh.faceOwner();
    const labelUList& nei = mesh.faceNeighbour();
    const vectorField& Sf = mesh.faceAreas();
    const vectorField& C = mesh.cellCentres();

    const labelUList& faceCells = patch().faceCells();
    const vectorField nf(patch().nf());
    const vectorField& Cf = patch().Cf();

    sampleCells_.setSize(patch().size());
    ySample_.setSize(patch().size());

    label nShort = 0;

    forAll(faceCells, i)
    {
        label celli = faceCells[i];
        label entryFace = patch().start() + i;

        for (label layer = 2; layer <= sampleLayer_; layer++)
        {
            // The face of celli, other than the one we came through, whose
            // outward normal is closest to the inward wall normal
            label best = -1;
            scalar bestCos = -great;
            forAll(cells[celli], cfi)
            {
                const label facei = cells[celli][cfi];
                if (facei == entryFace) continue;

                const vector n
                (
                    (own[facei] == celli ? 1 : -1)*Sf[facei]
                   /max(mag(Sf[facei]), vSmall)
                );
                const scalar c = -(n & nf[i]);
                if (c > bestCos)
                {
                    bestCos = c;
                    best = facei;
                }
            }

            if (best < 0 || !mesh.isInternalFace(best))
            {
                nShort++;
                break;
            }

            celli = own[best] == celli ? nei[best] : own[best];
            entryFace = best;
        }

        sampleCells_[i] = celli;
        ySample_[i] = max(-((C[celli] - Cf[i]) & nf[i]), small);
    }

    const label nFaces = returnReduce(patch().size(), sumOp<label>());
    reduce(nShort, sumOp<label>());

    if (nFaces > 0)
    {
        Info<< type() << " " << patch().name()
            << ": sample layer " << sampleLayer_
            << ", mean sample distance "
            << returnReduce(sum(ySample_), sumOp<scalar>())/nFaces
            << ", walk stopped early on " << nShort << " of " << nFaces
            << " faces" << endl;
    }
}


Foam::tmp<Foam::vectorField>
Foam::nutSpaldingLESWallFunctionFvPatchScalarField::sampleVelocity() const
{
    if (sampleCells_.size() != patch().size())
    {
        calcSampling();
    }

    const label patchi = patch().index();
    const momentumTransportModel& turbModel =
        db().lookupType<momentumTransportModel>(internalField().group());

    const volVectorField& U = turbModel.U();
    const fvPatchVectorField& Uw = U.boundaryField()[patchi];
    const vectorField nf(patch().nf());

    tmp<vectorField> tUs(new vectorField(patch().size()));
    vectorField& Us = tUs.ref();

    forAll(Us, i)
    {
        const vector Ur(U[sampleCells_[i]] - Uw[i]);
        Us[i] = Ur - (Ur & nf[i])*nf[i];
    }

    return tUs;
}


Foam::tmp<Foam::vectorField>
Foam::nutSpaldingLESWallFunctionFvPatchScalarField::inputVelocity() const
{
    if (Tf_ > 0 && Uf_.size() == patch().size())
    {
        return tmp<vectorField>(new vectorField(Uf_));
    }

    return sampleVelocity();
}


Foam::tmp<Foam::scalarField>
Foam::nutSpaldingLESWallFunctionFvPatchScalarField::uTau() const
{
    const label patchi = patch().index();
    const momentumTransportModel& turbModel =
        db().lookupType<momentumTransportModel>(internalField().group());

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const scalarField magUin(mag(inputVelocity()));

    tmp<scalarField> tuTau(new scalarField(patch().size(), 0.0));
    scalarField& ut = tuTau.ref();

    forAll(ut, facei)
    {
        ut[facei] = wallLaw::spaldingUTau
        (
            magUin[facei], ySample_[facei], nuw[facei], kappa_, E_
        );
    }

    return tuTau;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::nutSpaldingLESWallFunctionFvPatchScalarField::nut() const
{
    const label patchi = patch().index();
    const momentumTransportModel& turbModel =
        db().lookupType<momentumTransportModel>(internalField().group());

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magGradU(mag(Uw.snGrad()));
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    return max
    (
        scalar(0),
        sqr(uTau())/(magGradU + rootVSmall) - nuw
    );
}


void Foam::nutSpaldingLESWallFunctionFvPatchScalarField::writeLocalEntries
(
    Ostream& os
) const
{
    nutWallFunctionFvPatchScalarField::writeLocalEntries(os);
    writeEntry(os, "sampleLayer", sampleLayer_);
    writeEntry(os, "Tf", Tf_);
    if (Tf_ > 0 && Uf_.size() == patch().size())
    {
        writeEntry(os, "Uf", Uf_);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nutSpaldingLESWallFunctionFvPatchScalarField::
nutSpaldingLESWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, fvMesh>& iF,
    const dictionary& dict
)
:
    nutWallFunctionFvPatchScalarField(p, iF, dict),
    sampleLayer_(dict.lookupOrDefault<label>("sampleLayer", 1)),
    Tf_(dict.lookupOrDefault<scalar>("Tf", 0)),
    Uf_(),
    filterTimeIndex_(-1),
    sampleCells_(),
    ySample_()
{
    if (sampleLayer_ < 1)
    {
        FatalIOErrorInFunction(dict)
            << "sampleLayer " << sampleLayer_ << " < 1" << exit(FatalIOError);
    }

    if (Tf_ > 0 && dict.found("Uf"))
    {
        Uf_ = vectorField("Uf", dict, p.size());
    }
}


Foam::nutSpaldingLESWallFunctionFvPatchScalarField::
nutSpaldingLESWallFunctionFvPatchScalarField
(
    const nutSpaldingLESWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, fvMesh>& iF,
    const fieldMapper& mapper
)
:
    nutWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    sampleLayer_(ptf.sampleLayer_),
    Tf_(ptf.Tf_),
    Uf_(),
    filterTimeIndex_(-1),
    sampleCells_(),
    ySample_()
{
    if (ptf.Uf_.size() == ptf.patch().size())
    {
        Uf_.setSize(p.size());
        mapper(Uf_, ptf.Uf_);
    }
}


Foam::nutSpaldingLESWallFunctionFvPatchScalarField::
nutSpaldingLESWallFunctionFvPatchScalarField
(
    const nutSpaldingLESWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, fvMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(wfpsf, iF),
    sampleLayer_(wfpsf.sampleLayer_),
    Tf_(wfpsf.Tf_),
    Uf_(wfpsf.Uf_),
    filterTimeIndex_(wfpsf.filterTimeIndex_),
    sampleCells_(),
    ySample_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::nutSpaldingLESWallFunctionFvPatchScalarField::yPlus() const
{
    const label patchi = patch().index();
    const momentumTransportModel& turbModel =
        db().lookupType<momentumTransportModel>(internalField().group());

    const scalarField& y = turbModel.yb()[patchi];
    const tmp<scalarField> tnuw = turbModel.nu(patchi);

    if (sampleCells_.size() != patch().size())
    {
        calcSampling();
    }

    return y*uTau()/tnuw();
}


void Foam::nutSpaldingLESWallFunctionFvPatchScalarField::map
(
    const fvPatchScalarField& ptf,
    const fieldMapper& mapper
)
{
    nutWallFunctionFvPatchScalarField::map(ptf, mapper);

    const nutSpaldingLESWallFunctionFvPatchScalarField& tiptf =
        refCast<const nutSpaldingLESWallFunctionFvPatchScalarField>(ptf);

    if (tiptf.Uf_.size() == tiptf.patch().size())
    {
        Uf_.setSize(patch().size());
        mapper(Uf_, tiptf.Uf_);
    }
    else
    {
        Uf_.clear();
    }

    sampleCells_.clear();
    ySample_.clear();
}


void Foam::nutSpaldingLESWallFunctionFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    nutWallFunctionFvPatchScalarField::reset(ptf);

    const nutSpaldingLESWallFunctionFvPatchScalarField& tiptf =
        refCast<const nutSpaldingLESWallFunctionFvPatchScalarField>(ptf);

    Uf_ = tiptf.Uf_;
    sampleCells_.clear();
    ySample_.clear();
}


void Foam::nutSpaldingLESWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // The filter advances once per time step, with the velocity at the first
    // evaluation of the step (the solution of the previous step)
    if (Tf_ > 0 && db().time().timeIndex() != filterTimeIndex_)
    {
        const vectorField Us(sampleVelocity());

        if (Uf_.size() != patch().size())
        {
            Uf_ = Us;
        }
        else
        {
            const scalar eps = min(db().time().deltaTValue()/Tf_, scalar(1));
            Uf_ = (1 - eps)*Uf_ + eps*Us;
        }

        filterTimeIndex_ = db().time().timeIndex();
    }

    nutWallFunctionFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        nutSpaldingLESWallFunctionFvPatchScalarField
    );
}

// ************************************************************************* //
