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
    layerCells_.setSize(nWalk());
    layerY_.setSize(nWalk());
    forAll(layerCells_, l)
    {
        layerCells_[l] = labelList(patch().size(), -1);
        layerY_[l] = scalarField(patch().size(), 0.0);
    }

    label nShort = 0;

    forAll(faceCells, i)
    {
        label celli = faceCells[i];
        label entryFace = patch().start() + i;

        // Distance of a cell centre from the wall plane of face i
        auto yWall = [&](const label c)
        {
            return max(-((C[c] - Cf[i]) & nf[i]), small);
        };

        layerCells_[0][i] = celli;
        layerY_[0][i] = yWall(celli);
        sampleCells_[i] = celli;

        for (label layer = 2; layer <= nWalk(); layer++)
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
                if (layer <= sampleLayer_)
                {
                    nShort++;
                }
                break;
            }

            celli = own[best] == celli ? nei[best] : own[best];
            entryFace = best;

            layerCells_[layer - 1][i] = celli;
            layerY_[layer - 1][i] = yWall(celli);
            if (layer <= sampleLayer_)
            {
                sampleCells_[i] = celli;
            }
        }

        ySample_[i] = yWall(sampleCells_[i]);
    }

    sampled_ = true;

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


void Foam::nutSpaldingLESWallFunctionFvPatchScalarField::checkSampling() const
{
    if (!sampled_)
    {
        calcSampling();
    }
}


Foam::tmp<Foam::vectorField>
Foam::nutSpaldingLESWallFunctionFvPatchScalarField::sampleVelocity() const
{
    checkSampling();

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
    if (filtered() && Uf_.size() == patch().size())
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

    checkSampling();
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

    uTau_ = uTau();

    return max
    (
        scalar(0),
        sqr(uTau_)/(magGradU + rootVSmall) - nuw
    );
}


void Foam::nutSpaldingLESWallFunctionFvPatchScalarField::writeLocalEntries
(
    Ostream& os
) const
{
    nutWallFunctionFvPatchScalarField::writeLocalEntries(os);
    writeEntry(os, "sampleLayer", sampleLayer_);
    if (TfAuto_)
    {
        writeEntry(os, "Tf", word("auto"));
        writeEntry(os, "TfCrossings", TfCrossings_);
    }
    else
    {
        writeEntry(os, "Tf", Tf_);
    }
    writeEntry(os, "nutFloor", nutFloor_);
    if (filtered() && Uf_.size() == patch().size())
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
    Tf_(0),
    TfAuto_(!dict.found("Tf")),
    TfCrossings_(dict.lookupOrDefault<scalar>("TfCrossings", 4)),
    nutFloor_(dict.lookupOrDefault<bool>("nutFloor", false)),
    Uf_(),
    filterTimeIndex_(-1),
    sampled_(false),
    sampleCells_(),
    ySample_(),
    layerCells_(),
    layerY_(),
    uTau_()
{
    if (dict.found("Tf"))
    {
        ITstream& is = dict.lookup("Tf");
        const token t(is);
        if (t.isWord())
        {
            if (t.wordToken() != "auto")
            {
                FatalIOErrorInFunction(dict)
                    << "Tf " << t.wordToken() << ": a time scale or auto"
                    << exit(FatalIOError);
            }
            TfAuto_ = true;
        }
        else if (t.isNumber())
        {
            Tf_ = t.number();
        }
        else
        {
            FatalIOErrorInFunction(dict)
                << "Tf: a time scale or auto" << exit(FatalIOError);
        }
    }

    if (TfCrossings_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "TfCrossings " << TfCrossings_ << " <= 0" << exit(FatalIOError);
    }

    if (sampleLayer_ < 1)
    {
        FatalIOErrorInFunction(dict)
            << "sampleLayer " << sampleLayer_ << " < 1" << exit(FatalIOError);
    }

    if (filtered() && dict.found("Uf"))
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
    TfAuto_(ptf.TfAuto_),
    TfCrossings_(ptf.TfCrossings_),
    nutFloor_(ptf.nutFloor_),
    Uf_(),
    filterTimeIndex_(-1),
    sampled_(false),
    sampleCells_(),
    ySample_(),
    layerCells_(),
    layerY_(),
    uTau_()
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
    TfAuto_(wfpsf.TfAuto_),
    TfCrossings_(wfpsf.TfCrossings_),
    nutFloor_(wfpsf.nutFloor_),
    Uf_(wfpsf.Uf_),
    filterTimeIndex_(wfpsf.filterTimeIndex_),
    sampled_(false),
    sampleCells_(),
    ySample_(),
    layerCells_(),
    layerY_(),
    uTau_()
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

    checkSampling();

    return y*uTau()/tnuw();
}


Foam::label Foam::nutSpaldingLESWallFunctionFvPatchScalarField::limitNut
(
    scalarField& nut
) const
{
    checkSampling();

    if (!nutFloor_ || uTau_.size() != patch().size())
    {
        return 0;
    }

    const momentumTransportModel& turbModel =
        db().lookupType<momentumTransportModel>(internalField().group());
    const tmp<scalarField> tnuw = turbModel.nu(patch().index());
    const scalarField& nuw = tnuw();

    label n = 0;
    forAll(layerCells_, l)
    {
        const labelList& cells = layerCells_[l];
        const scalarField& y = layerY_[l];

        forAll(cells, i)
        {
            if (cells[i] < 0)
            {
                continue;
            }

            const scalar nutLaw = nuw[i]*wallLaw::spaldingNutByNu
            (
                y[i]*uTau_[i]/nuw[i], kappa_, E_
            );

            if (nutLaw > nut[cells[i]])
            {
                nut[cells[i]] = nutLaw;
                n++;
            }
        }
    }

    return n;
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

    sampled_ = false;
    sampleCells_.clear();
    ySample_.clear();
    layerCells_.clear();
    layerY_.clear();
    uTau_.clear();
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
    sampled_ = false;
    sampleCells_.clear();
    ySample_.clear();
    layerCells_.clear();
    layerY_.clear();
    uTau_.clear();
}


void Foam::nutSpaldingLESWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // The filter advances once per time step, with the velocity at the first
    // evaluation of the step (the solution of the previous step)
    if (filtered() && db().time().timeIndex() != filterTimeIndex_)
    {
        const vectorField Us(sampleVelocity());
        const scalar deltaT = db().time().deltaTValue();

        if (Uf_.size() != patch().size())
        {
            Uf_ = Us;
        }
        else if (TfAuto_)
        {
            // Tf = TfCrossings times the time the flow takes to cross the
            // wall cell, per face
            const scalarField& magSf = patch().magSf();
            forAll(Uf_, i)
            {
                const scalar U = max(mag(Us[i]), mag(Uf_[i]));
                const scalar eps = min
                (
                    deltaT*U/(TfCrossings_*sqrt(magSf[i]) + vSmall),
                    scalar(1)
                );
                Uf_[i] = (1 - eps)*Uf_[i] + eps*Us[i];
            }
        }
        else
        {
            const scalar eps = min(deltaT/Tf_, scalar(1));
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
