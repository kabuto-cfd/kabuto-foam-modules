/*---------------------------------------------------------------------------*\
  Kabuto CFD modules for OpenFOAM-14
  Copyright (C) 2026 Kabuto CFD
-------------------------------------------------------------------------------
License
    This file is part of kabuto-foam-modules, add-on libraries for OpenFOAM,
    distributed under the GNU General Public License version 3 or later.

Description
    Run-time selection of the RAS models of libkabutoMomentumTransportModels
    for incompressible flows.

\*---------------------------------------------------------------------------*/

#include "incompressibleMomentumTransportModels.H"
#include "makeIncompressibleMomentumTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "stressOmega.H"
makeRASModel(stressOmega);


// ************************************************************************* //
