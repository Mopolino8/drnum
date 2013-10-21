/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

Application
    sonicFoam

Description
    Transient solver for trans-sonic/supersonic, laminar or turbulent flow
    of a compressible gas.

\*---------------------------------------------------------------------------*/

#include <list>
#include <vector>

#include "fvCFD.H"
#include "psiThermo.H"
#include "turbulenceModel.H"
#include "pimpleControl.H"

#include "mpicommunicator.h"
#include "externalexchangelist.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


bool inBox(double x1, double y1, double z1, double x2, double y2, double z2, double x, double y, double z)
{
  if (x < x1 || x > x2) return false;
  if (y < y1 || y > y2) return false;
  if (z < z1 || z > z2) return false;
  return true;
}


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    pimpleControl pimple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "drnumInit.H"
    #include "drnumSync.H"
    #include "drnumOverwrite.H"
    
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;              

        #include "readTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        #include "rhoEqn.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"
            #include "drnumOverwrite.H"
            #include "EEqn.H"
            #include "drnumOverwrite.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
                #include "drnumOverwrite.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        rho = thermo.rho();

        #include "drnumSync.H"
        #include "drnumOverwrite.H"
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    #include "drnumStop.H"
    
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
