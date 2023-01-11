/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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

#include "codedFixedValueFvPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
//{{{ begin codeInclude

//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = 80160fc90c13a9211e4c64e62666d4c79d4172c7
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void parabolicVelocity_80160fc90c13a9211e4c64e62666d4c79d4172c7(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatchVectorField,
    parabolicVelocityFixedValueFvPatchVectorField
);


const char* const parabolicVelocityFixedValueFvPatchVectorField::SHA1sum =
    "80160fc90c13a9211e4c64e62666d4c79d4172c7";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

parabolicVelocityFixedValueFvPatchVectorField::
parabolicVelocityFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF)
{
    if (false)
    {
        Info<<"construct parabolicVelocity sha1: 80160fc90c13a9211e4c64e62666d4c79d4172c7"
            " from patch/DimensionedField\n";
    }
}


parabolicVelocityFixedValueFvPatchVectorField::
parabolicVelocityFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict)
{
    if (false)
    {
        Info<<"construct parabolicVelocity sha1: 80160fc90c13a9211e4c64e62666d4c79d4172c7"
            " from patch/dictionary\n";
    }
}


parabolicVelocityFixedValueFvPatchVectorField::
parabolicVelocityFixedValueFvPatchVectorField
(
    const parabolicVelocityFixedValueFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct parabolicVelocity sha1: 80160fc90c13a9211e4c64e62666d4c79d4172c7"
            " from patch/DimensionedField/mapper\n";
    }
}


parabolicVelocityFixedValueFvPatchVectorField::
parabolicVelocityFixedValueFvPatchVectorField
(
    const parabolicVelocityFixedValueFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF)
{
    if (false)
    {
        Info<<"construct parabolicVelocity sha1: 80160fc90c13a9211e4c64e62666d4c79d4172c7 "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

parabolicVelocityFixedValueFvPatchVectorField::
~parabolicVelocityFixedValueFvPatchVectorField()
{
    if (false)
    {
        Info<<"destroy parabolicVelocity sha1: 80160fc90c13a9211e4c64e62666d4c79d4172c7\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void parabolicVelocityFixedValueFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        Info<<"updateCoeffs parabolicVelocity sha1: 80160fc90c13a9211e4c64e62666d4c79d4172c7\n";
    }

//{{{ begin code
    #line 32 "/home/julius/Documents/NumSim/NumSim/task3_openfoam/karmanvortexstreet/0/U/boundaryField/inlet"
const vectorField& Cf = patch().Cf();
            vectorField& field = *this;

            const scalar Umax = 1.0;
            forAll(Cf, faceI)
            {
                const scalar y = Cf[faceI][1];
                field[faceI] = vector(Umax*(1-((pow(y-0.5)/0.5,2))),0,0);
            }
//}}} end code

    this->fixedValueFvPatchField<vector>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

