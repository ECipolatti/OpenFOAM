/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicAMILduInterfaceField(),
    coupledFvPatchField<Type>(p, iF),
    cyclicAMIPatch_(refCast<const cyclicAMIFvPatch>(p))
{}


template<class Type>
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    cyclicAMILduInterfaceField(),
    coupledFvPatchField<Type>(p, iF, dict, dict.found("value")),
    cyclicAMIPatch_(refCast<const cyclicAMIFvPatch>(p))
{
    if (!isA<cyclicAMIFvPatch>(p))
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }

    if (!dict.found("value") && this->coupled())
    {
        this->evaluate(Pstream::commsTypes::blocking);
    }
}


template<class Type>
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const cyclicAMIFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    cyclicAMILduInterfaceField(),
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    cyclicAMIPatch_(refCast<const cyclicAMIFvPatch>(p))
{
    if (!isA<cyclicAMIFvPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const cyclicAMIFvPatchField<Type>& ptf
)
:
    cyclicAMILduInterfaceField(),
    coupledFvPatchField<Type>(ptf),
    cyclicAMIPatch_(ptf.cyclicAMIPatch_)
{}


template<class Type>
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const cyclicAMIFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicAMILduInterfaceField(),
    coupledFvPatchField<Type>(ptf, iF),
    cyclicAMIPatch_(ptf.cyclicAMIPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::cyclicAMIFvPatchField<Type>::coupled() const
{
    return cyclicAMIPatch_.coupled();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::cyclicAMIFvPatchField<Type>::patchNeighbourField() const
{
//    const Field<Type>& iField = this->primitiveField();
//    const labelUList& nbrFaceCells =
//        cyclicAMIPatch_.cyclicAMIPatch().neighbPatch().faceCells();

//    Field<Type> pnf(iField, nbrFaceCells);

//    tmp<Field<Type>> tpnf;
//    if (cyclicAMIPatch_.applyLowWeightCorrection())
//    {
//        tpnf = cyclicAMIPatch_.interpolate(pnf, this->patchInternalField()());
//    }
//    else
//    {
//        tpnf = cyclicAMIPatch_.interpolate(pnf);
//    }

//    if (doTransform())
//    {
//        tpnf.ref() = transform(forwardT(), tpnf());
//    }

//    return tpnf;

//______________________________________________________________
    //me traigo los valores internos del campo
    //me traigo los nombres de todas las celdas de interes
    //pnf es cada celda con su field
    //creo tpnf que es el valor interpolado para la celda actual para actualizar result
    const Field<Type>& iField = this->primitiveField();

    const labelUList& nbrFaceCells =
            cyclicAMIPatch_.cyclicAMIPatch().neighbPatch().faceCells();

    Field<Type> pnf(iField, nbrFaceCells);
    tmp<Field<Type>> tpnf(new Field<Type>(this->size()));

    Field<Type>& result = tpnf.ref();

    vectorField SURF = cyclicAMIPatch().neighbPatch().Sf();
    scalarField magSURF = cyclicAMIPatch().neighbPatch().magSf();

    //Depende el tipo hago distintos ponderados----------

    if(isType<scalarField>(pnf)){

        Type sumaPnf = sum(pnf)*0;
        sumaPnf = sum(magSURF * pnf);
        //        forAll(pnf, index){
        //            sumaPnf += pnf[index] * magSURF[index];
        //        }
        scalar totalArea = sum(magSURF);
        result = sumaPnf/totalArea;
    }
//    if(isType<vectorField>(pnf)){
//        Info<<"pnf ES VECTORIAL--------"<<endl;
//        vector summ(0,0,0);// = vector::zero;
//        //Redefino pnf


//        //        vectorField=0;
//        //        vectorField proyection = SURF ^ pnf; //Cross Product

//        //        vectorField sumaPnf(0,0,0);// = sum(pnf)*0;
//        ////        summ = sum(proyection * pnf);
//        //        sumaPnf = sum(pnf);

//        //Test cross product, result of it (-1 11 -7)
//        vector S(1,2,3);
//        vector T(4,1,1);
//        vector R= S^T;
//        Info<<"S:"<<S<<endl;
//        Info<<"T:"<<T<<endl;
//        Info<<"R:"<<R<<endl;
//        Info<<"INFORMACION"<<endl;

//        //        result = summ / sumaPnf;
//    }
    if(isType<tensorField>(pnf)){
        Info<<"pnf es TENSORIAL--------"<<endl;
        Info<<"Estoy en un caso no soportado actualmente"<<endl;
    }

    return tpnf;
}
////-------------------------------------------------------------------
//    //me traigo los valores internos del campo
//    //me traigo los nombres de todas las celdas de interes
//    //pnf es cada celda con su field
//    //creo tpnf que es el valor interpolado para la celda actual para actualizar result
//    const Field<Type>& iField = this->primitiveField();
//    const labelUList& nbrFaceCells =
//            cyclicAMIPatch_.cyclicAMIPatch().neighbPatch().faceCells();


//    tmp<Field<Type>> tpnf(new Field<Type>(this->size()));



//    vectorField SURF = cyclicAMIPatch().neighbPatch().Sf();
//    scalarField magSURF = cyclicAMIPatch().neighbPatch().magSf();

//    //Depende el tipo hago distintos ponderados----------
////    if(isType<scalarField>(iField)){
////        Info<<"pnf ES ESCALAR--------"<<endl;
////        Field<Type> pnf(iField, nbrFaceCells);
////        Type sumaPnf = sum(pnf)*0;
////        sumaPnf = sum(magSURF * pnf);
//////        forAll(pnf, index){

//////            sumaPnf += pnf[index] * magSURF[index];
//////        }
////        Info<<"sumaPnf:------------"<<sumaPnf<<endl;
////        scalar totalArea = sum(magSURF);

////        Info<<"totalArea: "<<totalArea<<endl;
//////      Info<<"DIVISION: sumaPnf/totalArea = "<<sumaPnf/totalArea<<endl;

////        result = sumaPnf/totalArea;

////    }
//    if(isType<vectorField>(iField)){

//        Info<<"pnf ES VECTORIAL--------"<<endl;
//        vectorField result = tpnf.ref();
//        vectorField pnf(iField, nbrFaceCells);
//        vectorField summ=sum(pnf)*0;//(0,0,0);// = vector::zero;
//        //Redefino pnf

//        scalarField proyection = SURF ^ pnf; //Cross Product

//        vector sumaPnf(0,0,0);// = sum(pnf)*0;
//        summ = sum(proyection * pnf);
//        sumaPnf = sum(pnf);

//        result = summ / sumaPnf;
//    }


//    return tpnf;
//}


//Foam::tmp<Foam::scalarField>
//Foam::cyclicAMIFvPatchField<scalarField>::patchNeighbourField() const
//{
//    const scalarField & iField = this->primitiveField();

//    const labelUList& nbrFaceCells =
//            cyclicAMIPatch_.cyclicAMIPatch().neighbPatch().faceCells();

//    scalarField pnf(iField, nbrFaceCells);
//    tmp<scalarField> tpnf(new scalarField(this->size()));

//    scalarField & result = tpnf.ref();

//    vectorField SURF = cyclicAMIPatch().neighbPatch().Sf();
//    scalarField magSURF = cyclicAMIPatch().neighbPatch().magSf();

//    Type sumaPnf = sum(pnf)*0;
//    sumaPnf = sum(magSURF * pnf);

//    scalar totalArea = sum(magSURF);
//    result = sumaPnf/totalArea;

//    return tpnf;
//}



template<class Type>
const Foam::cyclicAMIFvPatchField<Type>&
Foam::cyclicAMIFvPatchField<Type>::neighbourPatchField() const
{
    const GeometricField<Type, fvPatchField, volMesh>& fld =
        static_cast<const GeometricField<Type, fvPatchField, volMesh>&>
        (
            this->primitiveField()
        );

    return refCast<const cyclicAMIFvPatchField<Type>>
    (
        fld.boundaryField()[cyclicAMIPatch_.neighbPatchID()]
    );
}


template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::updateInterfaceMatrix
(
    scalarField& result,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{

    const labelUList& nbrFaceCells =
        cyclicAMIPatch_.cyclicAMIPatch().neighbPatch().faceCells();

    scalarField pnf(psiInternal, nbrFaceCells);

    // Transform according to the transformation tensors
    transformCoupleField(pnf, cmpt);

    if (cyclicAMIPatch_.applyLowWeightCorrection())
    {
        scalarField pif(psiInternal, cyclicAMIPatch_.faceCells());
        pnf = cyclicAMIPatch_.interpolate(pnf, pif);
    }
    else
    {
        pnf = cyclicAMIPatch_.interpolate(pnf);
    }

    // Multiply the field by coefficients and add into the result
    const labelUList& faceCells = cyclicAMIPatch_.faceCells();

    forAll(faceCells, elemI)
    {
        result[faceCells[elemI]] -= coeffs[elemI]*pnf[elemI];
    }



       //__________________HECHO POR SANTIAGO_________

//    const labelUList& nbrFaceCells =
//        cyclicAMIPatch_.cyclicAMIPatch().neighbPatch().faceCells();

//    scalarField pnf(psiInternal, nbrFaceCells);

//    // Transform according to the transformation tensors
//    //transformCoupleField(pnf, cmpt);

//    //pnf = cyclicAMIPatch_.interpolate(pnf);

//    // Multiply the field by coefficients and add into the result
//    const labelUList& faceCells = cyclicAMIPatch_.faceCells();

//    scalarField tpnf(faceCells.size(),sum(pnf)/pnf.size());


//    forAll(faceCells, elemI)
//    {
//        result[faceCells[elemI]] -= coeffs[elemI]*tpnf[elemI];
//    }
}


template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::updateInterfaceMatrix
(
    Field<Type>& result,
    const Field<Type>& psiInternal,
    const scalarField& coeffs,
    const Pstream::commsTypes
) const
{

    const labelUList& nbrFaceCells =
        cyclicAMIPatch_.cyclicAMIPatch().neighbPatch().faceCells();

    Field<Type> pnf(psiInternal, nbrFaceCells);

    // Transform according to the transformation tensors
    transformCoupleField(pnf);

    if (cyclicAMIPatch_.applyLowWeightCorrection())
    {
        Field<Type> pif(psiInternal, cyclicAMIPatch_.faceCells());
        pnf = cyclicAMIPatch_.interpolate(pnf, pif);
    }
    else
    {
        pnf = cyclicAMIPatch_.interpolate(pnf);
    }

    // Multiply the field by coefficients and add into the result
    const labelUList& faceCells = cyclicAMIPatch_.faceCells();

    forAll(faceCells, elemI)
    {
        result[faceCells[elemI]] -= coeffs[elemI]*pnf[elemI];
    }


  //__________________HECHO POR SANTIAGO_________

//    const labelUList& nbrFaceCells =
//        cyclicAMIPatch_.cyclicAMIPatch().neighbPatch().faceCells();

//    Field<Type> pnf(psiInternal, nbrFaceCells);
//    //pnf se tiene que "transformar" de neighb a owner (3-->1 o vv).
//    //pnf = cyclicAMIPatch_.interpolate(pnf);

//    // Multiply the field by coefficients and add into the result
//    const labelUList& faceCells = cyclicAMIPatch_.faceCells();

//    Field<Type>tpnf(faceCells.size(),sum(pnf)/pnf.size());

//    forAll(faceCells, elemI)
//    {
//        result[faceCells[elemI]] -= coeffs[elemI]*tpnf[elemI];
//    }


}


template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    this->writeEntry("value", os);
}


// ************************************************************************* //
