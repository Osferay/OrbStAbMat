##############################################################################
# Given a set of matrices that act on \Q^d, returns the radical module
##############################################################################
InstallGlobalFunction( RadicalModuleAbelianMatGroup,
    function( mats )

    local   U,
            V,
            hom,
            act,
            found;
    
    U   := IdentityMat( Length( mats[1] ) );
    V   := [];

    repeat
        hom   := NaturalHomomorphismBySemiEchelonBases( U, V );
        act   := List( mats, x -> InducedActionFactorByNHSEB( x, hom ) );
        #Find if in the induced action there is a non-diagonalizable element.
        found := RadicalElementBasisMatrixAlegebra( act );
        
        if IsMatrix(found) then
            #If there is such an element, compute the module U*found and loop.
            found := JordanDecomposition( found )[2];
            V     := List( hom.factor, x -> x*found );
            V     := List( V, x -> PreimagesRepresentativeByNHSEB( x, hom ) );
            V     := Concatenation( V, hom.kernel );
            V     := VectorspaceBasis( V );
        fi;
        
    until not IsMatrix( found );

    return V;
end );

#####################################################################################
# Given a set of matrices that act on \Q^d, returns the radical series
#####################################################################################
InstallGlobalFunction( RadicalSeriesAbelianMatGroup,
    function( mats )

    local   act,
            U,
            ser,
            rad,
            i,
            V,W,
            hom;

    U   := IdentityMat( Length( mats[1] ) );
    ser := [ U ];
    act := ShallowCopy( mats );

    repeat
        rad  := RadicalModuleAbelianMatGroup( act );
        hom  := NaturalHomomorphismBySemiEchelonBases( U, rad );
        act  := List( mats, x -> InducedActionFactorByNHSEB( x, hom ) );
        Add( ser, rad );
        U    := IdentityMat( Length( mats[1] ) );
        
    until Length( rad ) = 0;

    if DEBUG@ then 
    
        for i in [ 1..Length(ser)-1 ] do
            act  := InducedActionFactor( mats, U, ser[i], ser[i+1] ).act;

            #The induced action on the factor has to be semisimple and abelian.
            if not IsSemisimpleAbelianMatrixGroup( act ) then
                Error( "Incorrect radical series." );
            fi;

        od;
    fi;
    return ser;

end );

##############################################################
# Given a list of semisimple commutative matrices it returns
# the homogeneous splitting of the list, that is the nullspace
# of each factor of the minimal polynomial of the primitive
# element
##############################################################
HomogeneousSplitting := function( mats )
    local p,f,K,split,i;

    p := PrimitiveElementMatrixAlgebra( mats );
    f := Factors( p.minpol );

    split := [ ];
    for i in [1..Length(f)] do
        K := Value( f[i], p.b );
        Add( split, NullspaceMat( K ) );
    od;
    
    return split;

end;

################################################################################
# Given a list of semisimple commutative matrices it returns the homogeneous 
# series
################################################################################
InstallGlobalFunction( HomogeneousSemisimpleSeriesAbelianMatGroup, 
    function( mats )
    local   mact, 
            rad,
            U, V, W,
            act, 
            split, 
            i, j, 
            ser, 
            tmp, 
            f, hom;

    mact := ShallowCopy( mats );
    rad  := RadicalSeriesAbelianMatGroup( mact );
    U    := rad[1];
    ser  := [];

    for i in [ 1..Length(rad)-1 ] do
        Add( ser, rad[i] );
        act  := InducedActionFactor( mats, U, rad[i], rad[i+1] );
        split := HomogeneousSplitting( act.act );

        for j in [2..Length(split)] do
            tmp := SumOfSplitting( split, j );
            tmp := List( tmp, x -> PreimagesRepresentativeByNHSEB( x, act.hom ) );
            tmp := Concatenation( tmp, rad[i+1] );
            Add( ser, tmp );
        od;
    od;
    Add( ser, rad[i+1] );

    return ser;

end );

##############################################################
# Given a list of semisimple commutative matrices it returns
# the irreducible splitting of the list, that is the nullspace
# of each factor of the minimal polynomial of the primitive
# element divided by generators
##############################################################
IrreducibleSplitting := function( mats )
    local mact, f, p, K, split, i, V, tmp;

    mact := ShallowCopy( mats );
    p := PrimitiveElementMatrixAlgebra( mact );
    f := Factors(p.minpol);
    
    split := [ ];
    for i in [1..Length(f)] do

        K := Value( f[i], p.b );
        K := NullspaceMat( K );
        tmp := [];
        i := 1;
        repeat
            #Find a vector that is not in span of sub
            if not IsInSpan( K[i], tmp ) then
                V := K[i];
                V := InSpan( List( mact, x -> V*x ) , [ V ] );
                Add( split, V );
                tmp := Concatenation( tmp, V );
            fi;
            i := i + 1;
        until Rank( tmp ) = Rank( K );


    od;
    
    return split;

end;

################################################################################
# Given a list of semisimple commutative matrices it returns the irreducible 
# series
################################################################################
InstallGlobalFunction( IrreducibleSeriesAbelianMatGroup,
    function( mats )
    local mact, act, rad, U, V, W, iact, split, i, j, ser, tmp, f, hom;

    mact := ShallowCopy( mats );
    rad  := RadicalSeriesAbelianMatGroup( mact );
    U    := rad[1];
    ser  := [];

    for i in [ 1..Length(rad)-1 ] do
        Add( ser, rad[i] );
        act  := InducedActionFactor( mats, U, rad[i], rad[i+1] );
        split := IrreducibleSplitting( act.act );

        for j in [2..Length(split)] do
            tmp := SumOfSplitting( split, j );
            tmp := List( tmp, x -> PreimagesRepresentativeByNHSEB( x, act.hom ) );
            tmp := Concatenation( tmp, rad[i+1] );
            Add( ser, tmp );
        od;
    od;
    Add( ser, rad[i+1] );

    return ser;

end );

#####################################################
# Given a lattice computes the dual lattice of it
#####################################################
DualLattice := function( U )
    local A,D;

    A := ShallowCopy( U );
    if IsEmpty(A) then
        return A;
    else
        D := DenominatorMat( A );
        A := NullspaceIntMat( TransposedMat( D*A ) );
    fi;

    return A;

end;

#####################################################
# Given a lattice computes the dual lattice of the
# dual lattice of it
#####################################################
DualDualLattice := function( U )
    local A,D;
    
    A := ShallowCopy( U );
    if IsEmpty(A) then
        return A;
    else
        D := DenominatorMat( A );
        A := NullspaceIntMat( TransposedMat( D*A ) );

        if IsEmpty(A) then
            return IdentityMat( Length( U[1] ) );
        fi;

        return NullspaceIntMat( TransposedMat( A ) );
    fi;
    
end;

################################################################################
# Given a list of semisimple commutative matrices it returns an integral  
# irreducible block flag
################################################################################
InstallGlobalFunction( IntegralIrreducibleBlock,
    function( mats )
    local ser, dbg, i;

    ser := IrreducibleSeriesAbelianMatGroup( mats );
    ser := List( ser, DualDualLattice );

    if DEBUG@ then
        for i in [1.. Length(ser)-1] do
            dbg := SmithNormalFormIntegerMat( ser[i] );
            dbg := DiagonalOfMat( dbg );
            if ForAny( dbg, x -> x > 1 ) then
                Error( "The factors of the series are not pure." );
            fi;
        od;
    fi;

    return ser;

end );