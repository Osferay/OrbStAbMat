#####################################################################
# Induces the derivation der to the quotient U/V
#####################################################################
InducedDerivationQuotient := function( der, U, V )

    local   nat;

    nat := NaturalHomomorphismByLattices( U, V );

    return List( der, x -> ImageByNHLB( x, nat ) );

end;

#####################################################################
# Induces the action on U to the quotient V/W
#####################################################################
InducedActionFactor := function( mact, U, V, W )

    local   hom,
            iact,
            A, B;

    hom  := NaturalHomomorphismBySemiEchelonBases( U, V );
    iact  := List( mact, x -> InducedActionSubspaceByNHSEB( x, hom ) );

    A    := IdentityMat( Length( iact[1] ) );
    B    := List( W, x -> MemberBySemiEchelonBase( x, V ) );

    hom  := NaturalHomomorphismBySemiEchelonBases( A, B );
    iact := List( iact, x -> InducedActionFactorByNHSEB( x, hom ) );

    return rec( act := iact, hom := hom );
end;

#################################################################
# Given a list of matrices computes a basis of the matrix algebra
#################################################################
InstallGlobalFunction( BasisMatrixAlgebra,
    function( mact )

    local   B,      	#Basis matrix algebra
            act,        #Underlying action of mact
            tmp,        #temporal variable
            b,
            mult,       #Multiplication of elements
            d,          #Dimension of the action
            i;          #Bucle variable

    act := ShallowCopy( mact );
    d   := Length( mact[1] );
    act := List( mact, Flat );
    B   := [];
    tmp := ShallowCopy( act );
    
    while Length( tmp ) > 0 do
        b := tmp[1];
        Remove( tmp, 1 );
        mult := List( act, x -> MatrixMultFlat( b, x, d ) );
        
        for i in [1..Length( mult )] do
            if not IsInSpan( mult[i], B ) then
                Add( B, mult[i]);   
                Add( tmp, Last(B) );
            fi;
        od;
    od;
    TriangulizeMat(B);

    #Quick check
    if DEBUG@ = true then
        if B <> InSpan( act , B) then
            Error( "Wrong basis for the matrix algebra Q[G_p].");
        fi;
    fi;

    return List( B, x -> RecoverMatrixByVector( x, d ) );
 
end );

#################################################################################
# Given a list of matrices computes a non-diagonalizable element in the matrix 
# algebra or a basis of the matrix algebra if there is no such element
#################################################################################
RadicalElementBasisMatrixAlegebra := function( mact )
    local   d,
            act,
            B,
            tmp,
            mult,
            b,
            i;

    d := Length( mact[1] );
    act := ShallowCopy( mact );
    act := List( act, Flat );
    B := [];
    tmp := TriangulizedMat( act );
    while Length( tmp ) > 0 do
        b := tmp[1];
        Remove( tmp, 1 );
        mult := List( act, x -> MatrixMultFlat( b, x, d ) );
        for i in [1..Length( mult )] do
            if not IsInSpan( mult[i], B ) then
                Add( B, mult[i]);
                TriangulizeMat(B);
                if IsSemisimpleMatrix( RecoverMatrixByVector( Last(B), d ) ) then
                    # Continue spinning over the basis
                    Add( tmp, Last(B) );
                else
                    # We found a non semisimple element
                    return RecoverMatrixByVector( Last(B), d );
                fi;
            fi;
        od;
    od;
       
    return List( B, x -> RecoverMatrixByVector( x, d ) );

end ;

######################################################################
# Given a list of abelian semisimple matrices returns a primitive 
# element that generates the matrix algebra of the set.
######################################################################
PrimitiveElementMatrixAlgebra := function( mact )

    local   B,
            d,
            b,
            c,
            f,
            i;
    
    B := BasisMatrixAlgebra( mact );
    d := Length( B );

    #First test the generators

    i := 1;
    repeat
        b := B[i];
        f := MinimalPolynomial( b );
        i := i + 1;
    until Degree( f ) = d or i > Length(B);

    if i > Length(B) then
        repeat

            c := List( [1..d], x -> Random( Integers ) );
            b := c*B;
            f := MinimalPolynomial( Rationals, b );

        until Degree( f ) = d;
    fi;

    if DEBUG@ then
        if not IsInSpan( Flat(b), List( B, Flat) ) then
            Error( "Wrong primitive element of the matrix algebra." );
        fi;
    fi;
    
    return rec( b := b, minpol := f );
end;

######################################################################
# Given a list of abelian semisimple matrices returns a basis 
# of the relattion lattice.
######################################################################
RelationLatticeAbelianMat := function( mact )

    local   F,
            B,
            p,
            rel;

    B   := BasisMatrixAlgebra( mact );
    p   := PrimitiveElementMatrixAlgebra( mact );
    F   := FieldByMatrixBasisNC( B );
    SetPrimitiveElement( F, p.b );
    SetDefiningPolynomial( F, p.minpol );
    
    return RelationLatticeOfTFUnits( F, mact );

end;

######################################################################
# Given a list of abelian semisimple matrices and a matrix returns a 
# vector expressing if the given matrix is in the group generated by
# the given list or false if is not the case
######################################################################
MembershipRelationLatticeAbelianMat := function( mact, M )

    local   B,
            p,
            F,
            act,
            rel;

    B   := BasisMatrixAlgebra( mact );
    p   := PrimitiveElementMatrixAlgebra( mact );
    F   := FieldByMatrixBasisNC( B );
    SetPrimitiveElement( F, p.b );
    SetDefiningPolynomial( F, p.minpol );

    if not IsUnitOfNumberField( F, M ) then
        return false;
    fi;

    act  := ShallowCopy( mact );
    act  := Concatenation( [M] , act );

    rel  := RelationLatticeOfUnits( F, act )[1];

    if PositionNonZero( rel ) > 1 or AbsInt( rel[1] ) <> 1 then
        return false;
    else
        return -rel{ [2..Length(rel)] } * rel[1];
    fi;

end;