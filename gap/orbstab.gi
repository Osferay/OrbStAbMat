#####################################################################
# Constructs the affine action associated to der using act
#####################################################################
AssociatedAffineAction := function( act, der ) 
    local   aff,
            tmp,
            i,j,
            e,
            d;

    aff := [];
    d   := Length( act[1] );
    for i in [1..Length( act )] do
        tmp := [];
        for j in [1..d] do
            Add( tmp, Concatenation( ShallowCopy(act[i][j]), [0] ) );
        od;
        
        Add( tmp, Concatenation( ShallowCopy(der[i]), [1] ) );
        Add( aff, tmp );
    od;

    e    := ListWithIdenticalEntries( d, 0 );
    Add( e, 1 );

    return rec( affine := aff, e := e );
end;

####################################################################
# Support function that returns the kernel of the function gamma
# the returned matrices are in the respective quotient.
####################################################################
GammaFOS := function( mats, rels, der, derC )

    local   gamma,
            base,
            aff,
            fos;

    base  := ShallowCopy( derC );
    base  := VectorspaceBasis( base );
    # Construct the affine action of mats and der
    aff   := AssociatedAffineAction( mats, der );
    #Constuct gamma as a operation that returns a vector modulo derC
    gamma := function( mat, v )
        local res;
        #The lenght of the affine action is one more than the lenght
        #of the elements in base.
        res := OnRight( v, mat );
        res := VectorModLattice( res{[1..Length(res)-1]}, base );
        Add( res, 1 );
        return res;
    end;

    #Compute the finite orbstab of the affine action for e
    fos  := MatrixOrbitStabilizer( aff.affine, rels, aff.e, gamma );

    return fos;

end;

####################################################################
# Support function that returns the matrices of the kernel of the 
# function gamma translated to the whole group. Also checks if this
# matrices are in the kernel.
####################################################################
CheckGammaFOS := function( mats, fos, delta, derC, U, V )

    local   ker,
            tmp,
            I,
            word,
            str,
            base,
            der;

    ker  := [];
    I    := IdentityMat( Length (mats[1] ) );
    # Translate the words in the stab to mats
    for word in fos.word do
        tmp := I;
        for str in word do
            tmp := tmp * mats[ str[1] ]^str[2];
        od;
        Add( ker, tmp );
    od;
    der  := List( ker, delta );
    der  := InducedDerivationQuotient( der, U, V );

    #Check if the obtained elms are in derC

    if DEBUG@ then
        base := ShallowCopy( derC );
        base := VectorspaceBasis( base );
        tmp := List( der, x -> MemberBySemiEchelonBase( x, base ) );
        if ForAny( tmp, x -> IsBool( x ) ) then
            Error( "Wrong kernel of gamma.");
        fi;
    fi;

    return rec( ker := ker, der := der );
end;

############################################################
# Support function that returns the kernel of the function 
# delta
############################################################
ExtendGammaFOS := function( fos, C, derC )
    
    local   ext,
            tmp,
            i;
        
    ext := [];

    for i in [1..Length( fos.der ) ] do
        tmp := -1*SolutionIntMat( derC, fos.der[i] );
        Add( ext, fos.ker[i] * ( MappedVector( tmp, C ) ) );
    od;

    return ext;
    
end;

############################################################
# Support function that returns the checks if the matrices
# of the kernel of delta are correct.
############################################################
CheckExtendGammaFOS := function( mats, delta, U, V )

    local   der;
    
    if DEBUG@ then
        der := List( mats, delta );
        der := InducedDerivationQuotient( der, U, V );

        if ForAny( der, x -> x <> x*0 ) then
            Error( "Wrong kernel of delta." );
        fi;
    fi;
end;

##################################################################################
# Function to compute the stablizer of the given vector v under the action of the
# list of matrices mats
##################################################################################
InstallGlobalFunction( StabilizerAbelianMatGroup,
    function( mats, v )

    local   ser,
            stab,
            imat,
            I,
            delta,
            i,
            rel,
            ider,
            C,
            derC,  
            fos;

    if not IsAbelianMatrixAction( mats ) then
        Error( "The group G has to be an abelian matrix group.");
    fi;

    ser    := IntegralIrreducibleBlock( mats );
    stab   := mats;
    I      := ser[1];

    delta := function( M )
        return v*(M-I);
    end;

    for i in [ 1..Length(ser)-1 ] do
    
        ider := List( stab, delta );
        ider := InducedDerivationQuotient( ider, ser[i], ser[i+1] );
        
        if ForAll( ider, x -> x = x*0 ) then
            #Do nothing as stab is preserved
        else

            imat := InducedActionFactor( stab, ser[1], ser[i], ser[i+1] ).act;
            rel  := RelationLatticeAbelianMat( imat );
            C    := List( rel, x -> MappedVector( x, stab ) );
            
            derC := List( C, delta );
            derC := InducedDerivationQuotient( derC, ser[i], ser[i+1] );
            
            if ForAll( derC, x -> x = x*0 ) then
                stab := C;
            else
                fos  := GammaFOS( imat, List( stab, Order ), ider, derC );
                fos  := CheckGammaFOS( stab, fos, delta, derC, ser[i], ser[i+1] );
                stab := ExtendGammaFOS( fos, C, derC );
                CheckExtendGammaFOS( stab, delta, ser[i], ser[i+1] );
            fi;

        fi;
    od;

    if DEBUG@ then
        if ForAny( stab, x -> v*x <> v ) then
            Error( "Wrong stabilizer" );
        fi;
    fi;

    return stab;
end );

#################################################################
# Given a list of matrices computes returns the exponents of the
# elements of a basis of the matrix algebra and the basis
#################################################################
ExponentsBasisMatrixAlgebra := function( mats )

    local   B,
            I,
            d,
            i,
            exp,
            act,
            coef,
            new;

    d   := Length( mats[1] );
    I   := IdentityMat( d );
    B   := [];
    exp := [];

    #First test the gens

    for i in [ 1..Length(mats) ] do
        new := Flat( mats[i] - I );
        if not IsInSpan( new, B ) then
            Add( B, new );
            coef := ListWithIdenticalEntries( Length(mats), 0 );
            coef[i] := 1;
            Add( exp, coef );
        fi;
    od;
    
    while Length( B ) < d do
        coef := List( [1..Length(mats)], x -> Random( [0..5] ) );
        new  := MappedVector( coef, mats ) - I;
        new  := Flat( new );
        if not IsInSpan( new, B ) then
            Add( B, new );
            Add( exp, coef );
        fi;
    od;
    
    if DEBUG@ then
        act := List( mats , Flat );
        if not ForAll( act, x -> IsInSpan( x, B ) ) then
            Error( "Wrong basis for the induced action." );
        fi;
    fi;

    B := List( B, x -> RecoverMatrixByVector( x, d ) );

    return rec( base := B, exp := exp );

end;

MembershipDelta := function( mats, imat, u, delta )

    local   I,
            B,
            exp,
            g;

    I   := IdentityMat( Length( imat[1] ) );
    # Compute a basis of Q[G] of the form g_i-1
    B   := ExponentsBasisMatrixAlgebra( imat );
    # Translate this basis to elements of G
    exp := List( B.exp, x -> MappedVector( x, mats ) );
    # Compute a basis of Q^d with vectors delta(g_i)
    exp := List( exp, delta );
    # Express u in the last basis
    exp := SolutionMat( exp, u );
    # Find the corresponding element in G
    g   := exp*B.base + I;
    # Test membership in G
    exp := MembershipRelationLatticeAbelianMat( imat, g );

    if IsBool(g) then 
        return rec( membership := false ); 
    else 
        g := MappedVector( exp, mats );

    fi;

    return rec( membership := true, y := g );
end;

MembershipDeltaC := function( u, U, V, derC)

    local   vec,
            sol;

    vec   := InducedDerivationQuotient( [u], U, V )[1];
    sol   := SolutionIntMat( derC, vec );

    if IsBool(sol) then
        return rec( membership := false );
    fi;

    return rec( membership := true, y := sol, vec := vec );

end;

MembershipGammaFOS := function( u, mats, rels, der, derC )

    local   I,
            fos,
            vec,
            pos,
            elm;

    I   := IdentityMat( Length( mats[1] ) );
    fos := GammaFOS( mats, rels, der, derC );
    vec := ShallowCopy(u);
    Add( u, 1 );
    
    pos := Position( fos.orbit, vec );

    if IsBool( pos ) then
        return rec( membership := false );
    fi;

    elm := TransversalElement( pos, fos, I );

    return rec( membership := true, y := elm, fos := fos );

end;

ExtendMembershipGammaFOS := function( fos, C, derC )
    local   ext;
    
    ext := ExtendGammaFOS( fos, C, derC );

end;

##################################################################################
# Function to compute the whether if the two given vectors v and w are in the same
# orbit under the action of the list of matrices mats. Also returns stab of v
##################################################################################
InstallGlobalFunction( OrbitAbelianMatGroup,
    function( mats, v, w )

    local   ser,
            stab,
            imat,
            I,
            mem,
            y,
            u,
            delta,
            i,
            rel,
            ider,
            C,
            derC,  
            fos;

    if not IsAbelianMatrixAction( mats ) then
        Error( "The group G has to be an abelian matrix group.");
    fi;

    ser    := IntegralIrreducibleBlock( mats );
    stab   := mats;
    I      := ser[1];
    y      := I;

    for i in [ 1..Length(ser)-1 ] do

        u := v*y;

        delta := function( M )
            return u*(M-I);
        end;
    
        ider := List( stab, delta );
        ider := InducedDerivationQuotient( ider, ser[i], ser[i+1] );
        
        if ForAll( ider, x -> x = x*0 ) then
            #u and w are in the same orbit iff w-u \in ser[i+1]
            if IsBool( MemberBySemiEchelonBase( w-u, ser[i+1] ) ) then return false; fi;
        else

            imat := InducedActionFactor( stab, ser[1], ser[i], ser[i+1] ).act;
            rel  := RelationLatticeAbelianMat( imat );
            C    := List( rel, x -> MappedVector( x, stab ) );
            
            derC := List( C, delta );
            derC := InducedDerivationQuotient( derC, ser[i], ser[i+1] );
            
            if ForAll( derC, x -> x = x*0 ) then
                mem  := MembershipDelta( stab, imat, w-u, delta );
                if mem.membership = false then return false; fi;

                y    := y*mem.y;
                stab := C;
            else
                mem  := MembershipDeltaC( w-u, ser[i], ser[i+1], derC );
                
                if mem.membership = false then 
                    fos  := MembershipGammaFOS( mem.vec, imat, List( stab, Order ), ider, derC );

                    if fos.membership = false then
                        return false;
                    fi;

                    fos  := CheckGammaFOS( stab, fos.fos, delta, derC, ser[i], ser[i+1] );
                    stab := ExtendMembershipGammaFOS( fos, C, derC );
                    CheckExtendGammaFOS( stab, delta, ser[i], ser[i+1] );
                else
                    y := y * ( MappedVector( mem.y, C ) );
                    
                    fos  := GammaFOS( imat, List( stab, Order ), ider, derC );
                    fos  := CheckGammaFOS( stab, fos, delta, derC, ser[i], ser[i+1] );
                    stab := ExtendGammaFOS( fos, C, derC );
                    CheckExtendGammaFOS( stab, delta, ser[i], ser[i+1] );
                fi;
            fi;

        fi;
    od;

    if DEBUG@ then
        if ForAny( stab, x -> v*x <> v ) then
            Error( "Wrong stabilizer." );
        fi;

        if v*y <> w then
            Error( "Wrong orbit element." );
        fi;
    fi;

    return rec( stab := stab, cong := y );
end);