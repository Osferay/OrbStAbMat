###########################################################################
# Given a list of matrices, its relative order and an vector v returns the
# orbit stabilizer via the operation op
###########################################################################
MatrixOrbitStabilizer := function( mats, rels, v, op )
    local   orbit,
            dict,
            trans,
            trels,
            tmats,
            tword,
            stab,
            word,
            i,f,j,n,t,s,k;

    # set up
    orbit := [v];
    dict := NewDictionary(v, true);
    AddDictionary(dict, v, 1);
    trans := [];
    trels := [];
    tmats := [];
    tword := [];
    stab  := [];
    word  := [];

    # construct orbit and stabilizer
    for i in Reversed( [1..Length(mats)] ) do

        # get new point
        f := op( v, mats[i] );
        j := LookupDictionary( dict, f );

        # if it is new, add all blocks
        n := orbit;
        t := [];
        s := 1;
        while IsBool( j ) do
            n := List( n, x -> op( x, mats[i] ) );
            Append( t, n );
            j := LookupDictionary( dict, op( n[1], mats[i] ) );
            s := s + 1;
        od;

        # add to orbit
        for k in [1..Length(t)] do
            AddDictionary( dict, t[k], Length(orbit) + k );
        od;
        Append( orbit, t );
        
        # add to transversal
        if s > 1 then
            Add( trans, mats[i]^-1 );
            Add( trels, s );
            Add( tword, i );
            Add( tmats, [ i, s ] );
        fi;
        
        # compute stabiliser element
        if rels[i] = infinity or s < rels[i] then
            if j = 1 then
                Add( stab, mats[i]^s );
                Add( word, [[i,s]] );
            else
                t := TransversalInverse(j, trels);
                Add( stab, mats[i]^s * SubsWord( t, trans ) );
                Add( word, Concatenation( [[i,s]], Translate( t, tword )));
            fi;
        fi;
    od;

    # return orbit and stabilizer
    return rec( orbit := orbit,
                trels := trels,
                trans := trans,
                tmats := tmats,
                stab  := Reversed(stab),
                word  := Reversed(word) );
end ;

StabilizerFiniteAbelianMatGroup := function( mats, v )

    local   rels,
            stab;

    rels := List( mats, Order );

    stab := MatrixOrbitStabilizer( mats, rels, v, OnRight ).stab;

    if DEBUG@ then
        if ForAny( stab, x -> v*x <> v ) then
            Error( "Wrong stabilizer." );
        fi;
    fi;

    return stab;
end;

OrbitFiniteAbelianMatGroup := function( mats, v, w )

    local   rels,
            stab,
            fos,
            I,
            pos,
            y;

    rels := List( mats, Order );
    I    := IdentityMat( Length( mats[1] ) );
    fos  := MatrixOrbitStabilizer( mats, rels, v, OnRight );

    pos  := Position( fos.orbit, w );

    if IsBool(pos) then 
        return false;
    fi;

    stab := fos.stab;
    y    := TransversalElement( pos, fos, I );

    if DEBUG@ then
        if ForAny( stab, x -> v*x <> v ) then
            Error( "Wrong stabilizer." );
        fi;

        if v*y <> w then
            Error( "Wrong orbit element." );
        fi;
    fi;

    return rec( stab := stab, y := y );

end;