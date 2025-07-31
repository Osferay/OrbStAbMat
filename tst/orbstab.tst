##### Testing the infinite orbit algorithm

gap> START_TEST("OrbStAbMat package: testall.tst");
gap> x:= X( Rationals );;
gap> pol := x^6 - x^5 + x^4 - x^3 + x^2 - x + 1;;
gap> a := CompanionMat(pol);;
gap> I := IdentityMat( 6 );;
gap> b := a-I;;
gap> c := a^4-a;;
gap> mats := [a,b,c];;
gap> v := [0,0,1,1,0,0];;
gap> w := v * a^2*b;;
gap> OrbitAbelianMatGroup( mats, v, w );
rec( 
  cong := [ [ 0, 0, 0, -1, 0, 1 ], [ 0, 0, 0, 1, -1, -1 ], 
      [ -1, 0, 0, -1, 1, 0 ], [ 1, -1, 0, 1, -1, 0 ], [ 0, 1, -1, -1, 1, 0 ], 
      [ 0, 0, 1, 0, -1, 0 ] ], 
  stab := 
    [ [ [ 1, 0, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0, 0 ], [ 0, 0, 1, 0, 0, 0 ], 
          [ 0, 0, 0, 1, 0, 0 ], [ 0, 0, 0, 0, 1, 0 ], [ 0, 0, 0, 0, 0, 1 ] ] 
     ] )
gap> StabilizerAbelianMatGroup( mats, v );
[ [ [ 1, 0, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0, 0 ], [ 0, 0, 1, 0, 0, 0 ], 
      [ 0, 0, 0, 1, 0, 0 ], [ 0, 0, 0, 0, 1, 0 ], [ 0, 0, 0, 0, 0, 1 ] ] ]
gap> STOP_TEST( "testall.tst", 10000 );