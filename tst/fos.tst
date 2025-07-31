##### Testing the finite orbit stablizer

gap> START_TEST("OrbStAbMat package: testall.tst");
gap> F := GF(7);;
gap> a := [[1,0,0], [0,6,0], [0,0,6]]*One(F);;
gap> b := [[6,0,0], [0,1,0], [0,0,1]]*One(F);;
gap> c := [[3,0,0], [0,2,0], [0,0,5]]*One(F);;
gap> mats := [a,b,c];;
gap> v := [1,0,3]*One(F);;
gap> StabilizerAbelianMatGroup( mats, v );
[ [ [ Z(7)^0, 0*Z(7), 0*Z(7) ], [ 0*Z(7), Z(7)^3, 0*Z(7) ], 
      [ 0*Z(7), 0*Z(7), Z(7)^0 ] ] ]
gap> w := v * a*c^3;;
gap> OrbitAbelianMatGroup( mats, v,w );
rec( 
  stab := 
    [ 
      [ [ Z(7)^0, 0*Z(7), 0*Z(7) ], [ 0*Z(7), Z(7)^3, 0*Z(7) ], 
          [ 0*Z(7), 0*Z(7), Z(7)^0 ] ] ], 
  y := [ [ Z(7)^3, 0*Z(7), 0*Z(7) ], [ 0*Z(7), Z(7)^0, 0*Z(7) ], 
      [ 0*Z(7), 0*Z(7), Z(7)^0 ] ] )
gap> STOP_TEST( "testall.tst", 10000 );
