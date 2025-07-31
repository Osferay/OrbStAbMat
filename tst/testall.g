#
# OrbStAbMat: Orbit stabilizer for abelian matrix groups
#
# This file runs package tests. It is also referenced in the package
# metadata in PackageInfo.g.
#
LoadPackage( "OrbStAbMat" );
dir := DirectoriesPackageLibrary( "OrbStAbMat", "tst" );
tst := [ "fos.tst", "orbstab.tst" ];

tst := List( tst, x -> Filename( dir, x ) );
TestDirectory( tst, rec(exitGAP := true));