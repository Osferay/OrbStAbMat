#
# OrbStAbMat: Orbit stabilizer for abelian matrix groups
# Reading global variables

if not IsBound(DEBUG@) then DEBUG@ := true; fi;

# Reading the implementation part of the package.
#
ReadPackage( "OrbStAbMat", "gap/general.gi");
ReadPackage( "OrbStAbMat", "gap/basis.gi");
ReadPackage( "OrbStAbMat", "gap/algebra.gi");
ReadPackage( "OrbStAbMat", "gap/series.gi");
ReadPackage( "OrbStAbMat", "gap/fos.gi");
ReadPackage( "OrbStAbMat", "gap/orbstab.gi");
