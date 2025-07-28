#
# OrbStAbMat: Orbit stabilizer for abelian matrix groups
#
# This file contains package meta data. For additional information on
# the meaning and correct usage of these fields, please consult the
# manual of the "Example" package as well as the comments in its
# PackageInfo.g file.
#
SetPackageInfo( rec(

PackageName := "OrbStAbMat",
Subtitle := "Orbit stabilizer for abelian matrix groups",
Version := "0.1",
Date := "28/07/2025", # dd/mm/yyyy format
License := "GPL-2.0-or-later",

Persons := [
  rec(
    FirstNames := "Oscar",
    LastName := "Fernandez Ayala",
    WWWHome := "https://osferay.github.io/",
    Email := "oscar00ayala@gmail.com",
    IsAuthor := true,
    IsMaintainer := true,
    #PostalAddress := TODO,
    Place := "Braunschweig",
    Institution := "TU Braunschweig",
  ),
],

SourceRepository := rec(
    Type := "git",
    URL := "https://github.com/Osferay/OrbStAbMat",
),
IssueTrackerURL := Concatenation( ~.SourceRepository.URL, "/issues" ),
PackageWWWHome  := "https://Osferay.github.io/OrbStAbMat/",
PackageInfoURL  := Concatenation( ~.PackageWWWHome, "PackageInfo.g" ),
README_URL      := Concatenation( ~.PackageWWWHome, "README.md" ),
ArchiveURL      := Concatenation( ~.SourceRepository.URL,
                                 "/releases/download/v", ~.Version,
                                 "/", ~.PackageName, "-", ~.Version ),

ArchiveFormats := ".tar.gz",

AbstractHTML   :=  "",

PackageDoc := rec(
  BookName  := "OrbStAbMat",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0_mj.html",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "Orbit stabilizer for abelian matrix groups",
),

Dependencies := rec(
  GAP := ">= 4.13",
  NeededOtherPackages := [ ],
  SuggestedOtherPackages := [ ],
  ExternalConditions := [ ],
),

AvailabilityTest := ReturnTrue,

TestFile := "tst/testall.g",

#Keywords := [ "TODO" ],

));


