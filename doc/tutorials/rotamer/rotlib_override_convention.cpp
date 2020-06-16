// Advanced!

Library::RotamerLibrary rot(ffps);
// Having a CYM around is a bit of an arse when it comes to making rotamer defs :-D
// A CYM is just a CYS with no terminal hydrogen!
rot.addIonisationAlias("CYS","CYM");
Library::RotConvention conv;
conv.setToDefaultConvention(); // load all normal amino acids
Library::ConventionDef cd;
cd.resName = "CYM";
cd.a1 = "CA";
cd.a2 = "N";
cd.a3 = "CB";
cd.chi.push_back(Library::ChiDef("N","CA","CB","SG"));
conv.addConvention( cd );
rot.overrideConventions( conv );

rot.readLib();