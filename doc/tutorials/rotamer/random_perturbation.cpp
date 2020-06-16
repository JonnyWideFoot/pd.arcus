// This should later be converted to a .py file and complementary info added

void randomRotamerApplication( bool _testSterics, bool importForeignLib )
{
	FFParamSet ffps;
	ffps.readLib("lib/amber03aa.ff");

	System sys( ffps );
	sys.add( NewProteinHelix(ffps,"*A-(CDEFGHIKLMNPQRSTVW)-Y*") );

	WorkSpace wspace( sys );
	wspace.info();
	wspace.printPDB("rotamer_perturb_random_start.pdb");
	wspace.addStdTra("rotamer_perturb_random_end");

	Library::RotamerLibrary rot(ffps);

	if(importForeignLib)
	{
		rot.convertLib("lib/rotlib/shetty/scl-B30-occ1.0-rmsd1.0-prop20.0",Library::RotLibConvert_Shetty());
		rot.writeLib("lib/rotlib/shetty.rotamer");
	}
	else
	{
		rot.readLib("lib/rotlib/shetty.rotamer");
	}

	StatClock timeMe;
	timeMe.Begin();

	RandomRotamerApplicator app( wspace, rot, ApplyCartesian );
	if( _testSterics )
		app.addFilter_DefaultSteric(); // Add the linear-time clash filter for all rotamer states
	app.test(500); // 500 random rotamer applications over the entire PickedResidueList

	timeMe.End();
	timeMe.ReportMilliSeconds();

	wspace.printPDB("rotamer_perturb_random_end.pdb");
}