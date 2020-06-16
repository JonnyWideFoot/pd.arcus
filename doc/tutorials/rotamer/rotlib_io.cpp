// The code below details possible import options.
// The user only ever wants one of these !!

// **NOTE**
// Readlib imports PD format rotamer libraries
// convertLib imports other format rotamer libraries by taking an additional Conversion class

using namespace Library;

RotamerLibrary rot(ffps);

// 1) Shetty format
// Import from a shetty file-pair
rot.convertLib("lib/rotlib/shetty/scl-B30-occ1.0-rmsd1.0-prop20.0",RotLibConvert_Shetty());
rot.writeLib("lib/rotlib/shetty.rotamer"); // write the lib in PD format to this path

// 2) Legacy PD format
// This will clear the previous import and start again
// Import an old-format PD library
rot.convertLib("lib/legacy.rotamer",RotLibConvert_OldPDFormat());

// 2b) Same for Bude; it used Legacy PD format
rot.convertLib("lib/rotlib/bude.legacy.rotamer",RotLibConvert_OldPDFormat());
rot.writeLib("bude.cart.rotamer");

// 3) Torsional:Cartesian interconversion! (Cartesian application has benchmarked as faster. Torsional models can be imported,
// and are, by default, applied as cartesian coordinated).
rot.writeLib("output.cart.rotamer"); // it was imported as a cartesian lib, so this will write cartesian definitions.
rot.writeLib("output.torsion.rotamer",Library::Torsional); // Force output in torsional format
// automatically clears the current library and imports again.
rot.readLib("output.torsion.rotamer"); // This time import is of torsional definitions!
// We can now force writing cartesian definitions. They should be the same, bar floating-point rounding
// and anomalies due to bond length and angle deviations between the original definitions and the idealised forms used in torsional conversion.
rot.writeLib("output2.cart.rotamer",Library::Cartesian); // (should match "output.torsion.rotamer")

// 4) Dunbrack format
// Read the Dunbrack backbone independent library
rot.convertLib("lib/rotlib/dunbrack/bbind02.May.rot",RotLibConvert_Dunbrack_BBInd());

// Read the Dunbrack backbone independent library
RotLibConvert_Dunbrack_BBDep dunbrackConv();
rot.convertLib("lib/rotlib/dunbrack/bbdep02.May.rot",dunbrackConv); // By-chi mode (best, default)
dunbrackConv.switchImportMode = true; // Flag a change to by-rotid mode
rot.convertLib("lib/rotlib/dunbrack/bbdep02.May.rot",dunbrackConv); // By-rotid mode (alternative)