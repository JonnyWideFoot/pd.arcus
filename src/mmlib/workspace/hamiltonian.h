#ifndef __HAMILTONIAN_H
#define __HAMILTONIAN_H







//-------------------------------------------------
//
/// \brief  This structure holds all the energies and their individual components
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author Mike Tyka & Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class PD_API Hamiltonian
{
public:
	double etot, epot, ekin;
	double epot_vdw, epot_elec, epot_bond, epot_angle, epot_torsion, epot_surf, 
		epot_pol, epot_pol_cross, epot_pol_self, ekin_linear, ekin_angular;

	double InternalVirial;
	// Scratch space
	double data1; // arbitary double data that are saved into the trajectory
	double data2;
	double data3;
	double data4;
	double data5;
	double data6;
	double data7;
	double data8;
	double data9;
	double data10;

	// statistics
	double dRMS;
	double cRMS;

	// new
	double epot_vdw_att;
	double epot_vdw_rep;

	Hamiltonian()
	{
		zero();
	}

	void zero()
	{
		etot = 0.0;
		epot = 0.0;
		epot_vdw = 0.0;
		epot_elec = 0.0;
		epot_bond = 0.0;
		epot_angle = 0.0;
		epot_torsion = 0.0;
		epot_surf = 0.0;
		epot_pol = 0.0;
		epot_pol_cross = 0.0;
		epot_pol_self = 0.0;
		ekin = 0.0;
		ekin_linear = 0.0;
		ekin_angular = 0.0;
		epot_vdw_att = 0.0;
		epot_vdw_rep = 0.0;
		InternalVirial = 0.0;

		data1 = 0.0;
		data2 = 0.0;
		data3 = 0.0;
		data4 = 0.0;
		data5 = 0.0;
		data6 = 0.0;
		data7 = 0.0;
		data8 = 0.0;
		data9 = 0.0;
		data10 = 0.0;

		dRMS = 0.0;
		cRMS = 0.0;
	}
};

#endif

