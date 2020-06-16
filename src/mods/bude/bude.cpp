#include "global.h"

#include "forcefields/bude.h"

namespace Physics
{
	BudeCustomProperties::BudeCustomProperties()
	{
		// Initialise variables - use these values to check that
		// the properties have been set properly in setup().

		ElectrostaticType = CHAR_MAX;
		Hardness = DBL_MAX;
		HydrophobicPotential = DBL_MAX;
		DistNpNp = DBL_MAX;
		DistNpP = DBL_MAX;
		RadiusScaling = DBL_MAX;
	}

	BudeForcefield::BudeForcefield(WorkSpace &newwspace)
		: ForcefieldBase(newwspace)
	{
		// Initialise member variables

		m_ReceptorIndex = -1;
		m_LigandIndex = -1;

		m_Cutoff = -1.0;
		m_Cutoff_elec_formal = -1.0; 
		m_Cutoff_elec_partial = -1.0;

		// Chose (PhysicsConst::econv)/22.5 because it's more easy to
		// see where the number came from than just putting 14.7485.
		// Alternatively could say 14.7485 (used to same number of
		// sig. figs as PhysicsConst::econv).
		m_Dielectric = (PhysicsConst::econv)/22.5;

		resetLocalEnergies();		
	}

	int BudeForcefield::setup()
	{
		// Proxies
		WorkSpace& wspace = getWSpace();
		size_t natom = wspace.natom();   ///< get the number of atoms in the workspace once in advance to save calculating it each loop

		if(OutputLevel) printf("\n\nBudeForcefield Setup: .....");

		resetLocalEnergies();           ///< Set all of the local energy stores to zero.
		m_BudeCustomProperties.clear(); ///< Clean this now in case it's already full.

		for( size_t i = 0; i < natom; i++ )
		{
			BudeCustomProperties _customProperties; // local store

			_customProperties.DistNpNp = wspace.atom[i].getCustomProperty_Double("DISTNPNP");
			_customProperties.DistNpP = wspace.atom[i].getCustomProperty_Double("DISTNPP");

			// Get the Electrostatic Type once and assign it to a local
			// variable rather than doing it for each if statement.
			std::string _elecType = wspace.atom[i].getCustomProperty_String("ELECTYPE");

			if( 0 == _elecType.compare("N") )
			{
				_customProperties.ElectrostaticType = m_elecType_N; // If it's N, store a 0.
			}
			else if( 0 == _elecType.compare("F") )
			{
				_customProperties.ElectrostaticType = m_elecType_F; // If it's F, store a 1.
			}
			else if( 0 == _elecType.compare("D") )
			{
				_customProperties.ElectrostaticType = m_elecType_D; // If it's D, store a 2.
			}
			else if( 0 == _elecType.compare("E") )
			{
				_customProperties.ElectrostaticType = m_elecType_E; // If it's E, store a 3.
			}
			else
			{
				throw ProcedureException("An atom's Electrostatic Type has been set incorrectly. Only N, F, D and E can currently exist");
			}

			_customProperties.Hardness = wspace.atom[i].getCustomProperty_Double("HARDNESS");
			_customProperties.HydrophobicPotential = wspace.atom[i].getCustomProperty_Double("HYDROPHOB");
			_customProperties.RadiusScaling = wspace.atom[i].getCustomProperty_Double("RADSCAL");

			m_BudeCustomProperties.push_back( _customProperties ); // Now each atom in the wspace will have all of its custom properties stored in m_BudeCustomProperties.

			/// Just make sure that all the custom properties for this atom have been added
			/// properly to the std::vector<BudeCustomProperties> m_BudeCustomProperties.
			if( (m_BudeCustomProperties[i].DistNpNp == DBL_MAX )
				|| (m_BudeCustomProperties[i].DistNpP == DBL_MAX )
				|| (m_BudeCustomProperties[i].ElectrostaticType == CHAR_MAX )
				|| (m_BudeCustomProperties[i].Hardness == DBL_MAX )
				|| (m_BudeCustomProperties[i].HydrophobicPotential == DBL_MAX )
				|| (m_BudeCustomProperties[i].RadiusScaling == DBL_MAX ) )
			{
				THROW(CodeException,"An error has occured during the entry of Custom Properties into an array. A property is still at its initialised value.");
			}

			/*printf("\nAtom: %4d,%3s, ElecType: %3d Radius: %7.4f Hydrophob: %7.4f Charge: %7.4f",
				i,
				wspace.atom[i].pdbname.c_str(),
				m_BudeCustomProperties[i].ElectrostaticType,
				wspace.atom[i].radius,
				m_BudeCustomProperties[i].HydrophobicPotential,
				wspace.atom[i].charge);*/

			/*printf("\nAtom: %4d,%3s, DistNpNp: %6.4f DistNpP: %6.4f ElecType: %3d Hardness: %6.4f Hydrophob: %7.4f RadScal: %6.4f",
				i,
				wspace.atom[i].pdbname.c_str(),
				m_BudeCustomProperties[i].DistNpNp,
				m_BudeCustomProperties[i].DistNpP,
				m_BudeCustomProperties[i].ElectrostaticType,
				m_BudeCustomProperties[i].Hardness,
				m_BudeCustomProperties[i].HydrophobicPotential,
				m_BudeCustomProperties[i].RadiusScaling);*/
		}

		needsetup = false;
		if(OutputLevel) printf("done\n\n");
		return 0;
	}

	void BudeForcefield::calcEnergies()
	{
		// Proxies
		WorkSpace& wspace = getWSpace();
		SnapShotAtom* pos = wspace.cur.atom;
		size_t natom = wspace.natom();       ///< get the number of atoms in the workspace once in advance to save calculating it each loop

		//if( needsetup )
		//{
		//	setup();
		//}

		resetLocalEnergies();    ///< Make sure all of the local energy stores are set to zero at the beginning.
		validateParams(wspace);  ///< Check the user has set all the necessary parameters.

		// Local variables
		const double _sqrCutoff = Maths::sqr( m_Cutoff ); // Square of m_Cutoff for comparison with _sqrDistij.
		
		/// Gives an answer of 22.5 for the conversion factor to fudge ~2.5kcal/mol H-bond.
		/// Calculated outside the loops as division is slow and it doesn't need to know
		/// what the current atoms are.
		double _kCalPerMolConversion = PhysicsConst::econv / m_Dielectric;
		
		size_t _AtomicMassS = 16; ///< Atomic mass of Sulphur.

		double _eSteric = 0.0;    ///< local store of steric energy. Initialised to 0.
		double _eElec = 0.0;      ///< local store of electrostatic energy. Initialised to 0.
		double _eDesolv = 0.0;    ///< local store of desolvation energy. Initialised to 0.

		const int _LigandStart = wspace.mol[m_LigandIndex].ifirst; ///< local store of the INDEX of the first ligand atom.
		const int _LigandEnd = wspace.mol[m_LigandIndex].ilast;    ///< local store of the INDEX of the last ligand atom.

		//printf("Index of Ligand first atom = %4d. Index of Ligand last atom = %4d\n",_LigandStart,_LigandEnd);

		// Loop over all atoms on the outside (i) and only ligand atoms on the inside (j)
		// thereby simply avoiding calculation of receptor-receptor interactions at all.
		for( size_t i = 0; i < natom; i++ )
		{
			for( int j = _LigandStart; j <= _LigandEnd; j++ )
			{
				
				double _sqrDistij = pos[i].p.sqrdist(pos[j].p); ///< square of the distance between atom i and j.
				
				// Don't calculate ligand internal energy (ie. IF both atoms are ligand atoms) if both atoms:
				//     - belong to the same residue,
				//     - are backbone atoms,
				//     - are sulphurs, atomic number 16 (to avoid problems with disulphide bonds).
				if( ( wspace.atom[i].imol == m_LigandIndex ) // Is i a ligand atom? We already know that j is only looping over ligand atoms.
					&& ( ( wspace.atom[i].ir == wspace.atom[j].ir )                                         // Are both atoms from the same residue?
					|| ( ( wspace.atom[i].isBackbone() ) && ( wspace.atom[j].isBackbone() ) )               // Are they both backbone atoms?
					|| ( ( wspace.atom[i].Z == _AtomicMassS ) && ( wspace.atom[j].Z == _AtomicMassS ) ) ) ) // Are they both Sulphurs?
				{
					continue;
				}
				
				// Don't do energy calculations if atoms are further apart than the cutoff.
				// (Use sqr distance and sqr cutoff because to calculate sqrt of distance is costly.)
				else if( _sqrDistij < _sqrCutoff )
				{
					// calculate distij and radij now to avoid repetition in all the energy calculations
					double distij = std::sqrt( _sqrDistij );                          ///< distance between atom i and atom j
					double radij = ( wspace.atom[i].radius + wspace.atom[j].radius ); ///< sum of the radii of atom i and atom j
					double _distijMinusRadij = distij - radij;                        ///< precalculated for speed.
					

					// --------------- Steric ---------------

					// Local stores
					_eSteric = 0.0; ///< Set local store of steric energy to 0 for each loop.
					
					// If  distij !< radij, don't do Steric calculation or store any energies.
					if( distij < radij )
					{
						double _HiPlusHjOver2 = ( m_BudeCustomProperties[i].Hardness + m_BudeCustomProperties[j].Hardness ) * 0.5; ///< pre-calculate this so only have to do it once.

						_eSteric = _HiPlusHjOver2 - ( _HiPlusHjOver2 * ( distij / radij ) );
					
						// --- store energies ---
						if( wspace.atom[i].imol == m_LigandIndex ) // If i and j are both ligand atoms.
						{
							m_EpotLigand_vdw += _eSteric; ///< add the Steric energy for the current atom pair to the total Ligand Internal Steric energy.
							m_Epot_vdw += _eSteric; ///< add the Steric energy for the current atom pair to the total Steric energy.
							//wspace.atom[j].epot += _eSteric;
						}
						else
						{
							m_EpotComplex_vdw += _eSteric; ///< add the Steric energy for the current atom pair to the total Complex (Ligand-Receptor) Steric energy.
							m_Epot_vdw += _eSteric; ///< add the Steric energy for the current atom pair to the total Steric energy.
							//wspace.atom[j].epot += _eSteric;
						}
					}

					// ------------ End of Steric ------------


					// ------------ Electrostatic ------------

					//Local stores
					_eElec = 0.0; ///< Set local store of electrostatic energy to 0 for each loop.
					char _ElecType_i = m_BudeCustomProperties[i].ElectrostaticType; ///< local store for atom i.
					char _ElecType_j = m_BudeCustomProperties[j].ElectrostaticType; ///< local store for atom j.

					// If both atoms have electrostatic type N then don't do
					// electrostatic calculation or store any energies.
					if( !( ( _ElecType_i == m_elecType_N ) || ( _ElecType_j == m_elecType_N ) ) )
					{
						// Precalculated values
						double _QiQjTimesConstant = (( wspace.atom[i].charge ) * ( wspace.atom[j].charge )) * _kCalPerMolConversion;


						// If distij < radij it doesn't matter whether an atom is F, D or E.
						if( distij < radij )
						{
							_eElec = _QiQjTimesConstant;
						}
						// If both atoms have electrostatic type F and distij < m_Cutoff_elec_formal
						else if( ( _ElecType_i == m_elecType_F ) && ( _ElecType_j == m_elecType_F ) && ( distij < m_Cutoff_elec_formal ) )
						{
							_eElec = ( _QiQjTimesConstant ) - ( ( _distijMinusRadij * _QiQjTimesConstant ) / ( m_Cutoff_elec_formal - radij ) );
						}
						else if( distij < m_Cutoff_elec_partial )
						{
							_eElec = ( _QiQjTimesConstant ) - ( ( _distijMinusRadij * _QiQjTimesConstant ) / ( m_Cutoff_elec_partial - radij ) );	
						}

						// If either atom is an E and the energy is positive then make
						// it negative because E can be H-bond donor or acceptor.
						if( (( _ElecType_i == m_elecType_E ) || ( _ElecType_j == m_elecType_E )) && ( _eElec > 0.0 ) )
						{
							_eElec = - _eElec;
						}

						// --- store energies ---
						if( wspace.atom[i].imol == m_LigandIndex ) // If i and j are both ligand atoms.
						{
							m_EpotLigand_elec += _eElec; ///< add the Electrostatic energy for the current atom pair to the total Ligand Internal Electrostatic energy.
							m_Epot_elec += _eElec; ///< add the Electrostatic energy for the current atom pair to the total Electrostatic energy.
							//wspace.atom[j].epot += _eElec;
							//printf("i: %4d j: %4d _eElec: %8.4f\n",i,j,_eElec);
						}
						else
						{
							m_EpotComplex_elec += _eElec; ///< add the Electrostatic energy for the current atom pair to the total Complex (Ligand-Receptor) Electrostatic energy.
							m_Epot_elec += _eElec; ///< add the Electrostatic energy for the current atom pair to the total Electrostatic energy.
							//wspace.atom[j].epot += _eElec;
							//printf("i: %4d j: %4d _eElec: %8.4f\n",i,j,_eElec);
						}
					}

					// --------- End of Electrostatic ---------


					// ------------- Desolvation -------------

					//Local stores
					_eDesolv = 0.0; ///< Set local store of desolvation energy to 0 for each loop.
					double _Ki = m_BudeCustomProperties[i].HydrophobicPotential;
					double _Kj = m_BudeCustomProperties[j].HydrophobicPotential;
					double _averageDistNpP = ( m_BudeCustomProperties[i].DistNpP + m_BudeCustomProperties[j].DistNpP ) * 0.5;
					double _averageDistNpNp = ( m_BudeCustomProperties[i].DistNpNp + m_BudeCustomProperties[j].DistNpNp ) * 0.5;

					// Don't do any desolvation calculations or store any energies if:
					//     - Ki or Kj are == 0
					//     - Ki and Kj are both > 0
					if( (_Ki != 0.0) && (_Kj != 0.0) && !( ( _Ki > 0.0 ) && ( _Kj > 0.0 ) ) )
					{
						/// If both atoms are hydrophobic ie. Ki < 0 and Kj < 0.
						if( ( _Ki < 0.0 ) && ( _Kj < 0.0 ) )
						{
							double _KiPlusKjOver2 = ( _Ki + _Kj ) * 0.5;

							if( distij < radij )
							{
								_eDesolv = _KiPlusKjOver2;
							}
							else if( distij < ( radij + _averageDistNpNp ) )
							{
								_eDesolv = _KiPlusKjOver2 - ( ( _distijMinusRadij * _KiPlusKjOver2 ) / _averageDistNpNp );
							}
						}
						/// Otherwise, if one atom is hydrophobic and one is hydrophilic.
						else
						{
							double _ModKi = std::abs(_Ki); ///< Modulus of Ki
							double _ModKj = std::abs(_Kj);///< Modulus of Kj
							double _ModKiPlusModKjOver2 = ( _ModKi + _ModKj ) * 0.5;

							if( distij < radij )
							{
								_eDesolv = _ModKiPlusModKjOver2;
								//printf("i: %4d j: %4d _eDesolv: %8.4f\n",i,j,_eDesolv);
							}
							else if( distij < ( radij + _averageDistNpP ) )
							{
								_eDesolv = _ModKiPlusModKjOver2 - ( ( _distijMinusRadij * _ModKiPlusModKjOver2 ) / _averageDistNpP );
							}
						}
					}
					
					// --- store energies ---
					if( wspace.atom[i].imol == m_LigandIndex ) // If i and j are both ligand atoms.
					{
						m_EpotLigand_desolv += _eDesolv; ///< add the Desolvation energy for the current atom pair to the total Ligand Internal Desolvation energy.
						m_Epot_desolv += _eDesolv; ///< add the Desolvation energy for the current atom pair to the total Desolvation energy.
						//wspace.atom[j].epot += _eDesolv;
					}
					else
					{
						m_EpotComplex_desolv += _eDesolv; ///< add the Desolvation energy for the current atom pair to the total Complex (Ligand-Receptor) Desolvation energy.
						m_Epot_desolv += _eDesolv; ///< add the Desolvation energy for the current atom pair to the total Desolvation energy.
						//wspace.atom[j].epot += _eDesolv;
					}
					
					/*

					// The following large commented out block is a version of the desolvation calculation exactly
					// as it is found in the old dockit. I've left it in here for the time being in case I want to
					// compare it for speed with the desolvation calculation above when optimising the code.

					double _DistNpP_i = m_BudeCustomProperties[i].DistNpP;
					double _DistNpP_j = m_BudeCustomProperties[j].DistNpP;
					double _DistNpNp_i = m_BudeCustomProperties[i].DistNpNp;
					double _DistNpNp_j = m_BudeCustomProperties[j].DistNpNp;

					// Pre-calculated to increase readability and speed.
					double _KiPlusKjOver2 = ( _Ki + _Kj ) * 0.5;
					double _KiMinusKjOver2 = ( _Ki - _Kj ) * 0.5;
					double _KjMinusKiOver2 = ( _Kj - _Ki ) * 0.5;
					double _averageDistNpP = ( m_BudeCustomProperties[i].DistNpP + m_BudeCustomProperties[j].DistNpP ) * 0.5;
					double _averageDistNpNp = ( m_BudeCustomProperties[i].DistNpNp + m_BudeCustomProperties[j].DistNpNp ) * 0.5;

					if( ( _Ki > 0.0 ) && ( _Kj < 0.0 ) )
					{
						if( distij < radij )
						{
							_eDesolv = _KiMinusKjOver2;
							printf("i: %4d j: %4d _eDesolv: %8.4f\n",i,j,_eDesolv);
						}
						else if( distij < ( radij + _averageDistNpP ) )
						{
							_eDesolv = _KiMinusKjOver2 - ( ( _distijMinusRadij * _KiMinusKjOver2 ) / _averageDistNpP );
						}
					}
					else if( ( _Ki < 0.0 ) && ( _Kj > 0.0 ) )
					{
						if( distij < radij )
						{
							_eDesolv = _KjMinusKiOver2;
							printf("i: %4d j: %4d _eDesolv: %8.4f\n",i,j,_eDesolv);
						}
						else if( distij < ( radij + _averageDistNpP ) )
						{
							_eDesolv = _KjMinusKiOver2 - ( ( _distijMinusRadij * _KjMinusKiOver2 ) / _averageDistNpP );
						}
					}

					else if( ( _Ki < 0.0 ) && ( _Kj < 0.0 ) )
					{
						if( distij < radij )
						{
							_eDesolv = _KiPlusKjOver2;
						}
						else if( distij < ( radij + _averageDistNpNp ) )
						{
							_eDesolv = _KiPlusKjOver2 - ( ( _distijMinusRadij * _KiPlusKjOver2 ) / _averageDistNpNp );
						}
					}

					// --- store energies ---
					if( wspace.atom[i].imol == m_LigandIndex ) // If i and j are both ligand atoms.
					{
						m_EpotLigand_desolv += _eDesolv; ///< add the Desolvation energy for the current atom pair to the total Ligand Internal Desolvation energy.
						m_Epot_desolv += _eDesolv; ///< add the Desolvation energy for the current atom pair to the total Desolvation energy.
					}
					else
					{
						m_EpotComplex_desolv += _eDesolv; ///< add the Desolvation energy for the current atom pair to the total Complex (Ligand-Receptor) Desolvation energy.
						m_Epot_desolv += _eDesolv; ///< add the Desolvation energy for the current atom pair to the total Desolvation energy.
					}

					*/

				// ---------- End of Desolvation ----------

				}
			}
		}

		// Will I need to multiply the calculated energies by (PhysicsConst::kcal2J / PhysicsConst::Na)
		// to get them into SI units since internally the energies are in J not kcal/mol?
		// THE DECISION IS 'NO' BECAUSE:
		//      - Electrostatic calculation already has a conversion factor included.
		//      - The atom hardness and hydrophobic (desolvation) potential used in
		//        Steric and desolvation calculations have been parameterised so that
		//        the energies are equivalent to kcal/mol.

		// Calculate the totals - sum of vdw, electrostatic and desolvation energies.
		m_EpotLigand_total = m_EpotLigand_vdw + m_EpotLigand_elec + m_EpotLigand_desolv;
		m_EpotComplex_total = m_EpotComplex_vdw + m_EpotComplex_elec + m_EpotComplex_desolv;
		m_Epot_total = m_Epot_vdw + m_Epot_elec + m_Epot_desolv;

		// Make sure that (m_Epot_vdw + m_Epot_elec + m_Epot_desolv) == m_EpotLigand_total + m_EpotComplex_total.
		D_ASSERT(Maths::SigFigEquality(m_Epot_total, ( m_EpotLigand_total + m_EpotComplex_total ), 9),
				CodeException, 
				"There is an error in the code, the condition 'Total energy == Ligand Internal Energy + Complex Energy' must be true.");

		// Store the energies from the local stores in the wspace Hamiltonian
		// individual potential energies vdw, elec and pol (desolv).
		wspace.ene.epot_vdw = m_Epot_vdw;
		wspace.ene.epot_elec = m_Epot_elec;
		wspace.ene.epot_pol = m_Epot_desolv;

		// and the total potential energy.
		wspace.ene.epot = m_Epot_total;

		// ----- DEBUGGING -----
		//for( int j = _LigandStart; j <= _LigandEnd; j++ )
		//{
		//	printf("At: %4d, %3s, res %3s Electype: %3d Charge: %8.4f atom[j].epot: %8.4f\n",j,wspace.atom[j].pdbname.c_str(),wspace.atom[j].parentl3name.c_str(),m_BudeCustomProperties[j].ElectrostaticType,wspace.atom[j].charge,wspace.atom[j].epot);
		//}
		// ---------------------
	}

	void BudeForcefield::calcEnergiesVerbose(Verbosity::Type level)
	{
		calcEnergies();
	}

	void BudeForcefield::resetLocalEnergies()
	{
		m_EpotLigand_vdw = 0.0;    
		m_EpotLigand_elec = 0.0;   
		m_EpotLigand_desolv = 0.0; 
		m_EpotLigand_total = 0.0;  

		m_EpotComplex_vdw = 0.0;   
		m_EpotComplex_elec = 0.0;  
		m_EpotComplex_desolv = 0.0;
		m_EpotComplex_total = 0.0; 

		m_Epot_vdw = 0.0;   
		m_Epot_elec = 0.0;  
		m_Epot_desolv = 0.0;
		m_Epot_total = 0.0;
	}
	
	void BudeForcefield::validateParams(const WorkSpace& _wspace) const
	{
		if( m_Cutoff >= 0 )
		{
			if( m_Cutoff_elec_formal >= 0 )
			{
				if( m_Cutoff_elec_partial >= 0 )
				{
					if( (m_LigandIndex >= 0) && (m_ReceptorIndex >= 0)
						&& (m_LigandIndex != m_ReceptorIndex)
						&& (m_LigandIndex < _wspace.nmol()) && (m_ReceptorIndex < _wspace.nmol()) )
					{
						return;
					}
					else
					{
						throw ProcedureException("Ligand Index and Receptor Index have been set incorrectly");
					}
				}
				else
				{
					throw ProcedureException("Cutoff for electrostatic interactions between partial charges (m_Cutoff_elec_partial) has been set incorrectly");
				}
			}
			else
			{
				throw ProcedureException("Cutoff for electrostatic interactions between formal charges (m_Cutoff_elec_formal) has been set incorrectly");
			}
		}
		else
		{
			throw ProcedureException("Global Cutoff (m_Cutoff) has been set incorrectly");
		}
	}

	// ---------- Set and Get functions ---------

	void BudeForcefield::setReceptorIndex(const int& _ReceptorIndex)
	{
		m_ReceptorIndex = _ReceptorIndex;
	}

	const int& BudeForcefield::getReceptorIndex() const
	{
		return m_ReceptorIndex;
	}

	void BudeForcefield::setLigandIndex(const int& _LigandIndex)
	{
		m_LigandIndex = _LigandIndex;
	}

	const int& BudeForcefield::getLigandIndex() const
	{
		return m_LigandIndex;
	}

	void BudeForcefield::setCutoff(const double& _cutoff)
	{
		m_Cutoff = _cutoff;
	}

	const double& BudeForcefield::getCutoff() const
	{
		return m_Cutoff;
	}

	void BudeForcefield::setCutoff_elec_formal(const double& _cutoff_elec_formal)
	{
		m_Cutoff_elec_formal = _cutoff_elec_formal;
	}

	const double& BudeForcefield::getCutoff_elec_formal() const
	{
		return m_Cutoff_elec_formal;
	}

	void BudeForcefield::setCutoff_elec_partial(const double& _cutoff_elec_partial)
	{
		m_Cutoff_elec_partial = _cutoff_elec_partial;
	}

	const double& BudeForcefield::getCutoff_elec_partial() const
	{
		return m_Cutoff_elec_partial;
	}

	// ---------- Info functions ---------
	
	void BudeForcefield::infoLine() const
	{
		printf( "%17s       %8.3f       %8.3f       %8.3f       %8.3f\n",
			"Ligand Internal: ",
			m_EpotLigand_vdw,
			m_EpotLigand_elec,
			m_EpotLigand_desolv,
			m_EpotLigand_total );

		printf( "%17s       %8.3f       %8.3f       %8.3f       %8.3f\n",
			"Complex: ",
			m_EpotComplex_vdw,
			m_EpotComplex_elec,
			m_EpotComplex_desolv,
			m_EpotComplex_total );

		printf( "%17s       %8.3f       %8.3f       %8.3f       %8.3f\n",
			"Total: ",
			m_Epot_vdw,
			m_Epot_elec,
			m_Epot_desolv,
			m_Epot_total );
	}

	void BudeForcefield::infoLineHeader() const
	{
		printf("\n%15s    %13s  %13s  %13s  %13s\n\n","Energy","Steric","Electrostatic","Desolvation","Total");
	}
}

