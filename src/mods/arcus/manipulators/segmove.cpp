#include "global.h"
#include "segmove.h"
#include "maths/fastrandom.h"
#include "library/rotamerlib.h"
#include "pickers/pickbase.h"


// Namespaces
using namespace Maths;

namespace Manipulator
{
	SegPerturbation::SegPerturbation( SegmentDef &_segdef ):
			SegmentDefUser(_segdef),
			MoveBase( _segdef.getWorkSpace() )
	{
		m_Rand = FastRandom::getInstance();
	}

	SegPerturbation::~SegPerturbation()
	{
		// nothing to clean up in the base class ...
	}

	// ------------------------------------
	// Begin: SegPerturbation_BBTorsional
	// ------------------------------------
	int SegPerturbation_BBTorsional::apply()
	{
		//int moveseverity = 0;
		for(size_t i = 0; i < segLength(); i++)
		{
			if( m_Rand->nextDouble() < m_Probability )
			{
				m_RotatePhi_AA[i].perturbRotation( m_Rand->nextNormal(m_MaxMagnitude) );
				m_RotatePsi_AA[i].perturbRotation( m_Rand->nextNormal(m_MaxMagnitude) );
			}
		}
		//return moveseverity; // can be '0' if no move has occured
		return 1;
	}

	// used to set the probabilities used by this class
	void SegPerturbation_BBTorsional::Set( double probability, double maxMagnitude )
	{
		m_Probability = probability;
		m_MaxMagnitude = maxMagnitude;
	}
	// ------------------------------------
	// End: SegPerturbation_BBTorsional
	// ------------------------------------

	// ------------------------------------
	// Begin: SegPerturbation_SCTorsional
	// ------------------------------------

	// used to set the probabilities used by this class
	void SegPerturbation_SCTorsional::Set( double probability, double maxMagnitude, double probability120 )
	{
		m_Probability120 = probability120;
		m_Probability = probability;
		m_MaxMagnitude = maxMagnitude;
	}

	int SegPerturbation_SCTorsional::apply()
	{
		RotBond& rotbond = wspace->rotbond();

		size_t j;
		for( size_t i = 0; i < segLength(); i++ )
		{
			RotBondCache& rotDef = m_RotBonds[i];
			for( j = 0; j < rotDef.size(); j++ )
			{
				if( m_Rand->nextDouble() < m_Probability120 )
				{
					// we want to do a 120 degree rotation, decide which direction...
					// angle to rotate (120 Degrees == 2.09439504 radians)
					if( m_Rand->nextBool() ) // use a bool, they're faster :-D
					{
						rotbond.rotate( rotDef.at(j), 2.09439504 );
					}
					else
					{
						rotbond.rotate( rotDef.at(j), -2.09439504 );
					}
				}
				else if( m_Rand->nextDouble() < m_Probability )
				{
					rotbond.rotate( rotDef.at(j), (double)m_Rand->nextNormal(m_MaxMagnitude) );
				}
			}
		}
		return 1; // 'Success': Perturbation magnitude is returned, we default to 1. We should return 0 if no perturbation occurs and 2 if a large perturbations occured, however there are not currently used.
	}

	// ------------------------------------
	// End: SegPerturbation_SCTorsional
	// ------------------------------------

	// ------------------------------------
	// Begin: SegPerturbation_SCRotamer
	// ------------------------------------

	SegPerturbation_SCRotamer::SegPerturbation_SCRotamer( 
		SegmentDef &_segdef,
		const Library::RotamerLibrary &_rotlib,
		double probability, 
		double forcedChangeCutoff,
		RotamerMode _mode ) 
		: SegmentDefUser( _segdef ),
		RotamerApplicatorBase( _segdef.getWorkSpace(), _rotlib, _mode ),
		rotlib(&_rotlib) 
	{
		m_Rand = FastRandom::getInstance();
		Set( probability, forcedChangeCutoff );
	}

	// public function to change the per-residue probability that a rotamer-change will occur
	void SegPerturbation_SCRotamer::Set( double probability, double forcedChangeCutoff )
	{
		m_Probability = probability;
		m_forcedchangecutoff = forcedChangeCutoff;
	}

	// Call this function to apply the best perturb rotamer states based on the internal 'm_Probability'
	// If a perturbation occurs then the environment is assessed to see if the application is applicable
	// for each rotamer state. The best one is chosen...
	int SegPerturbation_SCRotamer::apply()
	{
		int ir, irot;
		int nrots;
		int nvalidrots;
		double rotamersterics[100];
		int validrotamer[100];
		double stericlimit = 4.0;
		bool tryallrotamers;

		int endAt = segStartIndex() + segLength();
		for(ir = segStartIndex(); ir < endAt; ir++)
		{
			tryallrotamers = m_Rand->nextDouble() < m_Probability;

			nrots = getRotSterics( ir, &rotamersterics[0], tryallrotamers );

			if((!tryallrotamers) && (nrots <= 1))
				continue; // dont do anything if rotamer is ok and we dont force a change

			if(nrots <= 1)
				continue; // no rotamers to chose from

			// if(!changerotamer) printf("Forced change of rotamer ir=%d : steric score: %lf -> choices %d\n",
			// ir,rotamersterics[0],nrots);

			nvalidrots = 0;

			for(irot = 1; irot < nrots; irot++)
			{
				// printf("%lf \n",rotamersterics[irot]);
				if(rotamersterics[irot] < stericlimit)
				{
					validrotamer[nvalidrots] = irot;
					nvalidrots++;
				}
			}

			if(nvalidrots <= 0) // no rotamers to chose from, at least chose the one with the smallest steric score
			{
				double lowestSteric = DBL_MAX;
				int lowestStericIndex = 0;

				for(irot = 0; irot < nrots; irot++)
				{
					// printf("%lf \n",rotamersterics[irot]);
					if(rotamersterics[irot] < lowestSteric)
					{
						lowestSteric = rotamersterics[irot];
						lowestStericIndex = irot;
					}
				}

				if(lowestStericIndex == 0)
					continue; // dont do anything if the current conformation is already the lowest steric

				irot = lowestStericIndex; // otherwise change sidechain to the rotamer with the lowest steric value
			}
			else
			{
				// if several rotamers are fine choose randomly amongst the ok ones
				irot = validrotamer[rand() % nvalidrots];
			}

			// printf("Changing: %d \n",irot);
			applyRotamer(ir,irot - 1);
		}

		return 1;
	}

	// This function will do two related things:
	// a) it will check the steric situation of the current sidechain conformation or the
	// residue ir and suggest a change if above a certain threshold (scorecutoff) (using b) )
	// b) if tryallrotamers = true it will assess the entire rotamer distribution
	int SegPerturbation_SCRotamer::getRotSterics(int ir, double *rotsterics, bool tryallrotamers )
	{
		dvector *saveposition = NULL;
		int irot = 0;
		double sqrdistij;
		double radius = 3.4;
		double sqrradius = sqr(radius);
		double score;
		double scorei;
		double scoreij;

		const ParticleStore& atom = wspace->atom;
		const ResidueStore& res = wspace->res;
		SnapShotAtom *pos = wspace->cur.atom; // atom coordinate array

		// first save atom positions
		saveposition = new dvector[res[ir].ilast - res[ir].ifirst + 1];

		int iFirst =  res[ir].ifirst;
		for(int i = iFirst; i <= res[ir].ilast; i++)
		{
			saveposition[i-iFirst].setTo(pos[i].p);
		}

		// got through every possible rotamer, apply it, and calculate the
		// local steric impact
		// the first round through this seg accesses the current rotamer, its steric situation
		// information will be saved in rotsterics[0]
		for( size_t k = 0; k < nRot(ir)+1; k++ )
		{
			if( k != 0 ) applyRotamer(ir,k-1);

			// estimate the atomic overlaps;
			score = 0.0;

			for(int i = res[ir].ifirst; i <= res[ir].ilast; i++)
			{
				if(atom[i].isHydrogen())
					continue; // ignore hydrogens
				if(i == res[ir].iCA)
					continue; // skip backbone atoms
				if(i == res[ir].iC)
					continue; // skip backbone atoms
				if(i == res[ir].iN)
					continue; // skip backbone atoms
				if(i == res[ir].iH)
					continue; // skip backbone atoms
				if(i == res[ir].iO)
					continue; // skip backbone atoms

				scorei = 0;

				int nAtom = (size_t) atom.size();
				for(int j = 0; j < nAtom; j++)
				{
					// skip self interactions
					if((j >= res[ir].ifirst) && (j <= res[ir].ilast))
					{
						j = res[ir].ilast + 1;
						continue;
					}

					if(atom[j].isHydrogen())
						continue; // ignore hydrogens

					// now calculate the squaredistance between i & j
					sqrdistij = pos[i].p.sqrdist( pos[j].p );

					if(sqrdistij > sqrradius)
						continue;

					scoreij = sqr(1.0 - (sqrdistij / sqrradius));
					scorei += scoreij;
				}
				// printf("scorei: %lf \n",scorei);
				score += scorei;
			}

			//printf("Rot: %3d score = %lf \n",irot,score);
			//ene.epot = score * PhysicsConst::kcal2J / PhysicsConst::Na;
			//if(score < 3.5) outtra.append();

			rotsterics[irot] = score;
			irot++;

			if( (irot == 1) && (!tryallrotamers) && // if the current conformation (irot==1) is ok sterically
				(score < m_forcedchangecutoff)) // and caller has not explicitly asked for all rotameters, we quit
			{
				break;
			}
			// irot = 0 means current conformation, irot = 1 means the 0th rotamer,
			// hence the -1 in the line below
		}

		// restore positions
		for(int i = res[ir].ifirst; i <= res[ir].ilast; i++)
		{
			pos[i].p.setTo(saveposition[i - res[ir].ifirst]);
		}

		delete[] saveposition;

		return irot;
	}


	// ------------------------------------
	// End: SegPerturbation_SCRotamer
	// ------------------------------------

	// ------------------------------------
	// Begin: SegTweaker
	// ------------------------------------

	SegTweaker::SegTweaker()
	{
	}

	void SegTweaker::applyPerturbations()
	{
		for( size_t i = 0; i < size(); i++ )
		{
			element(i).apply();
		}
	}

	// ------------------------------------
	// End: SegTweaker
	// ------------------------------------
}

