// -*- C++ -*-
#include "Rivet/Analyses/MC_JetAnalysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "Rivet/Projections/HepMCHeavyIon.hh"

#include <fstream>
#include <map>

namespace Rivet {

/// @brief Analysis to study geometrical bias in heavy ion collisions
class GbiasInvestigatr : public MC_JetAnalysis {
public:

	GbiasInvestigatr()
	: MC_JetAnalysis("GbiasInvestigatr", 4, "Jets")
	{	}


	void init() {
		
		// project all final state particles
		FinalState fs(Cuts::abseta < 3.0);

		FastJets jetpro(fs, FastJets::ANTIKT, 0.4);
		declare(jetpro, "Jets");

		// variable initializations
		_eventNumber = 0; // keep track of the current event number
		_firstPrint = true; // in order to print the variables names on the file

		// init method from inherited class
		MC_JetAnalysis::init();


		// get file name
		string name;
		std::cout << std::endl << std::endl << "Out file name: ";
		std::cin >> name;
		if(name == "e\n") name = "eventdata.dat";
		std::cout << std::endl;

		// open file for writing
		_file.open(name, std::fstream::out);
		if(_file.is_open()) {
			_open = true; 
		} 
		else { 
			_open = false; 
			std::cout << "[WARNING] Could not open file." << std::endl;
		}

		return;
	}


	void analyze(const Event& event) {

		// read H line for collision variables
		const HepMC::HeavyIon *hion = event.genEvent()->heavy_ion();
		if(hion) {
			// JWL: polar angle for the jet creation relative to collision center
			_toWrite["Polar"] = hion->event_plane_angle();
			// JWL: distance from jet creation point to collision center
			_toWrite["JProdR"] = hion->eccentricity();
		}

		// find the two most energetic jets
		const Jets &jetlist = apply<FastJets>(event, "Jets").jetsByPt(Cuts::abseta < 2.0 && Cuts::pT > 20*GeV);

		if(jetlist.size() > 1) {
			
			Jet j1 = jetlist.front(); // first jet is highest energy jet
			if(j1.momentum().pT() < 80*GeV) vetoEvent; // reject events with jets whose pT < 80 GeV
			Jet j2; // second jet is the one with the angle vs j1 closest to pi

			// method for highest pT jet within PI/8 of 180 degrees
			bool jetFound = false;
			for(const Jet& j : jetlist) { 
				const double ang = fabs(j1.momentum().phi() - j.momentum().phi());
				if(fabs(ang-PI) < PI/8.) { // transversal angle
					j2 = j;
					// since the jets are ordered in pT, then the first jet that meets the criteria
					// of being around PI/8 from 180 degrees from the highest energy jet is the
					// one we want
					jetFound = true;
					break;
				}
			}
			if(!jetFound) vetoEvent;

			_toWrite["Jet1_pT"] = j1.momentum().pT();
			_toWrite["Jet2_pT"] = j2.momentum().pT();
			_toWrite["JetAngle"] = fabs(j1.momentum().phi() - j2.momentum().phi())*180/PI;
			_toWrite["Jet1_Ang"] = j1.momentum().phi();

			// calculate Aj
			double E1 = j1.momentum().pT();
			double E2 = j2.momentum().pT();
			double Aj = (E1 - E2) / (E1 + E2);
			_toWrite["Aj"] = Aj;
		}
		else {
			// no jets were found
			vetoEvent;
		}
		
		// event weight
		_toWrite["Weight"] = event.weights()[0];


		// write in file
		if(_firstPrint) { 
			// first print to write variable names in the table
			// jetsFound must be true otherwise not all Map variables were created
			print(true);
			_firstPrint = false;
		}
		++_eventNumber;
		print();

		// analyse method from inherited class
		MC_JetAnalysis::analyze(event);

		return;
	}


	void finalize() {

		// finalize method from inherited class
		MC_JetAnalysis::finalize();

		// close file
		if(_open) {
			_file.close();
			std::cout << "Everything written to file." << std::endl;
		}

		return;
	}


	void print(bool var = false) {
		
		if(!_file.is_open()) return;

		// print variable names
		if(var) {
			for(auto const& it : _toWrite) {
				_file << it.first << "	";
			}
			_file << '\n';
			return;
		}

		/*/ actually print stuff to terminal
		std::cout << "Event " << _eventNumber << std::endl;
		std::cout << "Polar: " << _toWrite["Polar"] << "    JProdR: " << _toWrite["JProdR"] << std::endl;
		std::cout << "Jet1 pT: " << _toWrite["Jet1_pT"] << "    Jet2 pT: " << _toWrite["Jet2_pT"];
		std::cout << "    Angle:" << _toWrite["JetAngle"] << std::endl;
		std::cout << "Weight: " << _toWrite["Weight"] << std::endl;
		std::cout << std::endl;//*/

		// actually print stuff to file
		for(auto const& it : _toWrite) {
			_file << it.second << "	";
		}
		_file << '\n';

		
		return;
	}


private:

	// class variables
	std::fstream _file;
	bool _open;

	int _eventNumber;
	bool _firstPrint;

	// variables to print from each event
	std::map<std::string,double> _toWrite;

	// other
	

};

// The hook for the plugin system
DECLARE_RIVET_PLUGIN(GbiasInvestigatr);

}