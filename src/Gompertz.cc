#include <module/Module.h> // include JAGS module base class
#include <distributions/DGomp.h> // include Gompertz distribution class

namespace Gompertz { // start defining the module namespace
	// Module class
	class GOMPModule : public jags::Module {
		public :
		GOMPModule (); // constructor
		~GOMPModule(); // destructor
	};

	// Constructor function
	GOMPModule::GOMPModule() : Module ("Gompertz") {
		insert(new DGomp); // inherited function to load objects into JAGS
	}
	// Destructor function
	GOMPModule::~GOMPModule() {
		std::vector<jags::Distribution*> const &dvec = distributions();
		for (auto *dist : dvec) {
			delete dist;
		}
	}
} // end namespace definition

Gompertz::GOMPModule _Gompertz_module;
