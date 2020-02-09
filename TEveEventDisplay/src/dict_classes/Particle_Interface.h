/*#ifndef Particle_Interface_h
#define Particle_Interface_h
#include <TObject.h>
#include <TSystem.h>
// ... libRIO
#include <TFile.h>
#include<string>

namespace mu2e{
	class Particle_Interface{
		public:
			#ifndef __CINT__
			explicit Particle_Interface();		
			virtual ~Particle_Interface();

			enum type{
				e_minus = 11,
				e_plus	= -11,
				mu_minus = 13,
				mu_plus = -13,
				pi_minus = -211,
				pi_plus = 211,
				proton = -2212,
				antiproton = 2212
			};

			enum linecolor{kRed, kGreen, kBlue, kYellow }; //e, mu, pi, p
			enum linestyle{kSolid, kDashed}; //+ or -

			double name() const { return _name;}
			double energy() const { return _energy;}
			double charge() const { return _charge;}
			double momentum() const { return _momentum;}
			double velocity() const { return _velocity;}


		private:
			double _energy;
			double _charge;
			double _momentum;
			double _velocity;
			double _name;

		#endif
		ClassDef(Particle_Interface, 0);
	};
}
#endif*/

