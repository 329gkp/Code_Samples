/*Code sample: software developed by the IceCube Neutrino Observatory is kept in private repos, so here I share a limited example of my work.
This program calculates the evolution of oscillation probabilities for particles ("neutrinos") that are born in the atmosphere and travel through Earth. 
An analysis I performed searched for signals of a phenomenon called decoherence that can affect neutrino propagation, which is the effect I implement here.
Note that the original code has been reduced for brevity.
The code below relies on nuSQuIDS and SQuIDS; nusquids is here https://github.com/arguelles/nuSQuIDS and squids is here https://github.com/jsalvado/SQuIDS/ .
Tools used to implement decoherence in nusquids can be found here https://github.com/ts4051/nuSQuIDS */ 

#include <vector>
#include <iostream>
#include <string>
#include <nuSQuIDS/nuSQuIDS.h>
#include <nuSQuIDS/marray.h>
#include <nuSQuIDS/nuSQuIDSDecoh.h>
#include <sstream>
#include <fstream>

using namespace nusquids;
int main(int argc, char* argv[]){

  // SPECIFY EARTH DENSITY MODEL AND INITIAL ATMOSPHERIC FLUX ("flux") LOCATION ("input_path") AND FINAL FLUX DESTINATION ("output_path")
  std::string input_path;
  std::string output_path;
  std::string flux;
  input_path    = argv[1];
  output_path   = argv[2];
  flux          = argv[3];

  const squids::Const units; // LOAD RELEVANT UNITS

  // SPECIFY DECOHERENCE STRENGTH "gamma_0", ENERGY POWER LAW VALUE "n_energy", AND DECOHERENCE MODEL "mode"
  double gamma_0 = atof(argv[4]);
  double n_energy = atof(argv[5]);
  std::string mode = argv[6];

  marray<double,1> diag({9},{0.,0.,0.,0.,0.,0.,0.,0.,0.}); // CREATE EMPTY DECOHERENCE OPERATOR 

  if(mode=="stateselection") { // FOR BREVITY, I'LL USE THIS TOY MODEL AS AN EXAMPLE
    diag[1] = gamma_0;
    } else{
        throw std::runtime_error("Did not specify available mode.");
    }

  if(output_path[output_path.length()-1]!='/'){ // CHECK OUTPUT DESTINATION FORMAT
    output_path = output_path +'/';}

  const unsigned int numneu = 3;

  // FORMAT: nus_atm( incident cosine zenith angle range, initial energies, number of flavors, neutrinos or antineutrinos, enable neutrino-matter interactions )
  nuSQUIDSAtm<nuSQUIDSDecoh> nus_atm(linspace(-1.,0.2,100),logspace(1.e2*units.GeV,1.e6*units.GeV,350),numneu,both,true);

  // INSERT DECOHERENCE WITH 1000 TEV PIVOT ENERGY FOR RUNTIME OPTIMIZATION
  for(nuSQUIDSDecoh& nsq : nus_atm.GetnuSQuIDS()){
  nsq.EnableDecoherence(true);
  nsq.Set_DecoherenceGammaMatrixDiagonal(diag);
  nsq.Set_DecoherenceGammaEnergyDependence(n_energy);
  nsq.Set_DecoherenceGammaEnergyScale(1000000.*units.GeV);
  }

  std::shared_ptr<EarthAtm> earth = std::make_shared<EarthAtm>(input_path +"/EARTH_MODEL_PREM.dat"); // LOAD EARTH MATTER DENSITY PROFILE
  nus_atm.Set_EarthModel(earth);

  double error = 1.0e-20; // SETUP INTEGRATION PRECISION
  nus_atm.Set_rel_error(error);
  nus_atm.Set_abs_error(error);

  // LOAD INITIAL ATMOSPHERIC NEUTRINO FLUX, CREATE NEUTRINO INITIAL STATE ARRAY
  marray<double,2> input_flux = quickread(input_path + "/" + flux);
  marray<double,4> inistate {nus_atm.GetNumCos(),nus_atm.GetNumE(),2,numneu};
  std::fill(inistate.begin(),inistate.end(),0);

  marray<double,1> cos_range = nus_atm.GetCosthRange();
  marray<double,1> e_range = nus_atm.GetERange();
  assert(input_flux.extent(0) == nus_atm.GetNumCos()*nus_atm.GetNumE()); // CHECK THAT THE INPUT FLUX SHAPE MATCHES NUSQUIDS OBJECT INPUT

  for (int ci = 0 ; ci < nus_atm.GetNumCos(); ci++){
    for (int ei = 0 ; ei < nus_atm.GetNumE(); ei++){
    double enu = e_range[ei]/units.GeV;
    double cth = cos_range[ci];
    assert( std::fabs(enu - input_flux[ci*e_range.size() + ei][1] ) < 1.e-4 ); // ENFORCE ENERGY AND ZENITH STEP SIZE TO NOT MISS "FAST OSCILLATION" BEHAVIOR
    assert( std::fabs(cth- input_flux[ci*e_range.size() + ei][0] ) < 1.e-4 );
    inistate[ci][ei][0][0] = input_flux[ci*e_range.size() + ei][2];
    inistate[ci][ei][1][0] = input_flux[ci*e_range.size() + ei][3];
    }
  }

  nus_atm.Set_initial_state(inistate,flavor);

  nus_atm.EvolveState(); // RUN THE EVOLUTION, SAVE THE FINAL STATE
  std::stringstream ss;
  ss<<gamma_0;
  nus_atm.WriteStateHDF5(output_path+"/atmospheric_"+ss.str()+"_"+std::to_string(n_energy)+"_"+mode+".hdf5");

  return 0;
}