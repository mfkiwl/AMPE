// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// Written by M.R. Dorr, J.-L. Fattebert and M.E. Wickett
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// - Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the disclaimer below.
// - Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the disclaimer (as noted below) in the
//   documentation and/or other materials provided with the distribution.
// - Neither the name of the LLNS/LLNL nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
// LLC, UT BATTELLE, LLC,
// THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
#include "CALPHADFreeEnergyFunctionsBinary.h"
#include "CALPHADMobility.h"
#include "CompositionStrategyMobilities.h"
#include "CALPHADFreeEnergyStrategyBinary.h"
#include "ConstantMolarVolumeStrategy.h"

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/Database.h"

#include <string>

using namespace SAMRAI;


int main(int argc, char *argv[])
{
   // Initialize MPI, SAMRAI

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   /* This extra code block is used to scope some temporaries that are
    * created, it forces the destruction before the manager is
    * shutdown.
    */
   {

      //-----------------------------------------------------------------------
      /*
       * Process command line arguments and dump to log file.
       * For non-restarted case, command line is:
       *
       *    executable <input file name>
       *
       */

      std::string input_filename;
      input_filename = argv[1];

      //-----------------------------------------------------------------------
      // Create input database and parse all data in input file.

      std::shared_ptr<tbox::MemoryDatabase> input_db(
          new tbox::MemoryDatabase("input_db"));
      tbox::InputManager::getManager()->parseInputFile(input_filename,
                                                       input_db);

      // make from input file name
      std::string run_name =
          input_filename.substr(0, input_filename.rfind("."));

      // Logfile
      std::string log_file_name = run_name + ".log";
      tbox::PIO::logOnlyNodeZero(log_file_name);

#ifdef GITVERSION
#define xstr(x) #x
#define LOG(x) tbox::plog << " AMPE: git version " << xstr(x) << std::endl;
      LOG(GITVERSION);
      tbox::plog << std::endl;
#endif

      tbox::plog << "input_filename = " << input_filename << std::endl;

      std::shared_ptr<tbox::Database> model_db =
          input_db->getDatabase("ModelParameters");

      double phase_well_scale = model_db->getDouble("phi_well_scale");

      std::string phase_well_func_type =
          model_db->getString("phi_well_func_type");

      std::string energy_interp_func_type = "pbg";
      std::string conc_interp_func_type = "pbg";
      std::string eta_interp_func_type = "pbg";

      std::shared_ptr<tbox::Database> temperature_db =
          model_db->getDatabase("Temperature");
      double temperature = temperature_db->getDouble("temperature");

      std::shared_ptr<tbox::Database> conc_db(
          model_db->getDatabase("ConcentrationModel"));
      std::string conc_avg_func_type =
          conc_db->getStringWithDefault("avg_func_type", "a");

      std::shared_ptr<tbox::Database> dcalphad_db =
          conc_db->getDatabase("Calphad");
      std::string calphad_filename = dcalphad_db->getString("filename");
      std::shared_ptr<tbox::MemoryDatabase> calphad_db(
          new tbox::MemoryDatabase("calphad_db"));
      tbox::InputManager::getManager()->parseInputFile(calphad_filename,
                                                       calphad_db);

      std::shared_ptr<tbox::Database> newton_db;
      if (conc_db->isDatabase("NewtonSolver"))
         newton_db = conc_db->getDatabase("NewtonSolver");

      bool with_third_phase = false;

      CALPHADFreeEnergyFunctionsBinary cafe(calphad_db, newton_db,
                                            energy_interp_func_type,
                                            conc_interp_func_type,
                                            with_third_phase);

      cafe.printEnergyVsComposition(temperature);

      cafe.preRunDiagnostics(303., 2899.);

      // choose pair of phases: phaseL, phaseA, phaseB
      const PhaseIndex pi0 = phaseL;
      const PhaseIndex pi1 = phaseA;

      // initial guesses
      double init_guess[2];
      model_db->getDoubleArray("initial_guess", &init_guess[0], 2);

      double lceq[2] = {init_guess[0], init_guess[1]};

      // compute equilibrium concentrations
      bool found_ceq = cafe.computeCeqT(temperature, pi0, pi1, &lceq[0]);
      if (lceq[0] > 1.) found_ceq = false;
      if (lceq[0] < 0.) found_ceq = false;
      if (lceq[1] > 1.) found_ceq = false;
      if (lceq[1] < 0.) found_ceq = false;

      if (found_ceq) {
         tbox::pout << "Found equilibrium concentrations: " << lceq[0]
                    << " and " << lceq[1] << "..." << std::endl;
      } else {
         tbox::pout << "WARNING: Equilibrium concentrations not found... "
                    << std::endl;
      }

      double molar_volume = conc_db->getDouble("molar_volume");

      tbox::plog << "ConstantMolarVolumeStrategy... " << std::endl;
      ConstantMolarVolumeStrategy mvstrategy(molar_volume, molar_volume,
                                             molar_volume);
      tbox::plog << "CALPHADFreeEnergyStrategy... " << std::endl;
      CALPHADFreeEnergyStrategyBinary free_energy_strategy(
          calphad_db, newton_db, energy_interp_func_type, conc_interp_func_type,
          &mvstrategy, -1, -1, -1, false);

      if (calphad_db->keyExists("MobilityParameters")) {
         tbox::plog << "CompositionStrategyMobilities... " << std::endl;
         int ncompositions = 1;
         bool with_third_phase = false;
         CompositionStrategyMobilities composition_strategy_mobilities(
             dcalphad_db, with_third_phase, ncompositions,
             &free_energy_strategy);

         composition_strategy_mobilities.printDiagnostics(temperature,
                                                          temperature);
      }
      cafe.energyVsPhiAndC(temperature, &lceq[0], found_ceq, phase_well_scale,
                           phase_well_func_type, 101, 100);

      input_db.reset();
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return (0);
}