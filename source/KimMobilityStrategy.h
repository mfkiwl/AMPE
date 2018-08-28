#ifndef included_KimMobilityStrategy
#define included_KimMobilityStrategy

#include "CALPHADFreeEnergyFunctions.h"
#include "SimpleQuatMobilityStrategy.h"

#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/hier/Patch.h"

class QuatModel;

#include <boost/make_shared.hpp>
#include <string>

using namespace SAMRAI;

/*
 * Based on S.G. Kim, Acta Mat. 55 (2007), p. 4391-4399
 */
class KimMobilityStrategy:
   public SimpleQuatMobilityStrategy
{
public:
   KimMobilityStrategy(
      QuatModel* quat_model,
      const int conc_l_id,
      const int conc_s_id,
      const int temp_id,
      const std::string& energy_interp_func_type,
      const std::string& conc_interp_func_type,
      boost::shared_ptr<tbox::Database> calphad_db,
      boost::shared_ptr<tbox::Database> newton_db,
      const unsigned ncompositions,
      const double DL,
      const double Q0);

   void computePhaseMobility(
      const boost::shared_ptr<hier::PatchHierarchy > hierarchy,
      int& phase_id,
      int& mobility_id,
      const double time,
      const CACHE_TYPE cache = CACHE );

   virtual double evaluateMobility(const double temp,
      const std::vector<double>&  phaseconc)=0;

protected:
   const int d_conc_l_id;
   const int d_conc_s_id;
   const int d_temp_id;

   const unsigned d_ncompositions;

   /*!
    * Diffusivity in liquid
    */
   const double d_DL;
   const double d_Q0;

   CALPHADFreeEnergyFunctions* d_calphad_fenergy;

private:

   void update(
      boost::shared_ptr< pdat::CellData<double> > cd_te,
      boost::shared_ptr< pdat::CellData<double> > cd_cl,
      boost::shared_ptr< pdat::CellData<double> > cd_cs,
      boost::shared_ptr< pdat::CellData<double> > cd_mobility,
      boost::shared_ptr<hier::Patch > patch );
};

#endif
