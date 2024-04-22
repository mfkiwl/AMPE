// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.

#include "RigidBodyForces.h"
#include "QuatFort.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"

using namespace SAMRAI;

RigidBodyForces::RigidBodyForces(const QuatModelParameters& model_parameters,
                                 const int phase_id, const int weight_id)
    : d_model_parameters(model_parameters),
      d_phase_id(phase_id),
      d_weight_id(weight_id)
{
   // forces between solid particles only, not liquid phase (last order parameter)
   const int nforces = d_model_parameters.norderp()-1;
   d_forces.resize(2 * nforces * nforces);
}

void RigidBodyForces::evaluatePairForces(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   for (auto& f : d_forces)
      f = 0.;

   int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           ++p) {
         std::shared_ptr<hier::Patch> patch = *p;
         evaluatePairForces(patch);
      }
   }

   sumReduction(&d_forces[0], d_forces.size());
}

void RigidBodyForces::printPairForces(std::ostream& os)
{
   os << "Grain pair forces:" << std::endl;
   int count = 0;
   for (auto& f : d_forces) {
      os << " " << f;
      count++;
      if (count % NDIM == 0) os << std::endl;
   }
   os << std::endl;
}

void RigidBodyForces::evaluatePairForces(std::shared_ptr<hier::Patch> patch)
{
   const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                              hier::PatchGeometry>(patch->getPatchGeometry()));
   const double* dx = patch_geom->getDx();

   const hier::Box& pbox = patch->getBox();
   const hier::Index& ifirst = pbox.lower();
   const hier::Index& ilast = pbox.upper();

   std::shared_ptr<pdat::CellData<double> > phase(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch->getPatchData(d_phase_id)));
   std::shared_ptr<pdat::CellData<double> > weight(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch->getPatchData(d_weight_id)));

   assert(phase);
   assert(weight);

   assert(weight->getGhostCellWidth() ==
          hier::IntVector(tbox::Dimension(NDIM), 0));

   // printf("orient_interp_func_type2 =%s\n",
   // d_model_parameters.orient_interp_func_type2().c_str()[0]);
   PHIPHI_FORCES(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                 ifirst(2), ilast(2),
#endif
                 dx, phase->getPointer(), phase->getGhostCellWidth()[0],
                 phase->getDepth()-1, d_model_parameters.gamma(),
                 d_model_parameters.m_moelans2011(), weight->getPointer(),
                 &d_forces[0]);
}
