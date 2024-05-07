// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
#ifndef included_RigidBodyMotionRHS
#define included_RigidBodyMotionRHS

#include "SAMRAI/hier/PatchHierarchy.h"

#include <vector>
#include <array>

using namespace SAMRAI;

class RigidBodyMotionRHS
{
 public:
   RigidBodyMotionRHS(const int data_scratch_id, const int weight_id,
                      const std::vector<double>& mobilities);

   void addRHS(std::shared_ptr<hier::PatchHierarchy> hierarchy,
               const int ydot_id,
               const std::vector<std::array<double, NDIM>>& forces);

 private:
   const int d_data_id;
   const int d_weight_id;

   /*!
    * Translational mobilities relating velocities to forces
    */
   const std::vector<double> d_mobilities;

   /*!
    * Volumes of each grain
    */
   std::vector<double> d_volumes;

   void computeVolumes(std::shared_ptr<hier::PatchHierarchy> hierarchy);
};

#endif