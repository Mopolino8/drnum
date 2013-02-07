#ifndef INTERPOLATION1_H
#define INTERPOLATION1_H

#include "math/coordtransformvv.h"
#include "patchgrid.h"
#include "cartesianpatch.h"

void run()
{
  // Testing transformations (debug only)
  CoordTransformVV ct;
  vec3_t b(1., 0., 0.);
  ct.setVector(b);
  vec3_t xyzo(0., 0., 0.);
  vec3_t xyz = ct.transform(xyzo);
  vec3_t xyzo_1 = ct.transformReverse(xyz);

  // Testing patchGrids
  PatchGrid patchGrid;
  patchGrid.setInterpolateData();
  patchGrid.setTransferPadded();

//  // Define and insert some patches
//  CartesianPatch patch_0;
//  CartesianPatch patch_1;
//  CartesianPatch patch_2;
//  patch_0.setupAligned(0.0, 0.0, 0.0, 10.0, 10.0, 10.0);
//  patch_0.resize(10,10,10);
//  patch_1.setupAligned(8.0, 0.0, 0.0, 18.0, 10.0, 10.0);
//  patch_1.resize(10,10,10);
//  patch_2.setupAligned(4.0, 9.0, 0.0, 14.0, 19.0, 10.0);
//  // patch_2.setupAligned(8.0, 0.0, 0.0, 18.0, 10.0, 10.0); // identical to patch_1
//  patch_2.resize(10,10,10);
//  patchGrid.insertPatch(&patch_0);
//  patchGrid.insertPatch(&patch_1);
//  patchGrid.insertPatch(&patch_2);

  // Define and insert lots of patches
  size_t nn_i = 20;
  size_t nn_j = 20;
  size_t nn_k = 20;
  CartesianPatch* patches;
  patches = new CartesianPatch[nn_i*nn_j*nn_k];
  for (size_t i = 0; i < nn_i; i++) {
    for (size_t j = 0; j < nn_j; j++) {
      for (size_t k = 0; k < nn_k; k++) {
        size_t l = i*nn_j*nn_k + j*nn_k + k;
        real x_low = 5*i;
        real y_low = 5*j;
        real z_low = 5*k;
        patches[l].setupAligned(x_low, y_low, z_low, x_low+10.0, y_low+10.0, z_low+10.0);
        patches[l].resize(10,10,10);
        patchGrid.insertPatch(&patches[l]);
      }
    }
  }

  // Compute dependencies
  patchGrid.computeDependencies(true);

//  patch_0.insertNeighbour(&patch_1);  /// will be automatic later
//  patch_0.insertNeighbour(&patch_2);  /// will be automatic later
//  patch_1.insertNeighbour(&patch_0);  /// will be automatic later
//  patch_1.insertNeighbour(&patch_2);  /// will be automatic later
//  patch_2.insertNeighbour(&patch_0);  /// will be automatic later
//  patch_2.insertNeighbour(&patch_1);  /// will be automatic later
//  patch_0.setTransferPadded();        /// will be envoqued by PatchGrid later
//  patch_1.setTransferPadded();        /// will be envoqued by PatchGrid later
//  patch_2.setTransferPadded();        /// will be envoqued by PatchGrid later
//  patch_0.finalizeDependencies();     /// will be automatic later
//  patch_1.finalizeDependencies();     /// will be automatic later
//  patch_2.finalizeDependencies();     /// will be automatic later
}



#endif // INTERPOLATION1_H
