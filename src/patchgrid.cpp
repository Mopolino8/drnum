#include "patchgrid.h"

PatchGrid::PatchGrid(size_t num_protectlayers, size_t num_overlaplayers)
{
  m_NumProtectLayers = num_protectlayers;
  m_NumOverlapLayers = num_overlaplayers;
  m_bbox_OK = false;
  //  m_hashBox = NULL;
}

void PatchGrid::setInterpolateData(bool interpolatedata)
{
  m_InterpolateData = interpolatedata;
}


void PatchGrid::setInterpolateGrad1N(bool interpolategrad1N)
{
  m_InterpolateGrad1N = interpolategrad1N;
}


void PatchGrid::setTransferPadded(bool trans_padded)
{
  m_TransferPadded = trans_padded;
}


void PatchGrid::setNumProtectLayers(size_t num_protectlayers)
{
  m_NumProtectLayers = num_protectlayers;
}


void PatchGrid::setNumOverlapLayers(size_t num_overlaplayers)
{
  m_NumOverlapLayers = num_overlaplayers;
}


PatchGrid::~PatchGrid()
{
  //  delete m_patchdependencies;  /// @todo delete method available for TInsectionList?
  //  delete m_patches;
  //  delete m_HashRaster;
}


void PatchGrid::computeDependencies(const bool& with_intercoeff)
{
  /// @todo function (or a similar one) might frequently be used for mesh generators later. Be fast!

  /** @todo I dislike the structure somehow, since there is no way for finding neighbours
   *  other than by checking connections on a patch-cell level.   */

  // Do the following
  //
  // 1) Build up a Splittree (or Octree or HashRaster), covering the whole mesh. Most likely a
  //    hash raster is OK, building a background grid.
  //    - Attach TInsectionList<size_t> to cells of background grid.
  //    - Fill in patch indicees in TInsectionList<size_t>. Save-side: use outer contours of
  //      patches as "belongs"-criterion. ATTENTION: must clear, if this is saveside.
  // 2) Loop through cells to find the ones that potentially overlap.
  //    - If more than one patch belongs to the same cell of the background grid, they are
  //      potentially dependent to each other. Store potential dependencies on patch::m_neighbours.
  // 3) Search interconects (e.g. interpolation partners) for all outer layer cells of all patches.
  //    - If any outer layer cells finds no donor at all: Error condition on execution programs.
  //    - If a potential neighbour-dependency does not serve any interpol request, erase it from
  //      patch::m_neighbours.

  // 1) Build a HashRaster for the mesh and fill in
  buildHashRaster(10000, true);

  //.. fill in
  for (size_t i_p = 0; i_p < m_patches.size(); i_p++) {
    // Save side inclusion: check region of bounding box of m_patches[l] only
    /// @todo If this is a performance issue later, do more severe filtering here
    //.. Find cell address range in hash raster covered by bounding box of patch
    vec3_t xyzo_min = m_patches[i_p]->accessBBoxXYZoMin();
    vec3_t xyzo_max = m_patches[i_p]->accessBBoxXYZoMax();
    size_t ic_min, jc_min, kc_min;
    size_t ic_max, jc_max, kc_max;
    bool inside_min, inside_max;
    inside_min = m_HashRaster.xyzToRefCell(xyzo_min[0], xyzo_min[1], xyzo_min[2],
                                           ic_min, jc_min, kc_min);
    inside_max = m_HashRaster.xyzToRefCell(xyzo_max[0], xyzo_max[1], xyzo_max[2],
                                           ic_max, jc_max, kc_max);
    //#ifdef DEBUG
    // Folowing should never happen, check anyway
    if(!inside_min || !inside_max) {
      BUG;
    }

    //#endif
    //.. Insert patch index into all cells of hash raster covered by bounding box of patch
    for (size_t ic_r = ic_min; ic_r <= ic_max; ic_r++) {
      for (size_t jc_r = jc_min; jc_r <= jc_max; jc_r++) {
        for (size_t kc_r = kc_min; kc_r <= kc_max; kc_r++) {
          m_HashRaster.insert(ic_r, jc_r, kc_r, i_p);
        }
      }
    }
  }

  // Find potential neighbours on hash raster.
  // If several patches have hit the same raster cell, these are potential neighbours
  vector<vector<size_t> > pot_neigh;
  pot_neigh.resize(m_patches.size());
  for (size_t l_r = 0; l_r < m_HashRaster.sizeL(); l_r++) {
    for (size_t ll_p1 = 0; ll_p1 < m_HashRaster.getNumItems(l_r); ll_p1++) {  // double loop for vice-versa insertion
      for (size_t ll_p2 = 0; ll_p2 < m_HashRaster.getNumItems(l_r); ll_p2++) {
        if(ll_p1 != ll_p2) {
          pot_neigh[ll_p1].push_back(ll_p2);
          pot_neigh[ll_p2].push_back(ll_p1);
        }
      }
    }
    // Keep mem low: remove duplicates
    for (size_t ll_p = 0; ll_p < m_HashRaster.getNumItems(l_r); ll_p++) {
      // unify(pot_neigh[ll_p]);
      sort(pot_neigh[ll_p].begin(), pot_neigh[ll_p].end());
      typename vector<size_t>::iterator it;
      it = unique(pot_neigh[ll_p].begin(), pot_neigh[ll_p].end());
      pot_neigh[ll_p].resize(it - pot_neigh[ll_p].begin());
    }
  }
  // Find true dependencies
  for (size_t i_p = 0; i_p < m_patches.size(); i_p++) {
    for (size_t ii_pn = 0; ii_pn < pot_neigh[i_p].size(); ii_pn++) {
      size_t i_pn = pot_neigh[i_p][ii_pn];
      m_patches[i_p]->insertNeighbour(m_patches[i_pn]);
    }
  }

  BUG;
}


void PatchGrid::buildHashRaster(size_t resolution, bool force)
{
  // Build a bounding box
  buildBoundingBox(force);

  // Define resolution. Aim approx. same delatas in x, y, z .
  vec3_t delta_xyzo = m_bbox_xyzo_max - m_bbox_xyzo_min;
  real delta_resolve = real(resolution) / (delta_xyzo[0] * delta_xyzo[1] * delta_xyzo[2]);
  real resolve_xo = delta_resolve * delta_xyzo[0];
  real resolve_yo = delta_resolve * delta_xyzo[1];
  real resolve_zo = delta_resolve * delta_xyzo[2];
  size_t i_np = size_t(resolve_xo);
  size_t j_np = size_t(resolve_yo);
  size_t k_np = size_t(resolve_zo);
  m_HashRaster.setUp(m_bbox_xyzo_min[0], m_bbox_xyzo_min[1], m_bbox_xyzo_min[2],
                     m_bbox_xyzo_max[0], m_bbox_xyzo_max[1], m_bbox_xyzo_max[2],
                     i_np, j_np, k_np);

}


size_t PatchGrid::insertPatch(Patch* new_patch)
{
  // Hand over general attributes
  setGeneralAttributes(new_patch);
  // Insert in list
  m_patches.push_back(new_patch);
  m_dependencies_OK = false;  /// @todo global logics on dependencies?
  return m_patches.size();
}


void PatchGrid::setGeneralAttributes(Patch* patch)
{
  patch->setNumOverlapLayers(m_NumOverlapLayers);
  patch->setNumProtectLayers(m_NumProtectLayers);
  patch->setInterpolateData(m_InterpolateData);
  patch->setInterpolateGrad1N(m_InterpolateGrad1N);
  patch->setTransferPadded(m_TransferPadded);
}


void PatchGrid::deletePatch(size_t i_patch)
{
  BUG;
  /** @todo Other patches may be dependent. These dependencies must be cleared.
    * Most likely this method will never be required, unless for mesh generation purposes.
    */
}


void PatchGrid::buildBoundingBox(const bool& force)
{
  if(!m_bbox_OK || force) {
    if(m_patches.size() > 0) {
      // prime
      m_bbox_xyzo_min = m_patches[0]->accessBBoxXYZoMin();
      m_bbox_xyzo_max = m_patches[0]->accessBBoxXYZoMax();
      for (size_t i_p = 0; i_p < m_patches.size(); i_p++) { // never mind checking m_patches[0] again
        vec3_t bbmin_h = m_patches[i_p]->accessBBoxXYZoMin();
        vec3_t bbmax_h = m_patches[i_p]->accessBBoxXYZoMax();
        if(m_bbox_xyzo_min[0] > bbmin_h[0]) m_bbox_xyzo_min[0] = bbmin_h[0];
        if(m_bbox_xyzo_min[1] > bbmin_h[1]) m_bbox_xyzo_min[1] = bbmin_h[1];
        if(m_bbox_xyzo_min[2] > bbmin_h[2]) m_bbox_xyzo_min[2] = bbmin_h[2];
        if(m_bbox_xyzo_max[0] < bbmax_h[0]) m_bbox_xyzo_max[0] = bbmax_h[0];
        if(m_bbox_xyzo_max[1] < bbmax_h[1]) m_bbox_xyzo_max[1] = bbmax_h[1];
        if(m_bbox_xyzo_max[2] < bbmax_h[2]) m_bbox_xyzo_max[2] = bbmax_h[2];
        ///< @todo A max/min function for vec3_t ?
        // m_bbox_xyzo_min.coordMin(bbmin_h);
        // m_bbox_xyzo_max.coordMin(bbmax_h);
      }
      m_bbox_OK = true;
    }
  }
}


