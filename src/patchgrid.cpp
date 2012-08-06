#include "patchgrid.h"

PatchGrid::PatchGrid()
{
  m_patches = NULL;
  m_patchdependencies = NULL;
  m_bbox_OK = false;
}

PatchGrid::~PatchGrid()
{
  delete m_patchdependencies;  /// @todo delete method available for TInsectionList?
  delete m_patches;
}

void PatchGrid::initLists(size_t size_default, size_t size_incr)
{
  m_patches = new TList<Patch*>(size_default, size_incr);
  m_patchdependencies = new TInsectionList<Patch*>(m_patches);
}

void PatchGrid::computeDependencies(bool with_intercoeff)
{
  /// @todo function (or a similar one) might frequently be used for mesh generators later. Be fast!

  /** @todo I dislike the structure somehow, since there is no way for finding neighbours
   *  other than by checking connections on a patch-cell level.   */

  // Do the following
  //
  // 1) Build up a Splittree (or Octree or CartRaster), covering the whole mesh. Most likely a
  //    Hash raster is OK, building a background grid. May use CartesianPatch or parts of it.
  //    - Attach TInsectionList<size_t> to cells of background grid.
  //    - Fill in patch indicees in TInsectionList<size_t>. Save-side: use outer contours of
  //      patches as "belongs"-criterion. ATTENTION: must clear, if this is saveside.
  // 2) Loop through cells to find the ones that potential overlap.
  //    - If more than one patch belongs to the same cell of the background grid, they are
  //      potentially dependent to each other. Store potential dependencies on patch::m_neighbours.
  // 3) Search interconects (e.g. interpolation partners) for all outer layer cells of all patches.
  //    - If any outer layer cells finds no donor at all: Error condition on execution programs.
  //    - If a potential neighbour-dependency does not serve any interpol request, erase it from
  //      patch::m_neighbours.

  // 1) Build a HashRaster for the mesh
  buildBoundingBox();



  BUG;
}


void PatchGrid::buildBoundingBox()
{
  if(!m_bbox_OK) {
    if(m_patches->NumEntries() > 0) {
      // prime
      m_bbox_xyzo_min = m_patches->At(0)->accessBBoxXYZoMin();
      m_bbox_xyzo_max = m_patches->At(0)->accessBBoxXYZoMax();
      FORALL(i_patch, m_patches->) { // never mind checking m_patches->At(0) again
        vec3_t bbmin_h = m_patches->At(i_patch)->accessBBoxXYZoMin();
        vec3_t bbmax_h = m_patches->At(i_patch)->accessBBoxXYZoMax();
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


