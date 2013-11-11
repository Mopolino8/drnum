// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of DrNUM.                                          +
// +                                                                      +
// + Copyright 2013 numrax GmbH, enGits GmbH                              +
// +                                                                      +
// + DrNUM is free software: you can redistribute it and/or modify        +
// + it under the terms of the GNU General Public License as published by +
// + the Free Software Foundation, either version 3 of the License, or    +
// + (at your option) any later version.                                  +
// +                                                                      +
// + DrNUM is distributed in the hope that it will be useful,             +
// + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
// + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
// + GNU General Public License for more details.                         +
// +                                                                      +
// + You should have received a copy of the GNU General Public License    +
// + along with DrNUM. If not, see <http://www.gnu.org/licenses/>.        +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include <unistd.h>

#include <vtkMultiBlockDataSet.h>
#include <vtkXMLMultiBlockDataWriter.h>

#include <QFile>

#include "patchgrid.h"
#include "stringtools.h"
#include "geometrytools.h"

PatchGrid::PatchGrid(size_t num_seeklayers, size_t num_addprotectlayers)
{
  // patchgroups of same type and solver codes
  m_PatchGroups = new PatchGroups(true); /// @todo need a handling for true or false

  // default overlap settings
  m_NumSeekLayers = num_seeklayers;
  m_NumAddProtectLayers = num_addprotectlayers;

  // initializations
  m_NumFields = 0;
  m_NumVariables = 0;
  m_InterpolateData = false;
  m_TransferType = "error";
  m_BboxOk = false;
}


void  PatchGrid::setNumberOfFields(size_t num_fields)
{
  m_NumFields = num_fields;
}


void  PatchGrid::setNumberOfVariables(size_t num_variables)
{
  m_NumVariables = num_variables;
}


void PatchGrid::setInterpolateData(bool interpolatedata)
{
  m_InterpolateData = interpolatedata;
}


void PatchGrid::defineVectorVar(const size_t& index_x)
{
  /// @todo No bounds or mem overlay checking !!!
  m_VectorVarIndices.push_back(index_x);
}


vector<size_t> PatchGrid::getVectorVarIndices()
{
  return m_VectorVarIndices;
}


void PatchGrid::setNumSeekLayers(size_t num_seeklayers)
{
  m_NumSeekLayers = num_seeklayers;
}


void PatchGrid::setNumAddProtectLayers(size_t num_addprotectlayers)
{
  m_NumAddProtectLayers = num_addprotectlayers;
}


void PatchGrid::setTransferType(string trans_type)
{
  if (trans_type != "ws" &&
      trans_type != "padded" &&
      trans_type != "padded_direct") {
    BUG;
  }

  m_TransferType = trans_type;

  if (m_TransferType == "padded" || m_TransferType == "padded_direct") {
    m_TransferPadded = true;
  }
}



PatchGrid::~PatchGrid()
{
  /// @todo missing destructor. Might block mem needed for later computations.
  //  delete m_patchdependencies;  /// @todo delete method available for TInsectionList?
  //  delete m_patches;
  //  delete m_HashRaster;
}


void PatchGrid::computeDependencies(const bool& with_intercoeff)
{
  /// @todo Some performance issues in this member function. Might be optimized, if needed.

  /// @todo Needs an error handling mechanism, at least a diagnose, if any receiving cell of a patch finds no donor at all.

  /** @todo Hash search relates to outer nodes of the patches rather than to the outer cell centers. This is not
    * a problem, but might causes some extra neighbour checks, see below 4) since the nodal range is larger. */

  // Do the following
  //
  // 1) Build up a hash raster VectorHashRaster<size_t> as background grid, covering the whole geometric
  //    region of the mesh.
  //    PERF-NOTE: a) Most likely Splittree or Octree search algorithm would enhance performace, if ever needed.
  //               b) VectorHashRaster might be slow, due to large vector<vector<T> > implementation.
  //
  // 2) Fill in patch indicees. Save side inclusion: mark all hash raster cells covered by the bounding box
  //    of respective patches.
  //    PERF-NOTE: Very save side inclusion of patches on boxes of hash raster. might be improved, if needed.
  //
  // 3) Find potential neighbours on hash raster. If several patches have hit the same raster box, these patches
  //    are potentially depenedent from each other. Store potential dependencies on pot_neigh.
  //
  // 4) Find real neighbour dependencies (e.g. interpolation partners) for all potentially dependend patches.
  //    - If a potential neighbour-dependency does not serve any interpol request, it will be excluded
  //      from Patch::m_neighbours.
  //
  // 5) Finish building up inter-patch transfer lists.
  //
  // 6) Write a log file for patch dependencies

  if(!with_intercoeff) {
    // with_intercoeff==false is an intended option for grid gen purposes to be used later.
    BUG;
  }

  for (size_t i_p = 0; i_p < m_Patches.size(); i_p++) {
    m_Patches[i_p]->extractSeekCells();
  }

  vector<vector<size_t> > pot_neigh;
  pot_neigh.resize(m_Patches.size());

  // 1) Build up a hash raster VectorHashRaster<size_t> as background grid, covering the whole geometric
  //    region of the mesh.
  { // start mem_block
    VectorHashRaster<size_t> m_HashRaster;
    buildHashRaster(10*m_Patches.size(), true,   // resolution chosen upon testing, may also use fixed size
                    m_HashRaster);

    // 2) Fill in patch indicees. Save side inclusion: mark all hash raster cells covered by the bounding box
    //    of respective patches.
    for (size_t i_p = 0; i_p < m_Patches.size(); i_p++) {
      //.. Find cell address range in hash raster covered by bounding box of patch
      //.... min-max in inertial coords
      vec3_t xyzo_min = m_Patches[i_p]->accessBBoxXYZoMin();
      vec3_t xyzo_max = m_Patches[i_p]->accessBBoxXYZoMax();
      //.... min-max in m_HashRaster coords in
      vec3_t xyz_min = m_HashRaster.getTransformI2T()->transform(xyzo_min);
      vec3_t xyz_max = m_HashRaster.getTransformI2T()->transform(xyzo_max);
      size_t ic_min, jc_min, kc_min;
      size_t ic_max, jc_max, kc_max;
      bool inside_min, inside_max;
      inside_min = m_HashRaster.xyzToRefNode(xyz_min[0], xyz_min[1], xyz_min[2],
                                             ic_min, jc_min, kc_min);
      inside_max = m_HashRaster.xyzToRefNode(xyz_max[0], xyz_max[1], xyz_max[2],
                                             ic_max, jc_max, kc_max);
      //#ifdef DEBUG
      // Folowing should never happen, check anyway
      if(!inside_min || !inside_max) {
        BUG;
      }
      //#endif
      //.. Insert patch index into all boxes of hash raster covered by bounding box of patch
      for (size_t ic_r = ic_min; ic_r <= ic_max; ic_r++) {
        for (size_t jc_r = jc_min; jc_r <= jc_max; jc_r++) {
          for (size_t kc_r = kc_min; kc_r <= kc_max; kc_r++) {
            m_HashRaster.insert(ic_r, jc_r, kc_r, i_p);
          }
        }
      }
    }

    // 3) Find potential neighbours on hash raster. If several patches have hit the same raster box, these patches
    //    are potentially depenedent from each other. Store potential dependencies on pot_neigh.
    // vector<vector<size_t> > pot_neigh;
    // pot_neigh.resize(m_patches.size());
    for (size_t l_r = 0; l_r < m_HashRaster.sizeL(); l_r++) {
      size_t n_patches_in_box = m_HashRaster.getNumItems(l_r);
      if (n_patches_in_box > 1) {
        for (size_t ll_p1 = 0; ll_p1 < (n_patches_in_box - 1); ll_p1++) {  // double loop for vice-versa insertion
          for (size_t ll_p2 = (ll_p1 + 1); ll_p2 < n_patches_in_box; ll_p2++) {
            size_t patch_1 = m_HashRaster.at(l_r, ll_p1);
            size_t patch_2 = m_HashRaster.at(l_r, ll_p2);
            if(patch_1 != patch_2) {
              pot_neigh[patch_1].push_back(patch_2);
              pot_neigh[patch_2].push_back(patch_1);
            } else {
              BUG;
            }
          }
        }
        //.. keep mem low: remove duplicates
        for (size_t ll_p = 0; ll_p < m_HashRaster.getNumItems(l_r); ll_p++) {
          size_t patch = m_HashRaster.at(l_r, ll_p);
          if(pot_neigh[patch].size() > 100) {  // try to avoid unifying over and over
            // unify(pot_neigh[patch]);
            sort(pot_neigh[patch].begin(), pot_neigh[patch].end());
            typename vector<size_t>::iterator it;
            it = unique(pot_neigh[patch].begin(), pot_neigh[patch].end());
            pot_neigh[patch].resize(it - pot_neigh[patch].begin());
          }
        }
      }
    }
  } // end mem_block
  // unify all potential neighbour lists
  for (size_t patch = 0; patch < m_Patches.size(); patch++) {
    sort(pot_neigh[patch].begin(), pot_neigh[patch].end());
    typename vector<size_t>::iterator it;
    it = unique(pot_neigh[patch].begin(), pot_neigh[patch].end());
    pot_neigh[patch].resize(it - pot_neigh[patch].begin());
    vector<size_t>(pot_neigh[patch]).swap(pot_neigh[patch]);
  }
  // 4) Find real neighbour dependencies (e.g. interpolation partners) for all potentially dependend patches.
  //    - If a potential neighbour-dependency does not serve any interpol request, it will be excluded
  //      from Patch::m_neighbours.
  for (size_t i_p = 0; i_p < m_Patches.size(); i_p++) {
    for (size_t ii_pn = 0; ii_pn < pot_neigh[i_p].size(); ii_pn++) {
      size_t i_pn = pot_neigh[i_p][ii_pn];
      m_Patches[i_p]->insertNeighbour(m_Patches[i_pn]);
    }
  }
  // 5) Finish building up inter-patch transfer lists.
  /// @todo Needs an error handling mechanism, at least a diagnose, if any receiving cell of a patch finds no donor at all.
  finalizeDependencies();
  m_DependenciesOk = true;

  // 6) Write a log file for patch dependencies
  writeDependenciesLog();
}


//void PatchGrid::checkSeekSuccess()
//{
//  for (size_t i_p = 0; i_p < m_Patches.size(); i_p++) {
//    m_Patches[i_p]->checkSeekSuccess();
//  }
//}


void PatchGrid::finalizeDependencies()
{
  for (size_t i_p = 0; i_p < m_Patches.size(); i_p++) {
    m_Patches[i_p]->finalizeDependencies();
  }
}


void PatchGrid::writeDependenciesLog(string dependencies_log)
{
  QFile log_file(dependencies_log.c_str());
  log_file.open(QIODevice::WriteOnly);
  QTextStream log(&log_file);
  log.setFieldWidth(0);

  // Number of patches
  log << "Number patches: " <<  m_Patches.size() << endl;

  // Transfer type
  log << "Transfer type:  " << m_TransferType.c_str() << endl;

  log << endl;

  // averaging loop over patches
  real total_stride = 0;
  size_t num_total = 0;
  size_t max_exchange = 0;
  size_t min_exchange = 10000000;
  size_t num_exchange = 0;
  for (size_t i_p = 0; i_p < m_Patches.size(); i_p++) {
    Patch* patch = m_Patches[i_p];
    for (size_t ii_n = 0; ii_n < patch->accessNumNeighbours(); ii_n++) {
      Patch* neighbour = patch->accessNeighbour(ii_n);

      //.. Number of cells to receive and serve (vice versa)
      bool vice_exists;
      size_t num_receiving;
      size_t receive_stride;
      bool versa_exists;
      size_t num_serving;
      size_t serve_stride;
      if (m_TransferType == "padded_direct") {
        patch->diagnoseViceVersaDependencies(neighbour,
                                             vice_exists, num_receiving, receive_stride,
                                             versa_exists, num_serving, serve_stride,
                                             true);
        total_stride += num_receiving*receive_stride;
        num_total += num_receiving;
        min_exchange = min(min_exchange, num_receiving);
        max_exchange = max(max_exchange, num_receiving);
        ++num_exchange;
      }
    }
  }
  log << "total number of cell transfers : " << num_total << endl;
  log << "average strid                  : " << total_stride/num_total << endl;
  log << "smallest exchange list         : " << min_exchange << endl;
  log << "longest exchange list          : " << max_exchange << endl;
  log << "average exchange list length   : " << num_total/max(size_t(1), num_exchange) << endl << endl;

  // output loop over patches
  for (size_t i_p = 0; i_p < m_Patches.size(); i_p++) {
    Patch* patch = m_Patches[i_p];
    log.setFieldWidth(0);
    log << "Patch " << i_p << ",  type: " << patch->accessPatchType() << endl;
    log << "  Number of neighbours: " << patch->accessNumNeighbours() << endl;

    log.setFieldWidth(10);
    log << "Neighbour"
        << "ID"
        << "type"
        << "receiving"
        << "stride"
        << "serving"
        << "stride"
        << endl;

    for (size_t ii_n = 0; ii_n < patch->accessNumNeighbours(); ii_n++) {
      Patch* neighbour = patch->accessNeighbour(ii_n);
      size_t n_id = patch->accessNeighbourIndex(ii_n);

      //.. Number of cells to receive and serve (vice versa)
      bool vice_exists;
      size_t num_receiving;
      size_t receive_stride;
      bool versa_exists;
      size_t num_seving;
      size_t serve_stride;
      if (m_TransferType == "padded_direct") {
        patch->diagnoseViceVersaDependencies(neighbour,
                                             vice_exists, num_receiving, receive_stride,
                                             versa_exists, num_seving, serve_stride,
                                             true);
      }
      log << ii_n
          << n_id
          << neighbour->accessPatchType()
          << num_receiving
          << receive_stride
          << num_seving
          << serve_stride
          << endl;
    }
    log << endl;
  }

}


void PatchGrid::findBoxOverlappingPatches(const vec3_t& cbox_min, const vec3_t& cbox_max,
                                          const bool& only_core,
                                          vector<size_t>& overlap_patches)
{
  overlap_patches.resize(0);
  for (size_t i_p = 0; i_p < m_Patches.size(); i_p++) {
    if (m_Patches[i_p]->checkBoxOverlap(cbox_min, cbox_max, only_core)) {
      overlap_patches.push_back(i_p);
    }
  }
}


void PatchGrid::accessAllDonorData(const size_t& field)
{
  if(m_TransferType == "ws") {
    accessAllDonorData_WS(field);
  }
  else if (m_TransferType == "padded") {
    accessAllDonorDataPadded(field, false);
  }
  else if (m_TransferType == "padded_direct") {
    accessAllDonorDataPadded(field, true);
  }
}


void PatchGrid::accessAllDonorData_WS(const size_t& field)
{
  for (size_t i_p = 0; i_p < m_Patches.size(); i_p++) {
    m_Patches[i_p]->accessDonorData_WS(field);
  }
}


void PatchGrid::accessAllDonorDataPadded(const size_t& field, const bool& direct)
{
  /** @todo There is a slight erraneous recursion in the case of protection
    * exceptions, allowing a vice-versa interpolation access. Due to this, it is
    * not possible to do sh-mem parallelisation on the loop for patches.
    */
  if (direct) { // use direct lists
    for (size_t i_p = 0; i_p < m_Patches.size(); i_p++) {
      m_Patches[i_p]->accessDonorDataDirect(field);
    }
  }
  else { // use m_InterCoeffData
    for (size_t i_p = 0; i_p < m_Patches.size(); i_p++) {
      m_Patches[i_p]->accessDonorDataPadded(field);
    }
  }
}

void PatchGrid::buildHashRaster(size_t resolution, bool force,
                                VectorHashRaster<size_t>& m_HashRaster)
{
  // Build a bounding box
  buildBoundingBox(force);

  // Define resolution. Aim approx. same delatas in x, y, z .
  vec3_t delta_xyzo = m_BboxXyzoMax - m_BboxXyzoMin;
  real delta_resolve = pow(real(resolution) / (delta_xyzo[0] * delta_xyzo[1] * delta_xyzo[2]), (1./3.));
  real resolve_xo = delta_resolve * delta_xyzo[0];
  real resolve_yo = delta_resolve * delta_xyzo[1];
  real resolve_zo = delta_resolve * delta_xyzo[2];
  size_t i_np = size_t(resolve_xo);
  size_t j_np = size_t(resolve_yo);
  size_t k_np = size_t(resolve_zo);
  if(i_np < 1) {i_np = 1;}
  if(j_np < 1) {j_np = 1;}
  if(k_np < 1) {k_np = 1;}
  m_HashRaster.setUp(m_BboxXyzoMin[0], m_BboxXyzoMin[1], m_BboxXyzoMin[2],
                     m_BboxXyzoMax[0], m_BboxXyzoMax[1], m_BboxXyzoMax[2],
                     i_np, j_np, k_np);

}


size_t PatchGrid::insertPatch(Patch* new_patch)
{
  // Hand over general attributes
  setGeneralAttributes(new_patch);
  // Insert in list
  m_Patches.push_back(new_patch);
  m_DependenciesOk = false;  /// @todo global logics on dependencies?
  return m_Patches.size() - 1;
}


void PatchGrid::setGeneralAttributes(Patch* patch)
{
  patch->setNumberOfFields(m_NumFields);
  patch->setNumberOfVariables(m_NumVariables);
  patch->setNumSeekLayers(m_NumSeekLayers);
  patch->setNumAddProtectLayers(m_NumAddProtectLayers);
  patch->setInterpolateData(m_InterpolateData);
  patch->setTransferPadded(m_TransferPadded);
}


void PatchGrid::readGrid(string gridfilename, real scale)
{
  // Get file path
  /** @todo Use a better file path definition. */
  char base_files_path[4096];
  char grid_file[4096];
  if(getcwd(base_files_path, sizeof(base_files_path)) == NULL) {
    cout << " Could not get current working directory" << endl;
    BUG;
  };
  strcpy(grid_file, base_files_path);
  strcat(grid_file, "/");
  strcat(grid_file, gridfilename.c_str());

  // Open grid file
  ifstream s_grid;
  s_grid.open(grid_file);
  if (s_grid.fail()) {
    cout << " faild to open grid file:" << grid_file << endl;
    BUG;
  };

  // Say something
  cout << "Reading PatchGrid::readGrid() from file " << grid_file << " ... ";
  /** @todo Preliminary format. */

  // Read file contents
  while(!s_grid.eof()) {
    //.. Read int-code of patch
    size_t patch_type;
    s_grid >> patch_type;
    // End marker
    if(patch_type == 0) {
      break;
    }
    //.. Read contents for patch between delimiting { }
    char c;
    string patchcomment;
    while(s_grid.get(c)) {
      if(c=='{') {
        break;
      }
      else {
        patchcomment.push_back(c);
      }
    }
    string line;
    {
      string read_line;
      getline(s_grid, read_line, '}');
      line.clear();
      for (size_t i_read_line = 0; i_read_line < read_line.size(); ++i_read_line) {
        char c = read_line[i_read_line];
        if (isspace(c)) {
          c = ' ';
        }
        line.push_back(c);
      }
    }
    istringstream iss(line);
    //.. Hand over to patches as needed
    //.... CartesianPatch
    if (patch_type == 1001) {
      //...... Create a new CartesianPatch
      CartesianPatch* new_patch;
      new_patch = new CartesianPatch(this);
      size_t index = insertPatch(new_patch);
      new_patch->setIndex(index);
      new_patch->setPatchComment(patchcomment);
      new_patch->readFromFile(iss, scale);
      m_PatchGroups->insertPatch(new_patch); /// @todo better outside of reading loop?
    }
    //.... Unstructured Patch
    else if(patch_type == 1010) {
      // ...
    }
    else {
      cout << "unknown patch type code" << patch_type << endl;
      BUG;
    }
  }
  cout << "done. " << endl;
}

void PatchGrid::writeData(size_t i_field, QString base_file_name, real time, int count)
{

  QString count_txt;
  count_txt.setNum(count);
  count_txt = count_txt.rightJustified(6, '0');
  QString file_name = base_file_name + "_" + count_txt + ".dnd";
  QFile file(file_name);
  file.open(QIODevice::WriteOnly);
  QDataStream stream(&file);
  if (sizeof(real) == 4) {
    stream.setFloatingPointPrecision(QDataStream::SinglePrecision);
    stream << int(4);
  } else {
    stream.setFloatingPointPrecision(QDataStream::DoublePrecision);
    stream << int(8);
  }
  stream << time;
  for (size_t i_p = 0; i_p < getNumPatches(); i_p++) {
    getPatch(i_p)->writeData(i_field, stream);
  }
}

real PatchGrid::readData(size_t i_field, QString file_name)
{
  if (file_name.right(4) != ".dnd") {
    file_name += ".dnd";
  }
  QFile file(file_name);
  file.open(QIODevice::ReadOnly);
  QDataStream stream(&file);
  int prec;
  stream >> prec;
  if (prec == 4) {
    stream.setFloatingPointPrecision(QDataStream::SinglePrecision);
  } else {
    stream.setFloatingPointPrecision(QDataStream::DoublePrecision);
  }
  real time;
  stream >> time;
  for (size_t i_p = 0; i_p < getNumPatches(); i_p++) {
    getPatch(i_p)->readData(i_field, stream);
  }
  return time;
}


void PatchGrid::scaleRefParental(real scfactor)
{
  for (vector<Patch*>::iterator p = m_Patches.begin(); p != m_Patches.end(); p++) {
    (*p)->scaleRefParental(scfactor);
  }
}


void PatchGrid::buildBoundingBox(const bool& force)
{
  if(!m_BboxOk || force) {
    if(m_Patches.size() > 0) {
      // Find max-min limits of all patches in this grid
      m_BboxXyzoMin = m_Patches[0]->accessBBoxXYZoMin();
      m_BboxXyzoMax = m_Patches[0]->accessBBoxXYZoMax();
      for (size_t i_p = 1; i_p < m_Patches.size(); i_p++) {
        vec3_t bbmin_h = m_Patches[i_p]->accessBBoxXYZoMin();
        vec3_t bbmax_h = m_Patches[i_p]->accessBBoxXYZoMax();
        m_BboxXyzoMin.minimisePerCoord(bbmin_h);
        m_BboxXyzoMax.maximisePerCoord(bbmax_h);
        //        if(m_BboxXyzoMin[0] > bbmin_h[0]) m_BboxXyzoMin[0] = bbmin_h[0];
        //        if(m_BboxXyzoMin[1] > bbmin_h[1]) m_BboxXyzoMin[1] = bbmin_h[1];
        //        if(m_BboxXyzoMin[2] > bbmin_h[2]) m_BboxXyzoMin[2] = bbmin_h[2];
        //        if(m_BboxXyzoMax[0] < bbmax_h[0]) m_BboxXyzoMax[0] = bbmax_h[0];
        //        if(m_BboxXyzoMax[1] < bbmax_h[1]) m_BboxXyzoMax[1] = bbmax_h[1];
        //        if(m_BboxXyzoMax[2] < bbmax_h[2]) m_BboxXyzoMax[2] = bbmax_h[2];
      }
      m_BboxOk = true;
    }
  }
}


void PatchGrid::setFieldToConst(size_t i_field, real *var)
{
  for (size_t i_p = 0; i_p < m_Patches.size(); i_p++) {
    m_Patches[i_p]->setFieldToConst(i_field, var);
  }
}


SinglePatchGroup* PatchGrid::getSinglePatchGroup(const size_t& ipg)
{
  return m_PatchGroups->accessSinglePatchGroup(ipg);
}

real PatchGrid::computeMinChLength()
{
  real min_ch_len_all = 0;
  bool first = true;
  for (size_t i_p = 0; i_p < m_Patches.size(); i_p++) {
    real min_ch_len = m_Patches[i_p]->computeMinChLength();
    if(first || min_ch_len<min_ch_len_all) {
      min_ch_len_all = min_ch_len;
    }
    first = false;
  }
  return min_ch_len_all;
}

void PatchGrid::writeToVtk(size_t i_field, string file_name, const PostProcessingVariables &proc_vars, int count)
{
  // Synchronise blocks before saving to ensure correct values in overlap
  /// @todo field handling needed. New field required, "0 = new" OK?
  accessAllDonorData(0);
  using namespace StringTools;
  vtkSmartPointer<vtkMultiBlockDataSet> multi_block = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  multi_block->SetNumberOfBlocks(getNumPatches());
  vector<vtkSmartPointer<vtkDataSet> > data_sets(getNumPatches());
  for (size_t i_patch = 0; i_patch < getNumPatches(); ++i_patch) {
    data_sets[i_patch] = getPatch(i_patch)->createVtkDataSet(i_field, proc_vars);
    multi_block->SetBlock(i_patch, data_sets[i_patch]);
  }
  if (count >= 0) {
    file_name += "_" + leftFill(toString(count), '0', 6);
  }
  file_name += ".vtm";
  vtkSmartPointer<vtkXMLMultiBlockDataWriter> vmb = vtkSmartPointer<vtkXMLMultiBlockDataWriter>::New();
  vmb->SetInput(multi_block);
  vmb->SetFileName(file_name.c_str());
  vmb->Write();
}

bool PatchGrid::findCell(vec3_t xo, int &id_patch, int &id_cell)
{
  id_patch = -1;
  id_cell = -1;
  real min_distance = MAX_REAL;
  for (size_t id_patch = 0; id_patch < getNumPatches(); ++id_patch) {
    int id_cell_patch = getPatch(id_patch)->findCell(xo);
    if (id_cell_patch) {
      vec3_t xo_cell = getPatch(id_patch)->xyzoCell(id_cell_patch);
      real distance = (xo - xo_cell).abs();
      if (distance < min_distance) {
        min_distance = distance;
        id_patch = id_patch;
        id_cell = id_cell_patch;
      }
    }
  }
  if (id_patch < 0) {
    return false;
  }
  return true;
}
