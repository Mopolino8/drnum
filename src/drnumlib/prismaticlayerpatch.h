#ifndef PRISMATICLAYERPATCH_H
#define PRISMATICLAYERPATCH_H

struct Node;
template <class T> class InsectionList;
struct Prisms2Nodes;
struct BoundaryCode;
class PrismaticLayerPatch;

//#include <cstddef>
//#include <string.h>
//#include <iostream>
//#include <fstream>
//#include <sstream>
#include "patch.h"

#ifdef WITH_VTK
#include <QString>
#include <QVector>
#include <vtkSmartPointer.h>
#include <vtkRectilinearGrid.h>
#include <vtkXMLRectilinearGridWriter.h>
#include <vtkFloatArray.h>
#include <vtkCellData.h>
#endif


/** Coordinate triple for nodes.
  * @todo may be replaced by a vec3_t, if more useful.
*/
struct Node
{
  real m_X;
  real m_Y;
  real m_Z;
};


/** Test version of insection list withou any automatic feeders, etc ...
  * @todo put in utility or similar location, or use fvmose::TInsectionList<T> and upgrade
  */
template <class T>
class InsectionList
{

protected:

  /** List of starting indicees in m_ItemList */
  vector<size_t> m_Start;

  /** List to hold item data */
  vector<T> m_ItemList;

public:
  /** Get the number of entries in the second dimension.
    *  @param i index in first dimension
    *  @return the number of items stored for i
    */
  size_t numItems(const size_t& i) const {
    if(i<(m_Start.size()-1)) {
      return (m_Start[i+1] - m_Start[i]);
    } else {
      return (m_ItemList.size() - m_Start[i]);
    }
  }

  /** Data access method (read/write)
    *  @param i index in first dimension
    *  @param k index in second dimension
    *  @return a reference to the entry
    */
  T& At(const size_t& i, const size_t& k) {
    return m_ItemList[m_Start[i]+k];
  }

  /** Get the index of an entry in the item list
    *  @param i index in first dimension
    *  @param k index in second dimension
    *  @return the address of the entry
    */
  size_t ItemIndex(const size_t& i, const size_t& k) const {
    return m_Start[i]+k;
  }

  /** Get the index of the first item belonging to a set of i.
    *  @param i index in first dimension
    *  @return the address of the entry
    */
  size_t ItemIndex(const size_t& i) const {
    return m_Start[i];
  }

  /** Sequential feed of a data set. Feeds data at next free position in
    * first dimension.
    *  @param num_j number of items in second dimension to feed.
    *  @param feed_set set of data items to feed.
    *  @return index in first dimension, where data were fed in.
    */
  size_t feedDataSet(const size_t& num_j,
                     const T* feed_set) {
    size_t i = m_Start.size();
    m_Start.push_back(m_ItemList.size());
    for(size_t jj = 0; jj < num_j; jj++) {
      m_ItemList.push_back(*(feed_set+jj));
    }
    return i;
  }

  /**
    * Get a reference to m_ItemList.
    *  @return reference to m_ItemList.
    */
  vector<T>& accessItemList() {
    return m_ItemList;
  }
};

/** Index sequences for nodes of prismatic cells.
  */
struct Prisms2Nodes
{
  /** Starting index of nodes of the bottom face of the cell.*/
  size_t m_StartBottom;

  /** Starting index of nodes of the top face of the cell.*/
  size_t m_StartTop;

  /** Number of nodes forming the bottom and top face each.*/
  size_t m_NumFaceNodes;
};


/** Boundary code reference.
  * @todo all
*/
struct BoundaryCode
{
  size_t m_ProtectionException;
};


/**
  * Class to manage prismatic layer patches.
  *
  * A prismatic layer patch is semi-structured as all layers share the same 2D cell arangements, topologically
  * based on a general 2D unstructured grid.
  *
  * Explicit data:
  *  - Topological set of a 2D generally unstructured polygonial faces in a reference surface, discribed as a
  *    set of reference node indicees each. Right hand rule applies to the node sequences, pointing inwards
  *    the mesh (towards higher cell layers).
  *  - Number of cell layers, all sharing the same layer topology of the reference surface.
  *  - Node coordinates for all nodes in the 3D-mesh, starting with the nodes in the reference surface.
  *  - Boundary code reference.
  *
  * Implicit data:
  *  - Reconstruction of node indicees of cells in the layers. Note: 0-th cell layers adjacent to reference surface.
  *    Nodes forming a cell(i_layer, i_2d):
  *      - bottom face nodes of cell: 2D-node-indicees of i_2d, incremented by num_2d_nodes * i_layer
  *      - top face nodes of cell:    2D-node-indicees of i_2d, incremented by num_2d_nodes * (i_layer + 1)
  *        (note reversed right hand rule for faces)
  *
  */
class PrismaticLayerPatch : public Patch
{
  //postponed #include "prismaticlayerpatch_common.h"

private:

protected: // attributes

  //TInsectionList<size_t> m_2DNodes;  /// @todo most functions in TInsectionList not needed. Take other data type?
  /** Indicees of nodes forming the 2D-cells of the reference surface. */
  InsectionList<size_t> m_RefFacesTo2DNodes;

  /** Coordinates of all nodes (3D indexing). */
  vector<Node> m_Coords;

  /** Max Index in 2D-nodes. I.g. equal to total number of nodes per 2D layer. */
  size_t m_Max2DNodeIndex;

  /** Flag to mark m_Num2DNodes as being available. */
  bool m_Max2DNodeIndexOK;

  /** Boundary condition reference on patch. */
  vector<BoundaryCode> m_BC;

  /** Epsilon for general purposes. */
  real m_Eps;

protected: // methods

  virtual void extractReceiveCells();
  virtual void buildBoundingBox();

public: // methods


  /**
   * Constructor
   * @param num_seeklayers default number of seeking element layers
   * @param num_addprotectlayers default number of additional protection layers
   */
  PrismaticLayerPatch(PatchGrid *patch_grid, size_t num_seeklayers = 2, size_t num_addprotectlayers = 0);


  /**
    * Read mesh data from file
    * @param s_mesh the stream to read from
    * @return true, if successful
    */
  virtual bool readFromFile(istringstream& iss_input);


  /**
    * Write mesh data to file
    * @param s_mesh the stream to write to
    * @return true, if successful
    */
  virtual bool writeToFile(ifstream& s_mesh);


  /**
    * Compute number of nodes in 2D layer.
    */
  void findMax2DNodeIndex();


  /**
    * Scale patch relative to origin of parental coordsyst.
    * NOTE: Affects reference position and physical patch size.
    * Virtual: base class Patch holds reference position. Derived class holds phys. patch size.
    * @param scfactor scaling factor.
    */
  virtual void scaleRefParental(real scfactor);


  /**
    * Set a single seek exception on a boundary identified by a bcode
    * @param i_bc the boundary code for which to set the exception
    * @param num_seeklayers number of seek layer to set for i_bc
    */
  void setSeekException(const size_t& i_bc,
                        const size_t& num_seeklayers);


  /**
    * Compute inter patch dependencies to a neighbour patch.
    * @param i_neighbour neighbour patch to consider.
    * @return true, if a dependency was found.
    */
  virtual bool computeDependencies(const size_t& i_neighbour);


  /**
   * @brief Get set of node indicees forming a cell.
   * @param i_layer index of the layer.
   * @param i_2d index in unstructured reference surface.
   * @param indicees_3D the struct of type Prisms2Nodes for full 3D node indexing in a primatic grid (return reference).
   */
  inline void prisms2Nodes(const size_t& i_layer, const size_t& i_2d,
                           Prisms2Nodes& indicees_3D);

  /**
   * Set up interpolation methods for giving data to foreign patches.
   * Example: Build up Split- or Octrees for search operations, etc ... depending on patch type.
   * @param num_protection number of overlap cell layers in which no data access is permissible.
   */
  virtual void setupInterpolators();

  /**
   * Get data interpolation coeff-sets.
   * @param x the x-value in the coords of the present patch
   * @param y the y-value in the coords of the present patch
   * @param z the z-value in the coords of the present patch
   * @param w_set WeightedSet<real> object on which to write data (return reference)
   * @return true, if interpol sets were found.
   */
  virtual bool computeCCDataInterpolCoeffs(real x, real y, real z,
                                           WeightedSet<real>& w_set);

  /**
   * Get directional derivative (grad*n) interpolation coeff-sets.
   * @param x the x-value in the coords of the present patch
   * @param y the y-value in the coords of the present patch
   * @param z the z-value in the coords of the present patch
   * @param nx the x-component of directional vector in the coords of the present patch
   * @param ny the y-component of directional vector in the coords of the present patch
   * @param nz the z-component of directional vector in the coords of the present patch
   * @param w_set WeightedSet<real> object on which to write data (return reference)
   * @return true, if interpol sets were found.
   */
  virtual bool computeCCGrad1NInterpolCoeffs(real x, real y, real z,
					     real nx, real ny, real nz,
                                             WeightedSet<real>& w_set);



#ifdef WITH_VTK
  template <typename TVariables> void writeToVtk(size_t i_field, QString file_name, const TVariables& variables);
#endif

};


void PrismaticLayerPatch::prisms2Nodes(const size_t& i_layer, const size_t& i_2d,
                                       Prisms2Nodes& nodes_3D)
{
  nodes_3D.m_StartBottom = i_layer * m_Max2DNodeIndex + m_RefFacesTo2DNodes.ItemIndex(i_2d);
  nodes_3D.m_StartTop    = (i_layer + 1) * m_Max2DNodeIndex + m_RefFacesTo2DNodes.ItemIndex(i_2d);
  nodes_3D.m_NumFaceNodes = m_RefFacesTo2DNodes.numItems(i_2d);
}


//#ifdef WITH_VTK
//template <typename TVariables>
//void PrismaticLayerPatch::writeToVtk(size_t i_field, QString file_name, const TVariables& variables)
//{
//  vtkSmartPointer<vtkRectilinearGrid> grid = vtkSmartPointer<vtkRectilinearGrid>::New();

//  vtkSmartPointer<vtkFloatArray> xc = vtkSmartPointer<vtkFloatArray>::New();
//  for (size_t i = 0; i < m_NumI + 1; ++i) {
//    xc->InsertNextValue(i*dx());
//  }

//  vtkSmartPointer<vtkFloatArray> yc = vtkSmartPointer<vtkFloatArray>::New();
//  for (size_t j = 0; j < m_NumJ + 1; ++j) {
//    yc->InsertNextValue(j*dy());
//  }

//  vtkSmartPointer<vtkFloatArray> zc = vtkSmartPointer<vtkFloatArray>::New();
//  for (size_t k = 0; k < m_NumK + 1; ++k) {
//    zc->InsertNextValue(k*dz());
//  }

//  grid->SetDimensions(m_NumI + 1, m_NumJ + 1, m_NumK + 1);
//  grid->SetXCoordinates(xc);
//  grid->SetYCoordinates(yc);
//  grid->SetZCoordinates(zc);

//  real* raw_var = new real [numVariables()];
//  for (int i_var = 0; i_var < variables.numScalars(); ++i_var) {
//    vtkSmartPointer<vtkFloatArray> var = vtkSmartPointer<vtkFloatArray>::New();
//    var->SetName(qPrintable(variables.getScalarName(i_var)));
//    var->SetNumberOfValues(variableSize());
//    grid->GetCellData()->AddArray(var);
//    vtkIdType id = 0;
//    for (size_t k = 0; k < m_NumK; ++k) {
//      for (size_t j = 0; j < m_NumJ; ++j) {
//        for (size_t i = 0; i < m_NumI; ++i) {
//          getVar(i_field, i, j, k, raw_var);
//          var->SetValue(id, variables.getScalar(i_var, raw_var));
//          ++id;
//        }
//      }
//    }
//  }
//  for (int i_var = 0; i_var < variables.numVectors(); ++i_var) {
//    vtkSmartPointer<vtkFloatArray> var = vtkSmartPointer<vtkFloatArray>::New();
//    var->SetName(qPrintable(variables.getVectorName(i_var)));
//    var->SetNumberOfComponents(3);
//    var->SetNumberOfTuples(variableSize());
//    grid->GetCellData()->AddArray(var);
//    vtkIdType id = 0;
//    for (size_t k = 0; k < m_NumK; ++k) {
//      for (size_t j = 0; j < m_NumJ; ++j) {
//        for (size_t i = 0; i < m_NumI; ++i) {
//          getVar(i_field, i, j, k, raw_var);
//          vec3_t v = variables.getVector(i_var, raw_var);
//          float vf[3];
//          vf[0] = v[0]; vf[1] = v[1]; vf[2] = v[2];
//          var->SetTuple(id, vf);
//          ++id;
//        }
//      }
//    }
//  }
//  delete [] raw_var;

//  vtkSmartPointer<vtkXMLRectilinearGridWriter> vtr = vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();
//  vtr->SetFileName(qPrintable(file_name + ".vtr"));
//  vtr->SetInput(grid);
//  vtr->Write();
//}
//#endif

#endif // PRISMATICLAYERPATCH_H
