#ifndef RASTER_H
#define RASTER_H

#include <cstddef>
#include "blockcfd.h"
#include "utility/weightedset.h"
#include "math/coordtransformvv.h"

//class Raster;

//#ifdef WITH_VTK
//#include <QString>
//#include <QVector>
//#include <vtkSmartPointer.h>
//#include <vtkRectilinearGrid.h>
//#include <vtkXMLRectilinearGridWriter.h>
//#include <vtkFloatArray.h>
//#include <vtkCellData.h>
//#endif

class Raster
{


protected: // attributes

  size_t m_NumL; ///< Number of cells in the raster;

  real m_Xo;  ///< xo-pos of reference point in intertial coords
  real m_Yo;  ///< yo-pos of reference point in intertial coords
  real m_Zo;  ///< zo-pos of reference point in intertial coords

  real m_Lx;   ///< size in x-direction (own coord syst)
  real m_Ly;   ///< size in y-direction (own coord syst)
  real m_Lz;   ///< size in z-direction (own coord syst)

  /// @todo change to vec3_t later
  //  SVec3 m_xyzCCMin;      ///< Cell center coords of lowest address cell
  //  SVec3 m_xyzInterMin;   ///< lower CC coord-limits to access data from foreign patches
  //  SVec3 m_xyzInterMax;   ///< upper CC coord-limits to access data from foreign patches

  /** @todo attention: m_NumProtectLayers has been defined in  I"patch".f ever "patch" will
    * be derived from this raster pattern, the following attributes will cause an inheritance conflict:
    * m_transformInertial2This, m_bbox_xyzo_min, m_bbox_xyzo_max, m_NumProtectLayers. */
  CoordTransformVV m_transformInertial2This; ///< transformation matrix to transform intertial coords into system of "this"
  vec3_t m_bbox_xyzo_min;                    ///< lowest coordinates of smallest box around patch in inertial coords.
  vec3_t m_bbox_xyzo_max;                    ///< highest coordinates of smallest box around patch in inertial coords.
  size_t m_NumProtectLayers;                 ///< Number of protection layers

  real m_Eps, m_EpsDX, m_EpsDY, m_EpsDZ;     ///< precision thresholds
  /** @todo We will need prec limits thoughout the program. It would be usefull to provide basic values depending on data
    *       type used, e.g. float basicEps = 1.e-6 and double basicEps = 1.e-12   This allows to set thresholds independent from
    *       actuial real type in use. */


protected: // methods


public: // methods

  Raster();

  virtual void resize(size_t num_l);

  void setNumProtectLayers(size_t num_protectlayers);  ///< Set number of protection layers

  size_t sizeL() const {return m_NumL;}

  /**
   * Build a bounding box in inertial coordinates around the  raster.
   * NOTE: Raster coords may be non-aligned to the xyzo-system.
   */
  virtual void buildBoundingBox()=0; // {BUG;};

  /**
   * Set up interpolation methods for giving data to foreign patches.
   * Example: Build up Split- or Octrees for search operations, etc ... depending on patch type.
   * @param num_protection number of overlap cell layers in which no data access is permissible.
   */
  virtual void setupInterpolators()=0; // {BUG;};

  /**
   * Get data interpolation coeff-sets.
   * @param x the x-value in the coords of the present patch
   * @param y the y-value in the coords of the present patch
   * @param z the z-value in the coords of the present patch
   * @param w_set WeightedSet<real> object on which to write data (return reference)
   * @return true, if interpol sets were found.
   */
  virtual bool computeCCDataInterpolCoeffs(real x, real y, real z,
                                           WeightedSet<real>& w_set)=0; // {BUG;}

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
                                             WeightedSet<real>& w_set)=0; // {BUG;}

  /**
   * Access m_transformInertial2This as pointer for extrnal use.
   */
    CoordTransformVV* getTransformI2T()
    {
      return &m_transformInertial2This;
    }
};

#endif // RASTER_H
