#include "rungekuttapg1.h"

RungeKuttaPG1::RungeKuttaPG1()
{
  m_PatchGrid = NULL;
}

RungeKuttaPG1::RungeKuttaPG1(PatchGrid* patch_grid)
{
  m_PatchGrid = patch_grid;
}

void RungeKuttaPG1::operator()(real dt)
{
  /** @todo Test version only: hard coded patch interactions. */

  // Prime transfer before copying new->old
  if (m_PatchGrid) {
    for (vector<size_t>::iterator sf = m_SyncField.begin(); sf != m_SyncField.end(); sf++) {
      m_PatchGrid->accessAllDonorData(*sf);

//      // test performance
//      size_t dummi_iter = 1000;
//      for (size_t ii=0; ii < dummi_iter; ii++) {
//        m_PatchGrid->accessAllDonorData(*sf);
//      }
    }
  }

  // copy new->old
  copyField(0, 1);

  // stage loop
  bool first_step = true;
  for (list<real>::iterator i = m_Alpha.begin(); i != m_Alpha.end(); ++i) {
    /** @todo Test version only: hard coded patch interactions. */
    //if (m_PatchGrid) {
    if (m_PatchGrid && !first_step) {
      for (vector<size_t>::iterator sf = m_SyncField.begin(); sf != m_SyncField.end(); sf++) {
        m_PatchGrid->accessAllDonorData(*sf);
      }
    }
    first_step = false;
    computeIterators((*i)*dt);
    runPostOperations();
    countFlops(1);
  }

}
