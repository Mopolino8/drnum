#include "rungekuttapg1.h"

RungeKuttaPG1::RungeKuttaPG1()
{
  m_PatchGrid = NULL;
}

void RungeKuttaPG1::operator()(real dt)
{
  copyField(0, 1);
  for (list<real>::iterator i = m_Alpha.begin(); i != m_Alpha.end(); ++i) {

    /** @todo Test version only: hard coded patch interactions. */
    if(m_PatchGrid) {
      for (vector<size_t>::iterator sf = m_SyncField.begin(); sf != m_SyncField.end(); sf++) {
        m_PatchGrid->accessAllDonorData_WS(*sf);
      }
    }

    computeIterators((*i)*dt);
    runPostOperations();
    countFlops(1);
  }
}
