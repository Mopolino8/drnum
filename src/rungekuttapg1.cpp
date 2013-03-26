#include "rungekuttapg1.h"

RungeKuttaPG1::RungeKuttaPG1()
{
    m_PatchGrid = NULL;
}

void RungeKuttaPG1::operator()(real dt)
{
    /** @todo Test version only: hard coded patch interactions. */

    // Prime transfer before copying new->old
    if(m_PatchGrid) {
        for (vector<size_t>::iterator sf = m_SyncField.begin(); sf != m_SyncField.end(); sf++) {
            m_PatchGrid->accessAllDonorData_WS(*sf);
        }
    }

    // copy new->old
    copyField(0, 1);

    // stage loop
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
