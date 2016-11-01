
/*
 *  @license GPL
 *  @author Benjamin Lefaudeux (blefaudeux at github)
 *  @file pyUKF.cpp
 *  @brief This is a wrapper to the UKF library, exposing it to Python using Boost.Python
 *  Very basic from now, just testing the idea..
 *
 *  @version 0.1
 *  @date 25-08-2015
 */

#include <stdlib.h>
#include <iostream>

// Boost / Python bindings
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <vector>
#include <list>

#include <unscented_KF.h>

using namespace std;

/*!
 * \brief The pyUKF struct
 * A structure declaring the UKF filter,
 * which we can access from Python.
 * To be completed...
 */
struct pyUKF {
public :
  UKF *filter;

  pyUKF()
  {
      // TBD
  }

  // Destructor : stop the brackground thread first !
  inline ~pyUKF()
  {
      // TBD
  }

  // TODO: Ben
};


/***********************************/
/* Define the Python access points */
BOOST_PYTHON_MODULE(libpyUKF)
{
  using namespace boost::python;

  class_<std::vector<double> >("double_vector")
      .def(vector_indexing_suite<std::vector<double> >());

  // Declare the Python API : constructor, feed the filter, get state, etc..
//  class_<pyPTAM>("pyUKF", init<std::string, bool>())
//      .def("Start",               &pyUKF::Start)
//      .def("GetPose",             &pyPTAM::GetPose)
//      .def("GetCurrentKeyframes", &pyPTAM::GetCurrentKeyframes)
//      .def("GetCurrentPoints",    &pyPTAM::GetCurrentPoints)
//      .def("GetDiscardedPoints",  &pyPTAM::GetDiscardedPoints)
//      .def("LoadARModel",         &pyPTAM::LoadARModel)
//      .def("ResetMap",            &pyPTAM::ResetMap)
//      .def("IsAlive",             &pyPTAM::IsAlive)
//      ;
}
