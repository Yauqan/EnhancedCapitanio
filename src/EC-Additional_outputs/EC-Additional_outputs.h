
#ifndef _aspect_enhanced_capitanio_additional_outputs_AO_h
#define _aspect_enhanced_capitanio_additional_outputs_AO_h

#include <aspect/material_model/interface.h>
#include <aspect/material_model/rheology/drucker_prager.h>
#include <aspect/material_model/equation_of_state/linearized_incompressible.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/simulator_access.h>

#include <vector>

namespace aspect
{
  using namespace dealii;
  namespace MaterialModel
  {
      template <int dim> class EC_Additional_outputs : public AdditionalMaterialOutputs<dim>
      {
        public:
          EC_Additional_outputs ( const unsigned int N ) : liquidus(N), solidus(N), Ta(N), Tss(N), Finter(N), Tstar(N) {}
          std::vector<double> liquidus;
          std::vector<double> solidus;
          std::vector<double> Ta;
          std::vector<double> Tss;
          std::vector<double> Finter;
          std::vector<double> Tstar;
      };
  }
}

#endif
