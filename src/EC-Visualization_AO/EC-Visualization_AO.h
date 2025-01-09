
#ifndef _aspect_enhanced_capitanio_visualization_h
#define _aspect_enhanced_capitanio_visualization_h

#include <aspect/postprocess/visualization.h>
#include <aspect/simulator_access.h>
#include <aspect/postprocess/interface.h>

#include <deal.II/numerics/data_postprocessor.h>

#include "../EC-Additional_outputs/EC-Additional_outputs.h"
#include "../EC-Material_model/EC-Material_model.h"

#include <vector>
#include <string>

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim> class EC_Visualization_AO : public DataPostprocessor<dim>, public SimulatorAccess<dim>, public Interface<dim>
      {
        public:
          std::vector<std::string> get_names () const override;
          std::vector<DataComponentInterpretation::DataComponentInterpretation> get_data_component_interpretation () const override;
          UpdateFlags get_needed_update_flags () const override;
          void evaluate_vector_field ( const DataPostprocessorInputs::Vector<dim> & input_data, std::vector<Vector<double>> & computed_quantities ) const override;
          static void declare_parameters ( ParameterHandler & prm );
          void parse_parameters ( ParameterHandler & prm ) override;
        private:
          std::vector<std::string> property_names;
      };
    }
  }
}

#endif
